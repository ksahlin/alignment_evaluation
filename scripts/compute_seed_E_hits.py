import os,sys
import argparse

import random
from collections import defaultdict, deque

import signal
from multiprocessing import Pool

from time import time


'''
    Below awesome fast[a/q] reader function taken 
    from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def get_minimizers(seq, k_size, w):
    minimizers = []
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    window_kmers = deque([seq[i:i+k_size] for i in range(w)])
    # print(len(window_kmers))
    curr_min = min(window_kmers)
    # seed_counts[curr_min] += 1
    minimizers.append(curr_min)

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = min(window_kmers)
            # seed_counts[curr_min] += 1
            minimizers.append(curr_min)

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            # seed_counts[curr_min] += 1
            minimizers.append(curr_min)
    return minimizers


def get_syncmers(seq, k, s, t):
    window_smers = deque([hash(seq[i:i+s]) for i in range(0, k - s + 1 )])
    curr_min = min(window_smers)
    pos_min =  window_smers.index(curr_min)
    syncmers = []
    if pos_min == t:
        kmer = seq[0 : k]
        # seed_counts[kmer] += 1
        syncmers.append(kmer)

    for i in range(k - s + 1, len(seq) - s):
        new_smer = hash(seq[i:i+s])
        # updating window
        discarded_smer = window_smers.popleft()
        window_smers.append(new_smer)

        # Make this faster by storing pos of minimum
        curr_min = min(window_smers)
        pos_min = window_smers.index(curr_min)
        if pos_min == t:
            kmer = seq[i - (k - s) : i - (k - s) + k]
            # seed_counts[kmer] += 1
            syncmers.append(kmer)

    return syncmers

def print_stats(method, k, results):
    seed_counts = defaultdict(int)
    for minm_list in results:
        for m in minm_list:
            seed_counts[m] += 1

    total_seed_count_sq = 0
    total_seed_count = 0
    total_seed_count_sq_1000_lim = 0
    total_seed_count_1000_lim = 0
    for seed_id, cnt in seed_counts.items():
        total_seed_count += cnt
        total_seed_count_sq += cnt**2
        if cnt <= 1000:
            total_seed_count_1000_lim += cnt
            total_seed_count_sq_1000_lim += cnt**2
    fraq_masked = 1 - total_seed_count_1000_lim/total_seed_count
    print("{0},{1},{2},{3},{4}".format(method, k, total_seed_count, int(round(total_seed_count_sq / total_seed_count,0)), fraq_masked ))
    # print("{0},{1},{2},{3},{4}".format(method, k, total_seed_count_1000_lim, int(round(total_seed_count_sq_1000_lim / total_seed_count_1000_lim,0)), 1000))


def min_single_helper(arguments):
    return get_minimizers(*arguments)

def syncmers_single_helper(arguments):
    return get_syncmers(*arguments)

def main(args):

    genome = {acc: seq for (acc, (seq, _)) in readfq(open(args.fasta, 'r'))}

    for acc,seq in genome.items():
        acc = acc.split()[0]
        # print(acc)
        genome[acc] = seq.replace("N", "") # remove Ns

    print(len(genome), sum([len(v) for k,v in genome.items()]))

    tot_seeding = 0.0
    tot_dict_filling  = 0.0
    w = 10
    for k in range(20, 101, 10):

        # compute minimizer stats
        start_seed = time()
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=int(args.n))
        try:
            res = pool.map_async(min_single_helper, [ (seq, k, w) for acc, seq in genome.items()] )
            results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            pool.close()
        pool.join()
        stop_seed = time()
        tot_seeding += stop_seed - start_seed
        # print("finished multi")


        # for acc, seq in genome.items():
        #     get_minimizers(seq, k, w, seed_counts)
        # seed_counts = defaultdict(int)
        start_fill = time()
        print_stats("minimizers", k, results)
        stop_fill = time()
        tot_dict_filling += stop_fill - start_fill

        # compute syncmer stats
        start_seed = time()

        s = k-4
        t = 2 # creates open syncmer with mid point with is used in strobealign
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGINT, original_sigint_handler)
        pool = Pool(processes=int(args.n))
        try:
            res = pool.map_async(syncmers_single_helper, [ (seq, k, s, t) for acc, seq in genome.items()] )
            results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            pool.close()
        pool.join()
        stop_seed = time()
        tot_seeding += stop_seed - start_seed
        # print("finished multi")

        start_fill = time()
        print_stats("syncmers", k, results)
        stop_fill = time()
        tot_dict_filling += stop_fill - start_fill

        # print("Total seeding:", tot_seeding)
        # print("Total dict filling:", tot_dict_filling)

        # seed_counts = defaultdict(int)
        # s = k-4
        # t = 2 # creates open syncmer with mid point with is used in strobealign
        # for acc, seq in genome.items():
        #     syncmers(seq, k, s, t, seed_counts )




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', type=str,  default=False, help='Path to genome')
    parser.add_argument('n', type=int,  default=4, help='Nr cores')

    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)