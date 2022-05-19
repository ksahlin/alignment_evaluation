import os,sys
import argparse

import random
from collections import defaultdict, deque


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


def minimizers(seq, k_size, w, seed_counts):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = min(window_kmers)

    seed_counts[curr_min] += 1
    # minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w+1,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = min(window_kmers)
            seed_counts[curr_min] += 1
            # minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            seed_counts[curr_min] += 1
            # minimizers.append( (curr_min, i) )

    # return minimizers


def syncmers(seq, k, s, t, seed_counts ):
    window_smers = deque([hash(seq[i:i+s]) for i in range(0, k - s + 1 )])
    curr_min = min(window_smers)
    pos_min =  window_smers.index(curr_min)
    syncmers = []
    if pos_min == t:
        kmer = seq[0 : k]
        # syncmers = [ (curr_min, 0) ]
        seed_counts[curr_min] += 1

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
            seed_counts[kmer] += 1
            # syncmers.append( (kmer,  i - (k - s) ) )

    # return syncmers

def print_stats(method, k, seed_counts):
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
    
    print("{0},{1},{2},{3},{4}".format(method, k, total_seed_count, int(round(total_seed_count_sq / total_seed_count,0)), -1 ))
    print("{0},{1},{2},{3},{4}".format(method, k, total_seed_count_1000_lim, int(round(total_seed_count_sq_1000_lim / total_seed_count_1000_lim,0)), 1000))


def main(args):

    genome = {acc: seq for (acc, (seq, _)) in readfq(open(args.fasta, 'r'))}

    for acc,seq in genome.items():
        acc = acc.split()[0]
        # print(acc)
        genome[acc] = seq.replace("N", "") # remove Ns

    print(len(genome), sum([len(v) for k,v in genome.items()]))


    w = 10
    for k in range(20, 101, 10):

        # compute minimizer stats
        seed_counts = defaultdict(int)
        for acc, seq in genome.items():
            minimizers(seq, k, w, seed_counts)
        print_stats("minimizers", k, seed_counts)


        # compute syncmer stats
        seed_counts = defaultdict(int)
        s = k-4
        t = 2 # creates open syncmer with mid point with is used in strobealign
        for acc, seq in genome.items():
            syncmers(seq, k, s, t, seed_counts )
        print_stats("syncmers", k, seed_counts)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', type=str,  default=False, help='Path to genome')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)