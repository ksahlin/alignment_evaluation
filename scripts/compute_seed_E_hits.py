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

def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'U':'A', 'u':'a', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def get_minimizers(seq, k_size, w, seed_counts):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    window_kmers = deque([hash(seq[i:i+k_size]) for i in range(w)])
    curr_min = min(window_kmers)
    j = list(window_kmers).index(curr_min)
    # minimizers = [ seq[j:j+k_size] ]

    kmer = seq[j:j+k_size]
    kmer_rc = reverse_complement(kmer)
    if kmer < kmer_rc:
        seed_counts[kmer] += 1
    else:
        seed_counts[kmer_rc] += 1


    for i in range(w+1,len(seq) - k_size):
        new_kmer = hash(seq[i:i+k_size])
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = min(window_kmers)
            j = list(window_kmers).index(curr_min) + i - w 
            # minimizers.append( seq[j:j+k_size]  )
            # seed_counts[seq[j:j+k_size]] += 1
            kmer = seq[j:j+k_size]
            kmer_rc = reverse_complement(kmer)
            if kmer < kmer_rc:
                seed_counts[kmer] += 1
            else:
                seed_counts[kmer_rc] += 1

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            # minimizers.append( seq[i:i+k_size] )
            # seed_counts[seq[i:i+k_size]] += 1
            kmer = seq[i:i+k_size]
            kmer_rc = reverse_complement(kmer)
            if kmer < kmer_rc:
                seed_counts[kmer] += 1
            else:
                seed_counts[kmer_rc] += 1
    # return minimizers


def get_syncmers(seq, k, s, t, seed_counts):
    window_smers = deque([hash(seq[i:i+s]) for i in range(0, k - s + 1 )])
    curr_min = min(window_smers)
    pos_min =  window_smers.index(curr_min)
    syncmers = []
    if pos_min == t:
        kmer = seq[0 : k]
        kmer_rc = reverse_complement(kmer)
        if kmer < kmer_rc:
            seed_counts[kmer] += 1
        else:
            seed_counts[kmer_rc] += 1

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
            kmer_rc = reverse_complement(kmer)
            if kmer < kmer_rc:
                seed_counts[kmer] += 1
            else:
                seed_counts[kmer_rc] += 1

    # return syncmers

def print_stats(method, k, seed_counts):
    # seed_counts = defaultdict(int)
    # for minm_list in results:
    #     for m in minm_list:
    #         seed_counts[m] += 1

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
    frac_masked = 1 - total_seed_count_1000_lim/total_seed_count
    print("{0},{1},{2},{3},{4}".format(method, k, total_seed_count, int(round(total_seed_count_sq / total_seed_count,0)), round(100*frac_masked, 1) ))
    # print("{0},{1},{2},{3},{4}".format(method, k, total_seed_count_1000_lim, int(round(total_seed_count_sq_1000_lim / total_seed_count_1000_lim,0)), 1000))


def min_single_helper(arguments):
    return get_minimizers(*arguments)

def syncmers_single_helper(arguments):
    return get_syncmers(*arguments)

def main(args):

    genome = {acc: seq.upper() for (acc, (seq, _)) in readfq(open(args.fasta, 'r'))}
    n = 10000000
    for acc,seq in list(genome.items()):
        acc = acc.split()[0]
        # print(acc)
        genome[acc] = seq.replace("N", "") # remove Ns

    k = args.k
    if args.type == "minimizers":
        w = 9
        seed_counts = defaultdict(int)
        for acc, seq in genome.items():
            M = get_minimizers(seq, k, w, seed_counts)
        print_stats("minimizers", k, seed_counts)

    elif args.type == "syncmers":
        s = k-4 
        t = 2 # creates open syncmer with mid point with is used in strobealign
        seed_counts = defaultdict(int)
        for acc, seq in genome.items():
            get_syncmers(seq, k, s, t, seed_counts)
        print_stats("syncmers", k, seed_counts)       






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', type=str,  default=False, help='Path to genome')
    # parser.add_argument('n', type=int,  default=4, help='Nr cores')
    parser.add_argument('k', type=int,  default=20, help='k-mer size')
    parser.add_argument('--type', type=str,  default=False, help='Either syncmenrs or minimizers')

    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)