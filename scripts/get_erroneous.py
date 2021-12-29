

import os,sys
import argparse

import pysam
from collections import defaultdict


'''
    Below readfq function taken from https://github.com/lh3/readfq/blob/master/readfq.py
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
        name, seqs, last = last[1:], [], None
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


def read_sam(sam_file, n):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    read_positions = {} # acc -> [ref_id, ref_start, refstop]


    for i, read in enumerate(SAM_file.fetch(until_eof=True)):
        if i > n:
            break
        if read.flag == 0 or read.flag == 16: # single end
            # print(read.query_name, len(read_positions))
            read_positions[read.query_name] = (read.reference_name, read.reference_start, read.reference_end, read)
        
        elif read.is_paired:
            if read.is_read1:
        # elif read.flag == 99 or  read.flag == 83: # Paired end first
                # if not (read.flag == 99 or  read.flag == 83):
                #     print(read.query_name, read.flag)
                if "/1" not in read.query_name[-4:]:
                    q_name = read.query_name + "/1"
                else:
                    q_name = read.query_name

        # elif read.flag == 147 or read.flag == 163: # Paired end second
            if read.is_read2:
                # if not (read.flag == 147 or read.flag == 163):
                #     print(read.query_name, read.flag)
                if "/2" not in read.query_name[-4:]:
                    q_name = read.query_name + "/2"
                else:
                    q_name = read.query_name

            if read.is_unmapped: 
                read_positions[q_name] = False
            else:
                read_positions[q_name] = (read.reference_name, read.reference_start, read.reference_end, read)
        
        elif (not read.is_paired) and read.is_unmapped: # single and unmapped
            read_positions[read.query_name] = False


    return read_positions     


def read_paf(paf_file):
    read_positions = {} # acc -> [ref_id, ref_start, refstop]
    mapped_to_multiple_pos = 0
    for line in open(paf_file, 'r'):
        vals = line.split()
        read_acc, ref_name, reference_start, reference_end  = vals[0], vals[5], int(vals[7]), int(vals[8])

        if read_acc in read_positions:
            mapped_to_multiple_pos += 1
            continue
        else:
            read_positions[read_acc] = (ref_name, reference_start, reference_end)
    return read_positions, mapped_to_multiple_pos


def overlap(q_a, q_b, p_a, p_b):
    assert q_a <= q_b and p_a <= p_b
    return  (p_a <= q_a <= p_b) or (p_a <= q_b <= p_b) or (q_a <= p_a <= q_b) or (q_a <= p_b <= q_b)

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def get_stats(truth, predicted1, predicted2, out_unaligned, logfile):

    nr_aligned_method1 = len(predicted1)
    nr_aligned_method2 = len(predicted2)
    nr_total = len(truth)
    
    good_method1 = 0
    bad_method1 = 0
    unaligned_method1 = 0

    good_method2 = 0
    bad_method2 = 0
    unaligned_method2 = 0

    to_improve = 0
    misaligned_dict = defaultdict(lambda: defaultdict(int))
    misaligned_read_acc = set()
    # missed = 0
    # unaligned = 0
    # good_ok = 0
    # missed_ok = 0
    # unaligned_ok = 0

    for read_acc in truth:
        if not truth[read_acc]:
            continue
        if (read_acc not in predicted1) or (read_acc not in predicted2):
            # excluding from analysis since subsampling didnt sample them in at least on of the methods
            continue 
        true_ref_id, true_start, true_stop, read = truth[read_acc]
        if not predicted1[read_acc]:
            unaligned_method1 += 1
        else:
            pred_ref_id, pred_start, pred_stop, read_p1 = predicted1[read_acc]
            if pred_ref_id == true_ref_id and overlap(pred_start, pred_stop, true_start, true_stop):
                good_method1 += 1
            else:
                bad_method1 += 1

            if not predicted2[read_acc]:
                unaligned_method2 += 1
                # print(read_acc, "UNALIGNED STROBEALIGN", true_ref_id, true_start, true_stop )
                out_unaligned.write(">{0}_{2}\n{1}\n".format(read.query_name, read.query_sequence, read.cigarstring))
                continue

            else:
                pred2_ref_id, pred2_start, pred2_stop, read_p2 = predicted2[read_acc]
                if not (pred2_start < pred2_stop):
                    print(read_acc, read_p2.reference_name, read_p2.reference_start, read_p2.reference_end,read_p2.cigarstring, true_ref_id)
                
                if (pred2_ref_id == true_ref_id) and overlap(pred2_start, pred2_stop, true_start, true_stop):
                    pass
                else:
                    to_improve +=1
                    # print(read_acc, "MISALIGNED STROBEALIGN" , pred_ref_id, pred2_start, pred2_stop, true_ref_id, true_start, true_stop )
                    if pred_ref_id == true_ref_id and overlap(pred_start, pred_stop, true_start, true_stop):
                        misaligned_read_acc.add(read_acc[:-2])
                        misaligned_dict[true_ref_id][pred_ref_id] += 1


        if not predicted2[read_acc]:
            unaligned_method2 += 1
        else:
            pred_ref_id, pred_start, pred_stop, read_p2 = predicted2[read_acc]
            if pred_ref_id == true_ref_id and overlap(pred_start, pred_stop, true_start, true_stop):
                good_method2 += 1
            else:
                bad_method2 += 1


    logfile.write("nr_aligned method1:{0}\n".format(nr_aligned_method1))
    logfile.write("good_method1:{0}\n".format(good_method1))
    logfile.write("bad_method1:{0}\n".format(bad_method1))
    logfile.write("unaligned_method1:{0}\n".format(unaligned_method1))
    
    logfile.write("nr_aligned method2:{0}\n".format(nr_aligned_method2))
    logfile.write("good_method2:{0}\n".format(good_method2))
    logfile.write("bad_method2:{0}\n".format(bad_method2))
    logfile.write("unaligned_method2:{0}\n".format(unaligned_method2))

    logfile.write("to_improve:{0}\n".format(to_improve))

    logfile.write("TRUE REF, PRED REF, ONLY STROBEALIGN MISALIGNED")
    for true_ref in misaligned_dict:
        s = 0
        for pred_ref in misaligned_dict[true_ref]:
            s+= misaligned_dict[true_ref][pred_ref]
            logfile.write("{0},{1}: {2}\n".format(true_ref, pred_ref, misaligned_dict[true_ref][pred_ref]))
        logfile.write("Only misaligned by strobealign {0}: {1}\n".format(true_ref, s))

    out_unaligned.close()
    return misaligned_read_acc
    # print("good (only method 2):", good_ok)
    # print("missed (only method 2):", missed_ok)
    # print("Unaligned (only method 2):", unaligned_ok) 


def get_stats_individual(truth, predicted, logfile, method):
    logfile.write("\n\n")
    logfile.write("METHOD: {0}\n".format(method))
    nr_aligned_method = len(predicted)
    nr_total = len(truth)
    
    good_method = 0
    bad_method = 0
    unaligned_method = 0
    to_improve = 0
    misaligned_dict = defaultdict(lambda: defaultdict(int))

    for read_acc in truth:
        if not truth[read_acc]:
            continue
        if (read_acc not in predicted):
            # excluding from analysis since subsampling didnt sample them in at least on of the methods
            continue 
        true_ref_id, true_start, true_stop, read = truth[read_acc]
        if not predicted[read_acc]:
            unaligned_method += 1
        else:
            pred_ref_id, pred_start, pred_stop, read_p1 = predicted[read_acc]
            if pred_ref_id == true_ref_id and overlap(pred_start, pred_stop, true_start, true_stop):
                good_method += 1
            else:
                bad_method += 1
                misaligned_dict[true_ref_id][pred_ref_id] += 1


    logfile.write("nr_aligned:{0}\n".format(nr_aligned_method))
    logfile.write("good_method:{0}\n".format(good_method))
    logfile.write("bad_method:{0}\n".format(bad_method))
    logfile.write("unaligned_method:{0}\n".format(unaligned_method))
    logfile.write("to_improve:{0}\n".format(to_improve))

    logfile.write("TRUE REF, PRED REF, MISALIGNED\n")
    tot_mis = {}
    for true_ref in misaligned_dict:
        s = 0
        for pred_ref in misaligned_dict[true_ref]:
            s+= misaligned_dict[true_ref][pred_ref]
            logfile.write("{0},{1}: {2}\n".format(true_ref, pred_ref, misaligned_dict[true_ref][pred_ref]))
        # logfile.write("Total uniquely misaligned on {0}: {1}\n".format(true_ref, s))
        tot_mis[true_ref] = s

    for ref, nr in sorted(tot_mis.items(), key = lambda x: x[1], reverse=True):
        logfile.write("Total uniquely misaligned, originally from {0}: {1}\n".format(ref, nr))

    return tot_mis


def main(args):

    truth = read_sam(args.truth, args.n)

    if args.predicted_sam_method1:
        predicted1 = read_sam(args.predicted_sam_method1, args.n)

    if args.predicted_sam_method2:
        predicted2 = read_sam(args.predicted_sam_method2, args.n)
    f = open(args.logfile, "w")
    tot_mis1 = get_stats_individual(truth, predicted1, f, "minimap2")
    tot_mis2 = get_stats_individual(truth, predicted2, f, "strobealign")

    delta_missed = {}
    for ref1, nr1 in sorted(tot_mis1.items(), key = lambda x: x[1], reverse=True):
        if ref1 in tot_mis2:
            nr2 = tot_mis2[ref1]
            delta_missed[ref1] = nr2 - nr1

    f.write("\n")
    for ref, nr in sorted(delta_missed.items(), key = lambda x: x[1], reverse=True):
        f.write("More misaligned strobealign, originally from {0}: {1}\n".format(ref, nr))


    all_misaligned = get_stats(truth, predicted1, predicted2, open(args.ou, "w"), f)

    om1 = open(args.om + "_1.fq", "w")
    om2 = open(args.om + "_2.fq", "w")
    all_r1 = { acc: (seq, qual) for acc, (seq, qual) in readfq(open(args.fq1, 'r'))}
    all_r2 = {acc: (seq, qual) for acc, (seq, qual) in readfq(open(args.fq2, 'r'))}
    for read_acc in all_misaligned:
        seq, qual = all_r1[read_acc + "/1"]
        om1.write(">{0}\n{1}\n+\n{2}\n".format(read_acc +'/1', seq, qual))
        seq, qual = all_r2[read_acc + "/2"]
        om2.write(">{0}\n{1}\n+\n{2}\n".format(read_acc +'/2', seq, qual))

    om1.close()
    om2.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--truth', type=str, default=False, help='True coordinates (SAM)')
    parser.add_argument('--predicted_sam_method1', type=str, default="", help='Predicted coordinates (SAM/PAF)')
    parser.add_argument('--predicted_sam_method2', type=str, default="", help='Predicted coordinates (SAM/PAF)')
    parser.add_argument('--fq1', type=str, default="", help='Predicted coordinates (SAM/PAF)')
    parser.add_argument('--fq2', type=str, default="", help='Predicted coordinates (SAM/PAF)')
    parser.add_argument('--om', type=str, default=None, help='Prefix path name to outfile misaligned file.')
    parser.add_argument('--ou', type=str, default=None, help='Path to outfile unaligned file.')
    parser.add_argument('--n', type=int, default=1000000, help='Number of reads to infer accuracy from (default is first 1M reads).')
    parser.add_argument('--logfile', type=str, default= "log.txt", help='Path to logfile with results.')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)