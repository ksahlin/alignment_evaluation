

import os,sys
import argparse

import pysam
from collections import defaultdict


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

def get_stats(truth, predicted1, predicted2, out_misaligned, out_unaligned, logfile):

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
                    out_misaligned.write(">{0}\n{1}\n".format(read.query_name, read.query_sequence))
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
        logfile.write("Total uniquely misaligned on {0}: {1}\n".format(true_ref, s))

    out_misaligned.close()
    out_unaligned.close()
    # print("good (only method 2):", good_ok)
    # print("missed (only method 2):", missed_ok)
    # print("Unaligned (only method 2):", unaligned_ok) 


def get_stats_individual(truth, predicted, logfile, method):
    logfile.write("")
    logfile.write("METHOD: {0}".format(method))
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
            pred_ref_id, pred_start, pred_stop, read_p1 = predicted1[read_acc]
            if pred_ref_id == true_ref_id and overlap(pred_start, pred_stop, true_start, true_stop):
                good_method += 1
            else:
                bad_method += 1
                misaligned_dict[true_ref_id][pred_ref_id] += 1


    logfile.write("nr_aligned:{0}\n".format(nr_aligned_method1))
    logfile.write("good_method:{0}\n".format(good_method1))
    logfile.write("bad_method:{0}\n".format(bad_method1))
    logfile.write("unaligned_method:{0}\n".format(unaligned_method1))
    logfile.write("to_improve:{0}\n".format(to_improve))

    logfile.write("TRUE REF, PRED REF, MISALIGNED")
    for true_ref in misaligned_dict:
        s = 0
        for pred_ref in misaligned_dict[true_ref]:
            s+= misaligned_dict[true_ref][pred_ref]
            # logfile.write("{0},{1}: {2}\n".format(true_ref, pred_ref, misaligned_dict[true_ref][pred_ref]))
        logfile.write("Total uniquely misaligned on {0}: {1}\n".format(true_ref, s))


def main(args):

    truth = read_sam(args.truth, args.n)

    if args.predicted_sam_method1:
        predicted1 = read_sam(args.predicted_sam_method1, args.n)

    if args.predicted_sam_method2:
        predicted2 = read_sam(args.predicted_sam_method2, args.n)

    get_stats_individual(truth, predicted1, "minimap2")
    get_stats_individual(truth, predicted2, "strobealign")

    # get_stats(truth, predicted1, predicted2, open(args.om, "w"), open(args.ou, "w"), open(args.logfile, "w"))





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--truth', type=str, default=False, help='True coordinates (SAM)')
    parser.add_argument('--predicted_sam_method1', type=str, default="", help='Predicted coordinates (SAM/PAF)')
    parser.add_argument('--predicted_sam_method2', type=str, default="", help='Predicted coordinates (SAM/PAF)')
    # parser.add_argument('--predicted_paf', type=str, default="", help='Predicted coordinates (SAM/PAF)')
    parser.add_argument('--om', type=str, default=None, help='Path to outfile misaligned file.')
    parser.add_argument('--ou', type=str, default=None, help='Path to outfile unaligned file.')
    parser.add_argument('--n', type=int, default=1000000, help='Number of reads to infer accuracy from (default is first 1M reads).')
    parser.add_argument('--logfile', type=str, default= "log.txt", help='Path to logfile with results.')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)