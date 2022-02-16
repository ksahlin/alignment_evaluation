


import os,sys
import argparse

import pysam

import re

def parse_gnu_time(stderr_file):
    lines = open(stderr_file, 'r').readlines()
    index_time_mm2 = False
    index_time_aa = False
    index_time_strobemap = False

    for l in lines:
        wct_match = re.search('[\d.]+ real', l) 
        mem_match = re.search('[\d]+  maximum resident set size', l) 

        # minimap2
        index_time_mm2_match = re.search('\[M::main::[\d.:]+\*[\d.:]+\] loaded/built the index for', l) 
        
        # strobemap

        index_time_strobemap_match = re.search('Total time indexing: [\d.:]+', l) 

        # accelalign
        index_time_accelalign_match = re.search('Setup reference in [\d.:]+ secs', l) 

        if wct_match:
            wallclocktime = float(wct_match.group().split()[0])
        if mem_match:
            mem_tmp = int(mem_match.group().split()[0])
            memory_gb = mem_tmp / 1000000.0 

        if index_time_mm2_match:
            prefix_cut = index_time_mm2_match.group().split('main::')[1]
            final_time = prefix_cut.split('*')[0]
            index_time_mm2 = float(final_time.strip())

        if index_time_strobemap_match:
            # print(index_time_strobemap_match)
            index_time_strobemap = float(index_time_strobemap_match.group().split(':')[1].strip())
        
        if index_time_accelalign_match:
            prefix_cut = index_time_accelalign_match.group().split('reference in ')[1]
            final_time = prefix_cut.split(' secs')[0]
            index_time_aa = float(final_time.strip())

    if index_time_mm2:
        tot_wallclock_secs = wallclocktime - index_time_mm2
    elif index_time_aa:
        tot_wallclock_secs = wallclocktime - index_time_aa
    elif index_time_strobemap:
        tot_wallclock_secs = wallclocktime - index_time_strobemap
    else:
        tot_wallclock_secs = wallclocktime


    return round(tot_wallclock_secs,2), round(memory_gb,2)

    # vals = list(map(lambda x: float(x), wallclocktime.split(".") ))
    # if len(vals) == 3:
    #     h,m,s = vals
    #     tot_wallclock_secs = h*3600.0 + m*60.0 + s
    #     if index_time_mm2:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_mm2
    #     elif index_time_aa:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_aa
    #     elif index_time_strobemap:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_strobemap

    # elif len(vals) == 2:
    #     s, = vals
    #     tot_wallclock_secs = m*60.0 + s

    #     if index_time_mm2:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_mm2
    #     elif index_time_strobemap:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_strobemap

        # return vals, memory_gb

    # else:
    #     return "BUG","BUG"


def read_sam(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    read_positions = {} # acc -> [ref_id, ref_start, refstop]


    for read in SAM_file.fetch(until_eof=True):
        if (read.flag == 0 or read.flag == 16) and (not read.is_secondary): # single end
            # print(read.query_name, len(read_positions))
            read_positions[read.query_name] = (read.reference_name, read.reference_start, read.reference_end)
        
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

            if not read.is_secondary and read.is_unmapped: 
                read_positions[q_name] = False
            elif not read.is_secondary:
                read_positions[q_name] = (read.reference_name, read.reference_start, read.reference_end)
        
        elif (not read.is_paired) and read.is_unmapped and (not read.is_secondary): # single and unmapped
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
    # if (q_a == q_b) or (p_a == p_b):
    #     print("Cigar bug")
    return  (p_a <= q_a <= p_b) or (p_a <= q_b <= p_b) or (q_a <= p_a <= q_b) or (q_a <= p_b <= q_b)

def get_stats(truth, predicted):

    # nr_aligned = len(predicted)
    nr_total = len(truth)
    unaligned = 0 
    nr_aligned = 0
    over_mapped = 0
    correct = 0
    # print("MISSING:", set(truth).difference(set(predicted)))
    for read_acc in predicted:
        if not truth[read_acc]:
            over_mapped += 1
            continue
        if not predicted[read_acc]:
            unaligned += 1
            continue

        nr_aligned += 1

        pred_ref_id, pred_start, pred_stop = predicted[read_acc]
        true_ref_id, true_start, true_stop = truth[read_acc]
        # print(read_acc, pred_start, pred_stop, true_start, true_stop)
        if pred_ref_id == true_ref_id and overlap(pred_start, pred_stop, true_start, true_stop):
            correct += 1
            # print(read_acc)
        else:
            pass
    #         print(read_acc, pred_ref_id, pred_start, pred_stop, true_ref_id, true_start, true_stop )
    # print(correct)
    # print(nr_aligned)
    aligned_percentage = 100*(nr_aligned/nr_total)
    accuracy = 0.0
    # if nr_aligned > 0:
    #     accuracy = 100*correct/nr_aligned
    accuracy = 100*correct/nr_total
    return aligned_percentage, accuracy, over_mapped


def main(args):

    truth = read_sam(args.truth)

    if args.predicted_sam:
        predicted = read_sam(args.predicted_sam)
    elif args.predicted_paf:
        predicted, mapped_to_multiple_pos = read_paf(args.predicted_paf)
        # print("Number of reads mapped to several positions (using first pos):", mapped_to_multiple_pos)

    percent_aligned, percent_correct, over_mapped = get_stats(truth, predicted)
    tot_wallclock_secs, memory_gb = parse_gnu_time(args.time_mem)
    out_str = "{0},{1},{2},{3},{4}".format(round(percent_aligned, 5), round(percent_correct, 5), over_mapped, tot_wallclock_secs, memory_gb)
    print(out_str)
    # print("Percentage aligned: {0}".format(round(percent_aligned, 3)))
    # print("Accuracy: {0}".format(round(percent_correct, 3)))
    # print("Over-aligned (unmapped in grough truth file): {0}".format(over_mapped))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--truth', type=str, default=False, help='True coordinates (SAM)')
    parser.add_argument('--predicted_sam', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    parser.add_argument('--predicted_paf', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    parser.add_argument('--time_mem', type=str, default="", help='time mem file)')
    # parser.add_argument('--outfile', type=str, default=None, help='Path to file.')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)