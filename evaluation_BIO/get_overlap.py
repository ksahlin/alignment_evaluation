


import os,sys
import argparse

import pysam


def read_sam(sam_file):
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    read_positions = {} # acc -> [ref_id, ref_start, refstop]


    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16: # single end
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

            if read.is_unmapped: 
                read_positions[q_name] = False
            else:
                read_positions[q_name] = (read.reference_name, read.reference_start, read.reference_end)
        
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


def overlap(t1, t2):
    q_id, q_a, q_b = t1
    p_id, p_a, p_b = t2
    assert q_a <= q_b and p_a <= p_b
    return (q_id == p_id) and  ( (p_a <= q_a <= p_b) or (p_a <= q_b <= p_b) or (q_a <= p_a <= q_b) or (q_a <= p_b <= q_b) )

def get_stats(sam1, sam2, sam3):
    total_aligned = {"sam1" : 0, "sam2": 0, "sam3": 0}
    overlaps = { "all":0, 
                 "sam1-sam2" : 0, "sam1-sam3" : 0, "sam2-sam3" : 0,
                 "sam1-unique" : 0, "sam2-unique" : 0, "sam3-unique" : 0}

    nr_total = len(sam1)
    print("TOT reads: ", nr_total)

    for read_acc in sam1:
        # check if mapped
        if sam1[read_acc]: # mapped in sam1
            total_aligned["sam1"] += 1
            sam1_rid, sam1_a, sam1_b = sam1[read_acc]
        if sam2[read_acc]: # mapped in sam2
            total_aligned["sam2"] += 1
            sam2_rid, sam2_a, sam2_b = sam2[read_acc]
        if sam3[read_acc]: # mapped in sam3
            total_aligned["sam3"] += 1
            sam3_rid, sam3_a, sam3_b = sam3[read_acc]

        # do the actual overlap analyses
        if sam1[read_acc] and sam2[read_acc] and sam3[read_acc]:
            o1 = overlap(sam1[read_acc], sam2[read_acc])
            o2 = overlap(sam1[read_acc], sam3[read_acc])
            o3 = overlap(sam2[read_acc], sam3[read_acc])

            if o1 and o2 and o3:



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
    if nr_aligned > 0:
        accuracy = 100*correct/nr_aligned
    return total_aligned, overlaps


def main(args):

    sam1 = read_sam(args.sam1)
    sam2 = read_sam(args.sam2)
    sam3 = read_sam(args.sam3)

    # if args.predicted_sam:
    #     predicted = read_sam(args.predicted_sam)
    # elif args.predicted_paf:
    #     predicted, mapped_to_multiple_pos = read_paf(args.predicted_paf)
    #     # print("Number of reads mapped to several positions (using first pos):", mapped_to_multiple_pos)

    total_aligned, overlaps = get_stats(sam1, sam2, sam3)

    for method, tot_aln in total_aligned.items():
        print("{0}: {1}".format(method, round(tot_aln/8000000,5)))

    for ovl, nr_in_common in overlaps.items():
        print("{0} : {1}".format(ovl, nr_in_common))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sam1', type=str, default=False, help='True coordinates (SAM)')
    parser.add_argument('--sam2', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    parser.add_argument('--sam3', type=str, default="", help='Perdicted coordinates (SAM/PAF)')

    # parser.add_argument('--paf1', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    # parser.add_argument('--paf2', type=str, default="", help='Perdicted coordinates (SAM/PAF)')
    # parser.add_argument('--paf3', type=str, default="", help='Perdicted coordinates (SAM/PAF)')

    parser.add_argument('--outfile', type=str, default=None, help='Path to file.')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)