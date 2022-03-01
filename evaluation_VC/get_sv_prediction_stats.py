import os,sys
import argparse

import gzip


def main(args):


    nr_true = 0
    for line in gzip.open(args.vcf_true,'r'):
        if line[0] != '#':
            nr_true += 1

    nr_pred = 0
    for line in open(args.vcf_pred, "r"):
        if line[0] != '#':
            nr_pred += 1

    nr_tp = 0
    for line in open(args.vcf_tp, "r"):
        if line[0] != '#':
            nr_tp += 1

    # print(nr_pred, nr_tp)
    precision = nr_tp / nr_pred
    recall = nr_tp / nr_true
    f_score = (2*precision*recall)/(precision + recall)

    # print(args.R, args.A, args.V)
    out_str = "{0},{1},{2},{3},{4},{5}".format(args.R, args.A, args.V, round(100*recall, 1), round(100*precision, 1), round(100*f_score, 1))
    # sys.stdout.write(out_str)
    print(out_str)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('vcf_true', type=str, default=False, help='VCF file')
    parser.add_argument('vcf_pred', type=str, default=False, help='VCF file')
    parser.add_argument('vcf_tp', type=str, default=False, help='VCF file')
    parser.add_argument('-A', type=str, default="", help='aligner')
    parser.add_argument('-R', type=str, default="", help='Read length')
    parser.add_argument('-V', type=str, default="", help='Variant type')
    # parser.add_argument('--outfile', type=str, default=None, help='Path to file.')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)