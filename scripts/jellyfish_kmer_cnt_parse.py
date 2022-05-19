
import sys, os, argparse

def main(args):
	tot_unique_seeds = 0
	more_than_10 = 0
	tot_kmer_count = 0
	count_squared = 0
	for line in open(args.kmer_count_file, 'r'):
		tmp = line.split() 
		nr_occ = int(tmp[0])
		cnt = int(tmp[1])

		tot_unique_seeds += cnt
		tot_kmer_count += cnt*nr_occ
		count_squared += cnt*(nr_occ**2)

		if nr_occ == 1:
			unique = cnt
		elif nr_occ >= 10:
			more_than_10 += cnt

	frac_unique = unique / tot_unique_seeds
	frac_more_than_10 = more_than_10 / tot_unique_seeds
	e_count = count_squared / tot_kmer_count
	# print(tot_kmer_count)
	# print(e_count)
	print("{0},{1},{2},{3}".format(args.genome, args.k, tot_kmer_count, e_count))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('kmer_count_file', type=str, default=False, help='Kmer counts histogram')
    parser.add_argument("k", type=int, help='Kmer size')
    parser.add_argument("genome", type=str, help='genome (e.g., chm13 or hg38')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)