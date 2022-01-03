import os,sys
import argparse
import errno

import random
from collections import defaultdict, Counter
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")


def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise



def plot_histogram_distance(x, outfolder, h, l, name, bins=50):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(x.keys(), x.values())
    ax.set_ylabel('Count')
    ax.set_xlabel('Offset of next strobe sampled')
    # plt.title('HASH: {0}, LINK: {1}'.format(h,l))
    # plt.xlim(0, 50)
    # plt.yscale('log')
    ax.set_xticks([1,5,10,15,20])
    ax.set_xticklabels(['1','5','10','15','20'])
    plt.tight_layout()
    outfile = os.path.join(outfolder, "{0}_{1}_{2}.pdf".format(h,l, name))
    plt.savefig(outfile)
    plt.close()
    plt.cla()
    plt.clf()


def main(args):

    distances_sampled = []
    for line in open(args.positions, "r"):
        # nohash,Sahlin1,contig-120_0,75,152
        h, l, ref, p1, p2 = line.strip().split(",")
        d = int(p2) - int(p1)
        distances_sampled.append(d)

    C2 = Counter(distances_sampled)
    # most_repetitive = C2.most_common(1)
    # print(h, l, "most_repetitive distance:", most_repetitive[1])
    plot_histogram_distance(C2, args.outfolder, h, l, "skew_distribution", bins=len(C2))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate sampling dispersity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('positions', type=str,  default=False, help='Input positions file')
    parser.add_argument('outfolder', type=str,  default=False, help='outfolder')


    args = parser.parse_args()
    mkdir_p(args.outfolder)

    main(args)


