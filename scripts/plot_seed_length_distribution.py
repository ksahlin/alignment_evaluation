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



def plot_histogram_seeds(data, outfolder, name, bins = 50):

    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.7)
    x = [d[0] for d in data]
    y1 = [d[1] for d in data]
    y2 = [d[2] for d in data]
    y3 = [d[3] for d in data]

    ax2 = ax.twinx()
    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.2))
    # ax.set_ylim(0, max(x)+10000)

    color=['blue','green','firebrick']
    # width = 0.3
    # print(x)
    # x1 = [d[0] - 0.5 for d in data]
    # x2 = [d[0] + width for d in data]
    # print(x1)
    # print(x2)

    print(y2)
    # p1 = ax.bar(x-(width*3), dfDataTest['ERRh'], width=width, color=color[0], align='center', label=r'$\eta_{ER,h}$')
    # p2 = ax.bar(x-(width*2), dfDataTest['HRRh'], width=width, color=color[1], align='center', label=r'$\eta_{ER,h}$')
    p1, = ax.plot(x, y1,  color=color[0], label=r'Count')
    p2, = ax2.plot(x, y2, color=color[1], label=r'Unique')
    p3, = ax3.plot(x, y3, color=color[2], label=">= 10 locations")

    # p5 = ax2.bar(x+(width*1), dfDataTest['HRRc'], width=width, color=color[4], align='center', label=r'$\eta_{HR,h}, \eta_{th,h}$')


    
    # ax.set_xlim(10, 100)
    # ax2.set_xlim(0, 100)

    ax.set_xlabel('Seed length', fontsize=16)
    ax.set_ylabel("Seed count (log scale)", fontsize=16)
    ax2.set_ylabel("Unique", fontsize=16)
    ax3.set_ylabel(">= 10 locations", fontsize=16)

    ax.set_xticks([i for i in range(20, 110, 10)])
    ax.set_xticklabels([str(i) for i in range(20, 110, 10)])

    ax.set_yscale('log')
    ax2.set_ylim(0.65, 1.0)
    y2s = [0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    ax2.set_yticks(y2s)
    ax2.set_yticklabels([str(int(100*i)) + "%" for i in y2s])

    ax3.set_ylim(0.0, 0.003)
    y3s = [0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003]
    # y3s = [0.002, 0.004, 0.006, 0.008, 0.01]
    ax3.set_yticks(y3s)
    ax3.set_yticklabels([str(round(100*i, 2)) + "%" for i in y3s])

    ax.yaxis.label.set_color(p1.get_color())
    ax2.yaxis.label.set_color(p2.get_color())
    ax3.yaxis.label.set_color(p3.get_color())

    tkw = dict(size=4, width=1.5)
    ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
    ax2.tick_params(axis='y', colors=p2.get_color(), **tkw)
    ax3.tick_params(axis='y', colors=p3.get_color(), **tkw)
    ax.tick_params(axis='x', **tkw)

    lns = [p1,p2,p3]
    ax.legend(handles=lns, loc='lower center')



    plt.tight_layout()
    outfile = os.path.join(outfolder, "{0}.pdf".format(name))
    plt.savefig(outfile)
    plt.close()
    plt.cla()
    plt.clf()


def main(args):

    data = []
    for line in open(args.positions, "r"):
        # nohash,Sahlin1,contig-120_0,75,152
        seed_len, count, frac_unique, frac_10_or_more = line.strip().split(",")
        if int(count) > 1000: #int(seed_len) < 150:
            data.append( (int(seed_len), int(count), float(frac_unique), float(frac_10_or_more)) )

    plot_histogram_seeds(data, args.outfolder, "seed_stats", bins = len(data) )



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate sampling dispersity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('positions', type=str,  default=False, help='Input positions file')
    parser.add_argument('outfolder', type=str,  default=False, help='outfolder')


    args = parser.parse_args()
    mkdir_p(args.outfolder)

    main(args)


