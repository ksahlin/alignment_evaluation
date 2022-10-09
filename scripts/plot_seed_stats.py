import os,sys
import argparse

import random

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot


def plot_e_hits(input_csv, outfolder, linewidth = 2.5):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.5)
    sns.set_style("whitegrid")
    indata = pd.read_csv(input_csv)
    g = sns.relplot(data=indata, x="median_seed_size", y="E_hits", hue="type", style="genome", linewidth = linewidth, kind="line", markers=True, markersize=8)
    g.set_axis_labels("Seed size", "E-hits")
    g.set(yscale="log")
    g.set( yticks= [i for i in range(10,99,10)] + [i for i in range(100,999,100)] + [i for i in range(1000,2999,1000)]) #, ylim=(0, 5200))

    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    # g.set(ylim=(94, 99), xticks=[50,75,100,150,200,250,300,500])
    # g.set_xticklabels(rotation=60, labels=[50,75,100,150,200,250,300,500])
    g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "e_hits.pdf"))
    plt.close()

def plot_frac_masked(input_csv, outfolder, linewidth = 2.5):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.5)
    sns.set_style("whitegrid")
    indata = pd.read_csv(input_csv)
    g = sns.relplot(data=indata, x="median_seed_size", y="fraction_masked_above_1000", hue="type", style="genome", linewidth = linewidth, kind="line", markers=True, markersize=8)
    g.set_axis_labels("Seed size", "Hard masked (%)")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    # g.set(ylim=(94, 99), xticks=[50,75,100,150,200,250,300,500])
    # g.set_xticklabels(rotation=60, labels=[50,75,100,150,200,250,300,500])
    g.tight_layout()
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "frac_masked.pdf"))
    plt.close()



def main(args):
    sns.set_style("whitegrid")

    plot_e_hits(args.csvfile, args.outfolder, linewidth = 2)
    plot_frac_masked(args.csvfile, args.outfolder, linewidth = 2)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('csvfile', type=str, help='results file')
    parser.add_argument('outfolder', type=str,  help='outfolder to plots.')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)