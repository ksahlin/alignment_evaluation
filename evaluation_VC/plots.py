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


def plot_sv_calling_results(input_csv, outfolder):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.6)
    sns.set_style("whitegrid")
    indata = pd.read_csv(input_csv)

    # Recall
    g = sns.relplot(data=indata, x="dataset", y="recall", hue="tool", style="type", kind="line", #dashes = dashes,
        col="type",  # hue="datastructure", style="datastructure",
        col_wrap=3, col_order=["SNV", "INDEL"])
    g.set_axis_labels("Dataset", "Recall")
    g.set(ylim=(94, 100), xticks=[100,150,200,250,300])
    plt.savefig(os.path.join(outfolder, "recall_plot.pdf"))
    plt.clf()

    # Precision
    g = sns.relplot(data=indata, x="dataset", y="precision", hue="tool", style="type", kind="line", #dashes = dashes,
    col="type",  # hue="datastructure", style="datastructure",
    col_wrap=3, col_order=["SNV", "INDEL"])
    g.set_axis_labels("Dataset", "Recall")
    g.set(ylim=(94, 100), xticks=[100,150,200,250,300])
    plt.savefig(os.path.join(outfolder, "precision_plot.pdf"))
    plt.clf()

    # F-Score
    g = sns.relplot(data=indata, x="dataset", y="fscore", hue="tool", style="type", kind="line", #dashes = dashes,
    col="type",  # hue="datastructure", style="datastructure",
    col_wrap=3, col_order=["SNV", "INDEL"])
    g.set_axis_labels("Dataset", "F-score")
    g.set(ylim=(94, 100), xticks=[100,150,200,250,300])
    plt.savefig(os.path.join(outfolder, "fscore_plot.pdf"))

    plt.close()

def plot_percentage_aligned(input_csv, outfolder):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.6)
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)

    # dashes = { "strobemap" : "", 
    #             "spaced_sparse": (5,5),
    #             "spaced_dense": (5,5),
    #              "minstrobes2" : (1,1),
    #              "minstrobes3" : (1,1),
    #              "randstrobes2" : (1,1),
    #              "randstrobes3" : (1,1),
    #              "hybridstrobes2" : (1,1),
    #              "hybridstrobes3" : (1,1)}
    # print(indata)
    g = sns.relplot(
        data=indata, x="read_length", y="aligned", hue="tool", style="type",
        col="dataset", kind="line",  #dashes = dashes, hue="datastructure", style="datastructure",
        col_wrap=3, col_order=["SIM1", "SIM2", "SIM3"])
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_axis_labels("Read length", "Percentage aligned")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(ylim=(94, 100), xticks=[100,150,200,250,300])
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "percentage_aligned_plot.eps"))
    plt.savefig(os.path.join(outfolder, "percentage_aligned_plot.pdf"))
    plt.close()

def plot_overaligned(input_csv, outfolder):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.6)
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)

    # dashes = { "strobemap" : "", 
    #             "spaced_sparse": (5,5),
    #             "spaced_dense": (5,5),
    #              "minstrobes2" : (1,1),
    #              "minstrobes3" : (1,1),
    #              "randstrobes2" : (1,1),
    #              "randstrobes3" : (1,1),
    #              "hybridstrobes2" : (1,1),
    #              "hybridstrobes3" : (1,1)}
    # print(indata)
    g = sns.relplot(
        data=indata, x="read_length", y="overaligned", hue="tool", style="type",
        col="dataset", kind="line",  #dashes = dashes, hue="datastructure", style="datastructure",
        col_wrap=3, col_order=["SIM1", "SIM2", "SIM3"])
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_axis_labels("Read length", "Overaligned")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set( xticks=[100,150,200,250,300]) #ylim=(40, 100),
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "overaligned_plot.eps"))
    plt.savefig(os.path.join(outfolder, "overaligned_plot.pdf"))
    plt.close()


def plot_memory_usage(input_csv, outfolder):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.6)
    # tool,dataset,read_length,time,memory
    indata = pd.read_csv(input_csv)

    sns.set_style("whitegrid")

    g = sns.relplot(
        data=indata, x="read_length", y="memory", hue="tool", style="type",
        col="dataset", kind="line",  #dashes = dashes, hue="datastructure", style="datastructure",
        col_wrap=3, col_order=["SIM1", "SIM2", "SIM3"])
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_axis_labels("Read length", "Memory usage (Gb)")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set( xticks=[100,150,200,250,300]) #ylim=(40, 100),
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "memory_plot.eps"))
    plt.savefig(os.path.join(outfolder, "memory_plot.pdf"))
    plt.close()

def plot_runtime(input_csv, outfolder):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.6)
    # tool,dataset,read_length,time,memory
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)
    g = sns.relplot(
        data=indata, x="read_length", y="time", hue="tool", style="type",
        col="dataset", kind="line",  #dashes = dashes, hue="datastructure", style="datastructure",
        col_wrap=3, col_order=["SIM1", "SIM2", "SIM3"])
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_axis_labels("Read length", "Time (sec)")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set( xticks=[100,150,200,250,300]) #ylim=(40, 100),
    g.set( yticks=[0,1000,2000,4000,6000,12000,18000,24000]) #ylim=(40, 100),

    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "time_plot.eps"))
    plt.savefig(os.path.join(outfolder, "time_plot.pdf"))
    plt.close()


def main(args):
    sns.set_style("whitegrid")

    plot_sv_calling_results(sv_csv, args.outfolder)
    plot_runtime(runtime_mem_csv, args.outfolder)
    plot_memory_usage(runtime_mem_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('sv_csv', type=str, help='results file')
    parser.add_argument('runtime_mem_csv', type=str, help='results file')
    parser.add_argument('outfolder', type=str,  help='outfolder to plots.')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)