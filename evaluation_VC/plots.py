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


def plot_sv_calling_results(input_csv, outfolder, palette, tools):
    sns.set(rc={'figure.figsize':(12,4)})
    matplotlib.rcParams.update({'font.size': 14})
    sns.set(font_scale=1.4)
    sns.set_style("whitegrid")
    # pyplot.figure(figsize=(4,16))
    indata = pd.read_csv(input_csv)
    snvs = indata[indata["type"] == "SNV"]
    indels = indata[indata["type"] == "INDEL"]

    f, axes = plt.subplots(1, 3)
    ax1 = sns.lineplot(data=snvs, x="dataset", y="recall", hue="tool", hue_order = tools, ax=axes[0], palette=palette, linewidth = 2.0)
    ax2 = sns.lineplot(data=snvs, x="dataset", y="precision", hue="tool", hue_order = tools, ax=axes[1], palette=palette, linewidth = 2.0)
    ax3 = sns.lineplot(data=snvs, x="dataset", y="fscore", hue="tool", hue_order = tools, ax=axes[2], palette=palette, linewidth = 2.0)
    
    # Shrink current axis by 20%
    # box = ax3.get_position()
    # ax3.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    # ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax1.get_legend().remove()
    ax2.get_legend().remove()
    ax3.get_legend().remove()

    # ax1.set_ylim(75, 100) # for BIO
    # ax2.set_ylim(75, 100) # for BIO
    # ax3.set_ylim(75, 100) # for BIO
    ax1.set_ylim(90, 100) # for SIM
    ax2.set_ylim(90, 100) # for SIM
    ax3.set_ylim(90, 100) # for SIM

    ax1.set_ylabel("Recall (%)")
    ax2.set_ylabel("Precision (%)")
    ax3.set_ylabel("F-score (%)")
    # ax1.tick_params(labelrotation=45)
    # ax2.tick_params(labelrotation=45)
    # ax3.tick_params(labelrotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder, "sv_snv.pdf"))
    plt.clf()

    f, axes = plt.subplots(1, 3)
    ax1 = sns.lineplot(data=indels, x="dataset", y="recall", hue="tool", hue_order = tools, ax=axes[0], palette=palette, linewidth = 2.0)
    ax2 = sns.lineplot(data=indels, x="dataset", y="precision", hue="tool", hue_order = tools, ax=axes[1], palette=palette, linewidth = 2.0)
    ax3 = sns.lineplot(data=indels, x="dataset", y="fscore", hue="tool", hue_order = tools, ax=axes[2], palette=palette, linewidth = 2.0)

    # Shrink current axis by 20%
    # box = ax3.get_position()
    # ax3.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    # ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax1.get_legend().remove()
    ax2.get_legend().remove()
    ax3.get_legend().remove()
    # ax1.set_ylim(2, 8.5) # for BIO
    # ax2.set_ylim(2, 8.5) # for BIO
    # ax3.set_ylim(2, 8.5) # for BIO
    ax1.set_ylim(40, 60) # for SIM
    ax2.set_ylim(40, 60) # for SIM
    ax3.set_ylim(40, 60) # for SIM
    ax1.set_ylabel("Recall (%)")
    ax2.set_ylabel("Precision (%)")
    ax3.set_ylabel("F-score (%)")
    # ax1.tick_params(labelrotation=45)
    # ax2.tick_params(labelrotation=45)
    # ax3.tick_params(labelrotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder, "sv_indel.pdf"))
    plt.clf()

    # # Recall
    # g = sns.relplot(data=indata, x="dataset", y="recall", hue="tool", kind="line", #dashes = dashes, style="type",
    #     col="type", col_order=["SNV", "INDEL"]) # hue="datastructure", style="datastructure",  col_wrap=3, )
    # g.set_axis_labels("Dataset", "Recall")
    # g.set(ylim=(0, 100))
    # plt.savefig(os.path.join(outfolder, "recall_plot.pdf"))
    # plt.clf()

    # # Precision
    # g = sns.relplot(data=indata, x="dataset", y="precision", hue="tool", kind="line", #dashes = dashes, style="type",
    #     col="type", col_order=["SNV", "INDEL"]) # hue="datastructure", style="datastructure",  col_wrap=3, )
    # g.set_axis_labels("Dataset", "Recall")
    # g.set(ylim=(0, 100))
    # plt.savefig(os.path.join(outfolder, "precision_plot.pdf"))
    # plt.clf()

    # # F-Score
    # g = sns.relplot(data=indata, x="dataset", y="fscore", hue="tool", kind="line", #dashes = dashes, style="type",
    #     col="type", col_order=["SNV", "INDEL"]) # hue="datastructure", style="datastructure",  col_wrap=3, )
    # g.set_axis_labels("Dataset", "F-score")
    # g.set(ylim=(0, 100))
    # plt.savefig(os.path.join(outfolder, "fscore_plot.pdf"))

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
    palette = {
    'minimap2': 'tab:blue',
    'strobealign': 'tab:green',
    'bwa_mem': 'tab:orange',
    'accelalign': 'tab:red',
    'bowtie2' : 'tab:purple',
    'urmap' : 'tab:grey',
    'snap' : 'pink'
    }
    tools =["minimap2", "bwa_mem", 'accelalign', "bowtie2", "snap", "urmap", "strobealign"]

    plot_sv_calling_results(args.sv_csv, args.outfolder, palette, tools)
    # plot_runtime(runtime_mem_csv, args.outfolder)
    # plot_memory_usage(runtime_mem_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('sv_csv', type=str, help='results file')
    # parser.add_argument('runtime_mem_csv', type=str, help='results file')
    parser.add_argument('outfolder', type=str,  help='outfolder to plots.')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)