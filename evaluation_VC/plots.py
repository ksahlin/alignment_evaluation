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


def plot_sv_calling_results(input_csv, outfolder, palette, tools, exp_type):
    sns.set(rc={'figure.figsize':(12,4)})
    matplotlib.rcParams.update({'font.size': 14})
    sns.set(font_scale=1.4)
    sns.set_style("whitegrid")
    # pyplot.figure(figsize=(4,16))
    indata = pd.read_csv(input_csv)
    snvs = indata[indata["variant"] == "SNV"]
    indels = indata[indata["variant"] == "INDEL"]

    f, axes = plt.subplots(1, 3)
    ax1 = sns.lineplot(data=snvs, x="dataset", y="recall", hue="tool", hue_order = tools, ax=axes[0], palette=palette, linewidth = 2.0)
    ax2 = sns.lineplot(data=snvs, x="dataset", y="precision", hue="tool", hue_order = tools, ax=axes[1], palette=palette, linewidth = 2.0)
    ax3 = sns.lineplot(data=snvs, x="dataset", y="f_score", hue="tool", hue_order = tools, ax=axes[2], palette=palette, linewidth = 2.0)
    
    # Shrink current axis by 20%
    # box = ax3.get_position()
    # ax3.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    # ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax1.get_legend().remove()
    ax2.get_legend().remove()
    ax3.get_legend().remove()

    if exp_type == 'bio':
        ax1.set_ylim(75, 100) # for BIO
        ax2.set_ylim(75, 100) # for BIO
        ax3.set_ylim(75, 100) # for BIO
    else:
        ax1.set_ylim(92, 100) # for SIM
        ax2.set_ylim(92, 100) # for SIM
        ax3.set_ylim(92, 100) # for SIM

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
    ax3 = sns.lineplot(data=indels, x="dataset", y="f_score", hue="tool", hue_order = tools, ax=axes[2], palette=palette, linewidth = 2.0)

    # Shrink current axis by 20%
    # box = ax3.get_position()
    # ax3.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    # ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax1.get_legend().remove()
    ax2.get_legend().remove()
    ax3.get_legend().remove()
    if exp_type == 'bio':
        ax1.set_ylim(2, 8.5) # for BIO
        ax2.set_ylim(2, 8.5) # for BIO
        ax3.set_ylim(2, 8.5) # for BIO
    else:
        ax1.set_ylim(34, 55) # for SIM
        ax2.set_ylim(34, 55) # for SIM
        ax3.set_ylim(34, 55) # for SIM
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
    # g = sns.relplot(data=indata, x="dataset", y="f_score", hue="tool", kind="line", #dashes = dashes, style="type",
    #     col="type", col_order=["SNV", "INDEL"]) # hue="datastructure", style="datastructure",  col_wrap=3, )
    # g.set_axis_labels("Dataset", "F-score")
    # g.set(ylim=(0, 100))
    # plt.savefig(os.path.join(outfolder, "fscore_plot.pdf"))

    plt.close()



def plot_memory_usage(input_csv, outfolder, palette, tools):
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

def plot_runtime(input_csv, outfolder, palette, tools):
    sns.set(rc={'figure.figsize':(12,4)})
    matplotlib.rcParams.update({'font.size': 14})
    sns.set(font_scale=1.4)
    sns.set_style("whitegrid")
    # pyplot.figure(figsize=(4,16))
    indata = pd.read_csv(input_csv)
    # f, axes = plt.subplots(1, 1)
    ax1 = sns.lineplot(data=indata, x="dataset", y="time", hue="tool", hue_order = tools, palette=palette, linewidth = 2.0)
    ax1.set_ylabel("Time (s)")
    ax1.set_yscale('log')
    ax1.set_yticks([i for i in range(1000,9999,1000)] + [i for i in range(10000,100001,10000)]) #, ylim=(0, 5200))
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder, "runtime.pdf"))
    plt.clf()

def main(args):
    sns.set_style("whitegrid")
    palette = {
    'minimap2': 'tab:blue',
    'strobealign': 'tab:green',
    'bwa_mem': 'tab:orange',
    'accelalign': 'tab:red',
    'bowtie2' : 'tab:purple',
    'urmap' : 'tab:grey',
    'snap' : 'pink',    
    'bwa_mem2' : 'black',
    }
    tools =["minimap2", "bwa_mem", 'accelalign', "bowtie2", "snap", "bwa_mem2", "strobealign"] # "urmap", remove from experiemtns because we get an error in multithreading mode - therefore not fair against urmap to compare runtime with only one core 

    # plot_sv_calling_results(args.sv_csv, args.outfolder, palette, tools, args.type)
    plot_runtime(args.runtime_mem_csv, args.outfolder, palette, tools)
    # plot_memory_usage(args.runtime_mem_csv, args.outfolder, palette, tools)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sv_csv', type=str, default= "", help='results file')
    parser.add_argument('--runtime_mem_csv', type=str, default= "", help='results file')
    parser.add_argument('--type', type=str, default= "sim", help='bio or sim')
    parser.add_argument('outfolder', type=str,  help='outfolder to plots.')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)