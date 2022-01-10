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


def plot_agreement(input_csv, outfolder, palette, tools, linewidth = 2.5):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.9)
    sns.set_style("whitegrid")
    indata = pd.read_csv(input_csv)
    g = sns.relplot(data=indata, x="read_length", y="accuracy", hue="tool", linewidth = linewidth, kind="line", #dashes = dashes,
        col="dataset", hue_order = tools,  palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_axis_labels("Read length", "Agreement with BWA and Bowtie2")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(ylim=(96, 99), xticks=[150,250])
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "BIO_accuracy_plot.eps"))
    plt.savefig(os.path.join(outfolder, "BIO_accuracy_plot.pdf"))
    plt.close()

def plot_percentage_aligned(input_csv, outfolder, palette, tools, linewidth = 2.5):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.9)
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)

    g = sns.relplot(
        data=indata, x="read_length", y="percent_aligned", hue="tool", linewidth = linewidth,
        col="dataset", kind="line",  hue_order = tools, palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_axis_labels("Read length", "Percentage aligned")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(ylim=(90, 100), xticks=[150,250])
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "BIO_percentage_aligned_plot.eps"))
    plt.savefig(os.path.join(outfolder, "BIO_percentage_aligned_plot.pdf"))
    plt.close()


def plot_runtime(input_csv, outfolder, palette, tools, linewidth = 2.5):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.9)
    # tool,dataset,read_length,time,memory
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)
    g = sns.relplot(
        data=indata, x="read_length", y="time", hue="tool", linewidth = linewidth,
        col="dataset", kind="line", hue_order = tools,  palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_axis_labels("Read length", "Time (sec)")
    g.set(yscale="log")
    g.set( yticks= [i for i in range(100,999,100)] + [i for i in range(1000,9999,1000)] + [i for i in range(10000,39999,10000)]) #, ylim=(0, 5200))
    # g.set_yticklabels( ["100"] + ["" for i in range(200,999,100)] + ["1000"] +  ["" for i in range(2000,9999,1000)] + ["10000"] +  [i for i in range(1000,2999,1000)]])

    g.set( xticks=[150,250]) #ylim=(40, 100),
    # g.set( yticks=[0,1000,2000,4000,6000,12000,18000,24000]) #ylim=(40, 100),

    # g.set(ylim=(95, 100))s
    # ax.set_xticks([18,24,30,36])
    # plt.savefig(os.path.join(outfolder, "BIOs_time_plot.eps"))
    plt.savefig(os.path.join(outfolder, "BIOtime_plot.pdf"))
    plt.close()

def add_column(infile):
    mod_outfile = open(infile + "_mod.csv", "w")
    for i, line in enumerate(open(infile,'r')):
        if i == 0:
            line = line.strip()+",type\n"
            mod_outfile.write(line)
            continue
        vals = line.strip().split(",")
        is_aln = True
        if "_map" in vals[0]:
            v_tmp = vals[0][:-4]
            vals[0] = v_tmp
            is_aln = False

        if is_aln:
            vals.append("align")
        else:
            vals.append("map")

        mod_line = ",".join([v for v in vals]) + "\n"
        mod_outfile.write(mod_line)
    mod_outfile.close()
    return mod_outfile.name



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
    tools =["minimap2", "bwa_mem", 'accelalign', "bowtie2",  "strobealign", "snap"] # "urmap",

    # accuracy_csv = add_column(args.accuracy_csv)
    # runtime_mem_csv = add_column(args.runtime_mem_csv)

    plot_percentage_aligned(args.csv, args.outfolder, palette, tools, linewidth = 2.5)
    plot_runtime(args.csv, args.outfolder, palette, tools, linewidth = 2.5)
    plot_agreement(args.csv, args.outfolder, palette, tools, linewidth = 2.5)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('csv', type=str, help='results file')
    parser.add_argument('outfolder', type=str,  help='outfolder to plots.')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)


# tool,dataset,read_length,time,memory,percent_aligned,agreement
# minimap2,BIO150,150,1591.1380000000001,11.431484,0.98765,7388427
# bwa_mem,BIO150,150,4676.0,5.468372,0.99362,0
# strobealign,BIO150,150,799.4620000000001,35.75272,0.9902,7409251
# accelalign,BIO150,150,738.45,19.305676,0.96928,7233470
# bowtie2,BIO150,150,7930.0,3.402236,0.99374,0
# snap,BIO150,150,1057.06,33.587392,0.97539,7384415
# urmap,BIO150,150,742.76,29.426892,0.9743,7400486
# minimap2,BIO250,250,2517.722,11.372176,0.99268,7431920
# bwa_mem,BIO250,250,7861.0,5.415552,0.99713,0
# strobealign,BIO250,250,970.3819999999998,35.917084,0.99258,7427035
# accelalign,BIO250,250,1077.33,19.30904,0.98596,7303527
# bowtie2,BIO250,250,19368.0,3.419988,0.99464,0
# snap,BIO250,250,2807.45,33.586308,0.93357,7166512
# urmap,BIO250,250,1422.84,29.427496,0.95306,7306319





# strobealign_map,BIO150,150,487.591,35.752552,
# minimap2_map,BIO150,150,1212.958,11.339444
# accelalign_map,BIO150,150,747.45,19.305656

# strobealign_map,BIO250,250,476.06199999999995,35.752712
# minimap2_map,BIO250,250,1636.795,11.296968
# accelalign_map,BIO250,250,1072.99,19.308