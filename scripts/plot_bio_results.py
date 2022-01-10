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
    g = sns.relplot(data=indata, x="read_length", y="agreement", hue="tool", linewidth = linewidth, kind="line", hue_order = tools,  palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_axis_labels("Read length", "Agreement with BWA and Bowtie2")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(xticks=[150,250])
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

    g = sns.relplot( data=indata, x="read_length", y="percent_aligned", hue="tool", linewidth = linewidth, kind="line",  hue_order = tools, palette=palette)
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
        data=indata, x="read_length", y="time", hue="tool", linewidth = linewidth, kind="line", hue_order = tools,  palette=palette)
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


def plot_all(input_csv, outfolder, palette, tools, linewidth = 2.5):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.9)
    sns.set_style("whitegrid")
    indata = pd.read_csv(input_csv)
    g = sns.relplot(data=indata, x="read_length", y="data", hue="tool", linewidth = linewidth, kind="line", hue_order = tools,  palette=palette,
        col="type", col_wrap=3, col_order=["aligned", "time", "agreement"], facet_kws={'sharey': False, 'sharex': True} )
    # axes = g.axes
    # g.set_axis_labels("Read length", "Agreement with BWA and Bowtie2")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    # axes.set_xticks([18,24,30,36] )
    # ax.set_ylim((75, 100))
    g.set(xticks=[150,250])
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    plt.savefig(os.path.join(outfolder, "BIO_all_plot.eps"))
    plt.savefig(os.path.join(outfolder, "BIO_all_plot.pdf"))
    plt.close()


def plot_all2(input_csv, outfolder, palette, tools, linewidth = 2.5):
    matplotlib.rcParams.update({'font.size': 18})
    sns.set(font_scale=1.9)
    sns.set_style("whitegrid")
    indata = pd.read_csv(input_csv)

    fig, axs = plt.subplots(ncols=3, figsize=(20,8))
    xx = sns.lineplot(x='read_length', y='time', hue="tool", markers = True, palette = palette,  linewidth = linewidth, data=indata, ax=axs[0])
    yy = sns.lineplot(x='read_length', y='percent_aligned', markers = True,  palette = palette,  hue="tool", linewidth = linewidth, data=indata, ax=axs[1],legend=0)
    zz = sns.lineplot(x='read_length', y='agreement', hue="tool", palette = palette, hue_order = ["minimap2", 'accelalign',  "strobealign", "snap", "urmap"],
                         markers = True,  linewidth = linewidth, data=indata, ax=axs[2],legend=0)
    # # handles, labels = axs.get_legend_handles_labels()
    # handles, labels = [(a + b + c) for a, b, c in zip(axs[0].get_legend_handles_labels(), axs[1].get_legend_handles_labels(),axs[2].get_legend_handles_labels())]
    # fig.legend(handles, labels, loc='upper center')
    # g = sns.relplot(data=indata, x="read_length", y="data", hue="tool", linewidth = linewidth, kind="line", hue_order = tools,  palette=palette,
    #     col="type", col_wrap=3, col_order=["aligned", "time", "agreement"], facet_kws={'sharey': False, 'sharex': True} )
    # axes = g.axes
    # g.set_axis_labels("Read length", "Agreement with BWA and Bowtie2")
    # g.set_xticklabels([18,24,30,36])
    # ax.set_ylabel("% unique")
    # ax.set_xlabel("k")
    axs[0].set_xticks([150,250] )
    # print(dir(axs[0]))
    axs[0].set_yscale("log")
    axs[0].set_yticks([i for i in range(500,999,100)] + [i for i in range(1000,9999,1000)] + [i for i in range(10000,39999,10000)])

    axs[0].set_ylabel("Time (log scale)")

    axs[1].set_xticks([150,250] )
    axs[1].set_ylabel("Aligned (%)")
    xt = (0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0)
    axs[1].set_yticks(xt)
    axs[1].set_yticklabels(([str(int(l*100)) for l in xt]))

    axs[2].set_xticks([150,250] )
    axs[2].set_ylabel("BWA & bowtie2 agreement (Million reads)")
    axs[2].set_ylim((7100000, 7500000))
    axs[2].set_yticks((7100000, 7200000, 7300000, 7400000, 7500000))
    axs[2].set_yticklabels(("7.1", "7.2", "7,3", "7.4", "7.5"))


    # box1 = xx.get_position()
    # xx.set_position([box1.x0, box1.y0, box1.width * 0.85, box1.height]) # resize position
    # box2 = yy.get_position()
    # yy.set_position([box2.x0, box2.y0, box2.width * 0.85, box2.height]) # resize position
    # box3 = zz.get_position()
    # zz.set_position([box3.x0, box3.y0, box3.width * 0.85, box3.height]) # resize position
    # zz.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

    # g.set(xticks=[150,250])
    # g.set(ylim=(95, 100))
    # ax.set_xticks([18,24,30,36])
    # fig.legend(loc=7)
    fig.tight_layout()
    # fig.subplots_adjust(right=0.75)   
    plt.savefig(os.path.join(outfolder, "BIO_all_plot.eps"))
    plt.savefig(os.path.join(outfolder, "BIO_all_plot.pdf"))
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
    tools =["minimap2", "bwa_mem", 'accelalign', "bowtie2",  "strobealign", "snap", "urmap"]
    tools_agree =["minimap2", 'accelalign', "strobealign", "snap"] # "urmap",

    # accuracy_csv = add_column(args.accuracy_csv)
    # runtime_mem_csv = add_column(args.runtime_mem_csv)
    plot_all2(args.csv, args.outfolder, palette, tools, linewidth = 2.5)

    # plot_percentage_aligned(args.csv, args.outfolder, palette, tools, linewidth = 2.5)
    # plot_runtime(args.csv, args.outfolder, palette, tools, linewidth = 2.5)
    # plot_agreement(args.csv, args.outfolder, palette, tools_agree, linewidth = 2.5)


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




# tool,dataset,read_length,memory,data,type
# minimap2,BIO150,150,11.431484,0.98765,aligned
# minimap2,BIO150,150,11.431484,7388427,agreement
# minimap2,BIO150,150,11.431484,1591.1380000000001,time
# bwa_mem,BIO150,150,5.468372,0.99362,aligned
# bwa_mem,BIO150,150,5.468372,0,agreement
# bwa_mem,BIO150,150,5.468372,4676.0,time
# strobealign,BIO150,150,35.75272,0.9902,aligned
# strobealign,BIO150,150,35.75272,7409251,agreement
# strobealign,BIO150,150,35.75272,799.4620000000001,time
# accelalign,BIO150,150,19.305676,0.96928,aligned
# accelalign,BIO150,150,19.305676,7233470,agreement
# accelalign,BIO150,150,19.305676,738.45,time
# bowtie2,BIO150,150,3.402236,0.99374,aligned
# bowtie2,BIO150,150,3.402236,0,agreement
# bowtie2,BIO150,150,3.402236,7930.0,time
# snap,BIO150,150,33.587392,0.97539,aligned
# snap,BIO150,150,33.587392,7384415,agreement
# snap,BIO150,150,33.587392,1057.06,time
# urmap,BIO150,150,29.426892,0.9743,aligned
# urmap,BIO150,150,29.426892,7400486,agreement
# urmap,BIO150,150,29.426892,742.76,time
# minimap2,BIO250,250,11.372176,0.99268,aligned
# minimap2,BIO250,250,11.372176,7431920,agreement
# minimap2,BIO250,250,11.372176,2517.722,time
# bwa_mem,BIO250,250,5.415552,0.99713,aligned
# bwa_mem,BIO250,250,5.415552,0,agreement
# bwa_mem,BIO250,250,5.415552,7861.0,time
# strobealign,BIO250,250,35.917084,0.99258,aligned
# strobealign,BIO250,250,35.917084,7427035,agreement
# strobealign,BIO250,250,35.917084,970.3819999999998,time
# accelalign,BIO250,250,19.30904,0.98596,aligned
# accelalign,BIO250,250,19.30904,7303527,agreement
# accelalign,BIO250,250,19.30904,1077.33,time
# bowtie2,BIO250,250,3.419988,0.99464,aligned
# bowtie2,BIO250,250,3.419988,0,agreement
# bowtie2,BIO250,250,3.419988,19368.0,time
# snap,BIO250,250,33.586308,0.93357,aligned
# snap,BIO250,250,33.586308,7166512,agreement
# snap,BIO250,250,33.586308,2807.45,time
# urmap,BIO250,250,29.427496,0.95306,aligned
# urmap,BIO250,250,29.427496,7306319,agreement
# urmap,BIO250,250,29.427496,1422.84,time

