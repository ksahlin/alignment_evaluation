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



def plot_runtime(input_csv, outfolder):
    sns.set(rc={'figure.figsize':(9,6)})
    matplotlib.rcParams.update({'font.size': 20})
    sns.set(font_scale=1.6)
    # tool,dataset,read_length,time,memory
    sns.set_style("whitegrid")

    indata = pd.read_csv(input_csv)

    # One liner to create a stacked bar chart.
    ax = sns.histplot(indata, x='dataset', hue='stage', weights='time',
                 multiple='stack', palette='tab10', shrink=0.8)
    ax.set_ylabel('Time')
    ax.set_xlabel('Dataset')
    # Fix the legend so it's not on top of the bars.
    legend = ax.get_legend()
    legend.set_bbox_to_anchor((1, 1))
    degrees = 70
    plt.xticks(rotation=degrees)
    plt.tight_layout()
    plt.rc('axes', labelsize=30)    # fontsize of the x and y labels
    plt.savefig(os.path.join(outfolder, "time_stages_plot.pdf"))
    plt.close()


def main(args):
    sns.set_style("whitegrid")
    plot_runtime(args.time_csv, args.outfolder)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('time_csv', type=str, help='results file')
    parser.add_argument('outfolder', type=str,  help='outfolder to plots.')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)