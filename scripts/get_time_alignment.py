
import os,sys
import argparse

import re

def parse_gnu_time(stderr_file):
    lines = open(stderr_file, 'r').readlines()
    index_time_mm2 = False
    index_time_aa = False
    index_time_strobemap = False

    for l in lines:
        wct_match = re.search('[\d.]+ real', l) 
        mem_match = re.search('[\d]+  maximum resident set size', l) 

        # minimap2
        index_time_mm2_match = re.search('\[M::main::[\d.:]+\*[\d.:]+\] loaded/built the index for', l) 
        
        # strobemap

        index_time_strobemap_match = re.search('Total time indexing: [\d.:]+', l) 

        # accelalign
        index_time_accelalign_match = re.search('Setup reference in [\d.:]+ secs', l) 

        if wct_match:
            wallclocktime = float(wct_match.group().split()[0])
        if mem_match:
            mem_tmp = int(mem_match.group().split()[0])
            memory_gb = mem_tmp / 1000000.0 

        if index_time_mm2_match:
            prefix_cut = index_time_mm2_match.group().split('main::')[1]
            final_time = prefix_cut.split('*')[0]
            index_time_mm2 = float(final_time.strip())

        if index_time_strobemap_match:
            # print(index_time_strobemap_match)
            index_time_strobemap = float(index_time_strobemap_match.group().split(':')[1].strip())
        
        if index_time_accelalign_match:
            prefix_cut = index_time_accelalign_match.group().split('reference in ')[1]
            final_time = prefix_cut.split(' secs')[0]
            index_time_aa = float(final_time.strip())

    if index_time_mm2:
        tot_wallclock_secs = wallclocktime - index_time_mm2
    elif index_time_aa:
        tot_wallclock_secs = wallclocktime - index_time_aa
    elif index_time_strobemap:
        tot_wallclock_secs = wallclocktime - index_time_strobemap

    return round(tot_wallclock_secs,2), round(memory_gb,2)

    # vals = list(map(lambda x: float(x), wallclocktime.split(".") ))
    # if len(vals) == 3:
    #     h,m,s = vals
    #     tot_wallclock_secs = h*3600.0 + m*60.0 + s
    #     if index_time_mm2:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_mm2
    #     elif index_time_aa:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_aa
    #     elif index_time_strobemap:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_strobemap

    # elif len(vals) == 2:
    #     s, = vals
    #     tot_wallclock_secs = m*60.0 + s

    #     if index_time_mm2:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_mm2
    #     elif index_time_strobemap:
    #         tot_wallclock_secs = tot_wallclock_secs - index_time_strobemap

        # return vals, memory_gb

    # else:
    #     return "BUG","BUG"

def main(args):
    tot_wallclock_secs, memory_gb = parse_gnu_time(args.infile)


    print("{0},{1}".format(tot_wallclock_secs, memory_gb))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('infile', type=str, default=None, help='Path to file.')
    # parser.set_defaults(which='main')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)