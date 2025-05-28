#!/usr/bin/env python3
# coding: utf-8
import argparse

def main(args):
    import os
    import sys
    import numpy as np
    import pandas as pd
    import subprocess

    # Get header length
    header_length = int(subprocess.run(
        'cat {} |head -1000 | grep "^@" | wc -l'.format(args.input_tab),
        shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))
    # Read input table
    tab = pd.read_table(args.input_tab, header=header_length)

    # Remove Mt and Pt
    tab = tab.drop(tab[(tab.CONTIG == 'Pt') | (tab.CONTIG == 'Mt')].index)
    len_tab_orig = len(tab)
    # Only count non-zero occurences
    tab_nonzero = tab.drop(tab[(tab.ALT_COUNT + tab.REF_COUNT) == 0].index)
    # Compute mean and sd
    mean = np.nanmean(tab_nonzero.REF_COUNT + tab_nonzero.ALT_COUNT)
    sd = np.std(tab_nonzero.REF_COUNT + tab_nonzero.ALT_COUNT)
    # Apply filtering
    tab = tab.drop(tab[(tab.REF_COUNT + tab.ALT_COUNT) >= mean + args.sd_factor * sd].index)

    # Print out statistics
    len_tab = len(tab)
    len_tab_nonzero = len(tab_nonzero)
    header_string = "#SAMPLE\tCOUNT\tNON_ZERO_COUNT\tFILTERED_COUNT\tMEAN_NON_ZERO\tSD_NON_ZERO"
    stat_string = "{sample}\t{len_tab_orig:d}\t{len_tab_nonzero:d}\t{len_tab:d}\t{mean:.1f}\t{sd:.1f}".format(
        sample=args.name, len_tab_orig=len_tab_orig, len_tab_nonzero=len_tab_nonzero, len_tab=len_tab, mean=mean, sd=sd)
    for s in [sys.stdout, sys.stderr]:
        print(header_string, file=s)
        print(stat_string, file=s)

    # Write output
    tab.to_csv(args.output_tab,
               sep='\t',
               header=True,
               index=False,
    )

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter allelic counts to remove outliers")

    parser.add_argument(
        "-n",
        "--name",
        help="sample name",
        required=False, default="",
    )
    parser.add_argument(
        "-i",
        "--input_tab",
        help=".tab file produced by gatk CollectAllelicCounds",
        required=True)
    parser.add_argument(
        "-o",
        "--output_tab",
        help="Output tab file",
        required=True)
    parser.add_argument(
        "-s",
        "--sd_factor",
        help="SD factor for filtering. Counts which are above sd_factor * mean are removed",
        default=None,
        required=True,
        type=float
    )

    args = parser.parse_args()
    main(args)
