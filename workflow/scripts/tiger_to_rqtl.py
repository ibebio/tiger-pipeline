#!/usr/bin/env python
# coding: utf-8



# Manish Goel, MPIPZ
# fixed by Ilja Bezrukov and Uli Lutz, MPIDB, June 2020


import os
import numpy as np
import pandas as pd
import argparse
from collections import deque

def read_chromosome_size_file(chromosome_size_file):
    # Read the chromosome sizes file and return a dictionary
    # Assume they are integers!
    chromosome_sizes = {}
    with open(chromosome_size_file) as f:
        for line in f:
            (chr, size) = line.split()
            chromosome_sizes[int(chr)] = size
    return chromosome_sizes
    


# chr_size = {1:30427671,
#             2:19698289,
#             3:23459830,
#             4:18585056,
#             5:26975502
#             }

# total size = 119146348, 119.146348 Mb


# In[141]:

# Gather our code in a main() function
def main(args):
    ## CHANGE THIS PATH
    # samples_file='/Users/ulutz/owncloud/MPI/CRISPR/FLC/Experiments_final_set_63accessions/F2_mapping_AmelieSpazierer_2020/TIGER/path_to_TIGER_F2-41-1.txt'

    chr_size = read_chromosome_size_file(args.chrom_size_file)

    with open(args.samples_file) as f:
        samples = f.readlines()
    samples = [{'name': s[0], 'path': s[1]} for s in [s.split() for s in samples]]
    allcospos = {}
    co_breaks = {1: deque(),
                 2: deque(),
                 3: deque(),
                 4: deque(),
                 5: deque()
                 }

    # In[142]:


    for sample in samples:
        co_pos = pd.read_table(sample['path'], header=None, sep = "\s+|\t+|\s+\t+|\t+\s+", engine='python')
        co_pos['sample'] = sample['name']

        print(co_pos)
        corrected_co_pos = deque()
        # col 0 = chromosome
        for chr_co in co_pos.groupby(by=[0]):
            brks = [np.mean([chr_co[1].iat[i, 2], chr_co[1].iat[i+1, 1]], dtype=int) for i in range(chr_co[1].shape[0]-1)]
            co_breaks[chr_co[0]].extend(brks)
            chr_co[1].iloc[:-1, 2] = brks
            chr_co[1].iloc[1:, 1] = brks
            corrected_co_pos.append(chr_co[1])
        allcospos[sample['name']] = pd.concat(corrected_co_pos)


    # In[143]:


    co_df = pd.DataFrame([(k, v) for k, value in co_breaks.items() for v in np.unique(value)] + [(i,1) for i in range(1,6)] + [(i,chr_size[i]) for i in range(1,6)])
    co_df.sort_values([0, 1], inplace=True)
    cos_reg = deque()
    for chr_co in co_df.groupby([0]):
        cos_reg.append(pd.DataFrame({0:chr_co[0], 1:list(chr_co[1][1][:-1]), 2:list(chr_co[1][1][1:])}))
    cos_reg = pd.concat(cos_reg, ignore_index=True)
    cos_reg = cos_reg.astype('int64') # NEW

    # In[144]:


    genotypes = {}
    for sample in samples:
        gt = np.array(['-']*cos_reg.shape[0])

        for row in allcospos[sample['name']].itertuples(index=False):
            # import pdb ; pdb.set_trace()
            reg_loci = cos_reg.loc[(cos_reg[0]==row[0]) & (cos_reg[1]>=row[1]) & (cos_reg[2]<=row[2])].index.values
            if row[3]=='CC':
                gt[reg_loci] = 'A'
            elif row[3]=='CL':
                gt[reg_loci] = 'H'
            elif row[3]=='LL':
                gt[reg_loci] = 'B'
        genotypes[sample['name']] = gt


    # In[145]:


    genotypes = pd.DataFrame(genotypes)
    genotypes['reg']=cos_reg[0].astype(str)+':'+cos_reg[1].astype(str)+'_'+cos_reg[2].astype(str)
    genotypes['chr']=cos_reg[0]
    genotypes = genotypes[['reg', 'chr'] + list([s['name'] for s in samples])]
    genotypes = genotypes.transpose()


    # In[146]:


    genotypes.to_csv(args.output, header=False)


    # In[ ]:





    # In[ ]:

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert TIGER pipeline output to rQTL mapping input",
        epilog="As an alternative to the commandline, params can be placed in a file, one per line, and specified on the commandline like '%(prog)s @params.conf'.",
        fromfile_prefix_chars='@')


    parser.add_argument(
        "-s",
        "--samples_file",
        help="File name of the samples file. Format: SAMPLE_NAME\tTIGER_output_file",
        required=True)
    parser.add_argument(
        "-c",
        "--chrom_size_file",
        help="File with chromosome sizes. Format: Chromosome_number chromosome_size",
        required=True)
    parser.add_argument(
        "-o",
        "--output",
        help="Output file, in csv format.",
        default=None,
        required=True
    )

    args = parser.parse_args()
    main(args)
