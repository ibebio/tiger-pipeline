#!/usr/bin/env python3
# coding: utf-8
import argparse

SD_FACTOR1 = 2.5
SD_FACTOR2 = 4

def main(args):
    import os
    import numpy as np
    import allel
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    def read_vcfs(vcf_list, fields):
        callset = {f: np.array([],dtype='float32') for f in fields }
        for vcf in vcf_list:
            current_callset = allel.read_vcf(vcf, fields=fields)
            for f in fields:
                callset[f] = np.append(callset[f], current_callset[f])
        return callset
        
    

    # Helper functions
    def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
        return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)
    def gauss(x,mu,sigma,A):
        return A*np.exp(-(x-mu)**2/2/sigma**2)


    def plot_dist_and_cutoff(field, cutoff, plot_prefix,title):
        out = plt.hist(callset['variants/{}'.format(field)], bins=100, label=field)
        plt.axvline(cutoff, color='black', linestyle=':', label='Cut-off')
        plt.title(title)
        plt.savefig('{}.{}.png'.format(plot_prefix,field))
        plt.close()

    def get_cutoff_stddev(field, sd_factor ,minus=False):

        if minus:
            cutoff = np.round(np.nanmean(callset['variants/{}'.format(field)]) - \
                          np.sqrt(np.nanmean( (np.nanmean(callset['variants/{}'.format(field)]) - \
                                               callset['variants/{}'.format(field)])**2)) * sd_factor)
        else:
            cutoff = np.round(np.nanmean(callset['variants/{}'.format(field)]) + \
                          np.sqrt(np.nanmean( (np.nanmean(callset['variants/{}'.format(field)]) - \
                                               callset['variants/{}'.format(field)])**2)) * sd_factor)
        print(field + ' ' + str(cutoff))
        return cutoff
                


    
    callset = read_vcfs(args.vcf, fields=['variants/QD', 'variants/FS', \
                                          'variants/MQ', 'variants/MQRankSum', \
                                          'variants/ReadPosRankSum'])

    # QD
    out = plt.hist(callset['variants/QD'], bins=100, label='QD')
    x = out[1]
    y = out[0]
    x= [ x[i] + (x[i + 1] - x[i])/2 for i in range(0,len(x) - 1)]
    params,cov = curve_fit(bimodal,x,y)
    plt.plot(x,bimodal(x,*params),color='red',lw=3,label='Fit')
    
    mean_1 = params[0]
    sd_1 = abs(params[1])
    mean_2 = params[3]
    sd_2 = abs(params[4])

    print("mean_1: {}, sd_1: {}, mean_2: {}, sd_2: {}".format(mean_1, sd_1, mean_2, sd_2))
                          
    if mean_1 >= mean_2:
        cutoff = mean_1 - SD_FACTOR1 * sd_1
    else:
        cutoff = mean_2 - SD_FACTOR1 * sd_2

    cutoff = np.floor(cutoff)
    plt.axvline(cutoff, color='black', linestyle=':', label='Cut-off')
    plt.title('QD < {QD:.1f}'.format(QD=cutoff))
    plt.savefig('{}.QD.png'.format(args.plot))
    plt.close()
    QD = cutoff

                
    # FS
    FS = get_cutoff_stddev('FS', SD_FACTOR1)
    plot_dist_and_cutoff('FS', FS, args.plot, title='FS > {FS:.1f}'.format(FS=FS))
              
    # MQ
    MQ=50.0

    # MQRankSum
    MQRankSum = np.nanquantile(callset['variants/MQRankSum'], 0.1)
    plot_dist_and_cutoff('MQRankSum', MQRankSum, args.plot, title='MQRankSum < {MQRankSum:.1f}'.format(MQRankSum=MQRankSum))
    
    # ReadPosRankSum
    ReadPosRankSum = get_cutoff_stddev('ReadPosRankSum', SD_FACTOR2, minus=True)
    plot_dist_and_cutoff('ReadPosRankSum', ReadPosRankSum, args.plot, title='ReadPosRankSum < {ReadPosRankSum:.1f}'.format(ReadPosRankSum=ReadPosRankSum))

    # Create filter string
    filter = "QD < {QD:.1f} || FS > {FS:.1f} || MQ < {MQ:.1f} || MQRankSum < {MQRankSum:.1f} || ReadPosRankSum < {ReadPosRankSum:.1f}".format(QD=QD, FS=FS, MQ=MQ, MQRankSum=MQRankSum, ReadPosRankSum=ReadPosRankSum)

    out_file = open(args.filter_file, 'w')
    out_file.write(filter)
    out_file.close()
    print(filter)


    
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Estimate filtering parameters from VCF file")

    parser.add_argument(
        "-vcf",
        "--vcf",
        help="VCF file. All variants (PASS or not PASS), are used! Multiple --vcf arguments are possible",
        required=True, action="append")
    parser.add_argument(
        "-f",
        "--filter_file",
        help="Output file with the filter line",
        required=True)
    parser.add_argument(
        "-p",
        "--plot",
        help="Plot file prefix",
        default=None,
        required=True
    )

    args = parser.parse_args()
    main(args)
