import allel
from matplotlib_venn import venn2, venn2_circles
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches


# Load SNPs
marker_excl_snp = allel.read_vcf(snakemake.input['isec_output_dir'] + '/' + '0000.vcf')
ref_excl_snp = allel.read_vcf(snakemake.input['isec_output_dir'] + '/' + '0001.vcf')
joined_snp = allel.read_vcf(snakemake.input['isec_output_dir'] + '/' + '0002.vcf')


# Function for counting variants
cal = lambda f: len(allel.GenotypeArray(f['calldata/GT']).count_alleles())


# construct plot
p = venn2(subsets=(cal(marker_excl_snp), cal(ref_excl_snp), cal(joined_snp)),
          set_labels=(marker_excl_snp['samples'][0],ref_excl_snp['samples'][0]),
          set_colors=('r', 'w'))
p = venn2_circles(subsets=(cal(marker_excl_snp), cal(ref_excl_snp), cal(joined_snp)),
                  linewidth=2)
legend_patch = mpatches.Patch(color='red', label='Complete marker positions', alpha=0.5)
plt.legend(handles=[legend_patch], bbox_to_anchor=(1, 1),
           bbox_transform=plt.gcf().transFigure, framealpha=0.5)
plt.title('Crossing id: {}'.format(snakemake.params['crossing_id']))
plt.savefig(snakemake.output['parental_snps_fig'])

# TODO: corrected markers, SNP/Marker density plots
