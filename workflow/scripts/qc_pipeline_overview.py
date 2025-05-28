#!/usr/bin/env python3
# coding: utf-8
import jinja2
import subprocess
import cairosvg
# Input/output files
src_excl_snp = snakemake.input['isec_output_dir'] + '/' + '0000.vcf'
ref_excl_snp = snakemake.input['isec_output_dir'] + '/' + '0001.vcf'
joined_snp = snakemake.input['isec_output_dir'] + '/' + '0002.vcf'
noTEnoTLR_snp = snakemake.input['no_te_no_tlr_vcf']
overview_image = snakemake.output['overview_image']
# TODO: corrected markers, SNP/Marker density plots

# Collect names
crossing_id = snakemake.params['crossing_id']
with open(snakemake.input['src_parent']) as f:
    src_name = f.readline().strip()
with open(snakemake.input['ref_parent']) as f:
    ref_name = f.readline().strip()

# Collect statistics
# Get count of marker (src) exclusive SNPs
# We have biallelic SNPs only, thus, counting the line numbers is enough,
# and simpler than loading and requiring the clunky scikit allel library
src_ex_count = int(subprocess.run(
    'cat {} |grep -v ^# | wc -l'.format(src_excl_snp),
    shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))
ref_ex_count = int(subprocess.run(
    'cat {} |grep -v ^# | wc -l'.format(ref_excl_snp),
    shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))
joined_count = int(subprocess.run(
    'cat {} |grep -v ^# | wc -l'.format(joined_snp),
    shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))
noTEnoTLR_count = int(subprocess.run(
    'cat {} |grep -v ^# | wc -l'.format(noTEnoTLR_snp),
    shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))
f2_num = len([c["f2_samples"] for c in snakemake.config["crossings"] if c["id"] == crossing_id][0])



f2_filt = int(subprocess.run(
    "cat results/tiger_analysis/F2.{crossing_id}/beta_mixture_models/*.bmm.intersections.txt.done |"
    "grep '^0$' | wc -l".format(crossing_id=snakemake.params['crossing_id']),
    shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8'))

# Render svg output
env = jinja2.Environment(loader=jinja2.FileSystemLoader('resources/'))
t = env.get_template('tiger_template.svg')
svg_out = t.render(crossing_id=crossing_id,
         src_name=src_name,
         ref_name=ref_name,
         src_count=src_ex_count + joined_count,
         src_ex_count=src_ex_count,
         ref_count=ref_ex_count + joined_count,
         noTEnoTLR_count=noTEnoTLR_count,
         f2_num=f2_num,
         f2_filt=f2_filt)
cairosvg.svg2png(bytestring=svg_out, write_to=overview_image, scale=1.5)
