#!/bin/bash
python charles_filter_n.py
bedtools intersect -wa -header -sorted -f 0.5 -r -a charles_pass3.vcf -b output_sorted.vcf -g ../../merging_new_callings/human.hg19.genome > comparison_charles.vcf
uniq comparison_charles.vcf > comparison_charles_uniq.vcf
grep -vc '^#' comparison_charles_uniq.vcf
grep -vc '^#' output_sorted.vcf
grep -vc '^#' charles_pass3.vcf