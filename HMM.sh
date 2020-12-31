#! /usr/bin

python all_chr_p_0.01.py
awk -F'\t' '{if($3=="gene") print $1,$4,$5,$7}' Sc.gff > Sc_gene.gff
python Distance.py