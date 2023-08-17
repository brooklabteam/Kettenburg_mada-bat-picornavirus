Building a maximum-likelihood Bayesian timetree using BEAST2
---














```
cd Desktop/developer/mada-bat-picornavirus/TreeTime/picornavirales

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_rooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_rooted_skyline_out --coalescent skyline

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_out

treetime --aln picornavirales_all_align.fasta --tree picornavirales_all_unrooted --dates picornavirales_all_treetime_date.csv --outdir picornavirales_all_unrooted_skyline_out --coalescent skyline

```
