# Group_project_3

GATK variants to table commands were used to create individual population tables as instructed on the GitHub https://github.com/SailerChristian/Divergence_Scan – furthermore the three population tables of both population groups were concatenated together to create two end files, one for each population (saline or non-saline)

./gatk VariantsToTable \
-V /Users/lewiswood/Desktop/project3datanew/non_saline/Bruc2_populations.vcf \
-F CHROM -F POS -F AC -F AN \
-O /Users/lewiswood/Desktop/project3datanew/non_saline/bruc2_raw.table

First script to run was the outlier script – DS1_outliers.py. This is part one of the divergence scan pipeline. The script calls functions written in the classes.py script. Therefore, its best to keep all the scripts together in the same directory. The popgen metrics are based on non-overlapping SNPs to avoid inflating the fst values. The outlier cut-off in the exam command was set to 0.5 and the tolerance of missing data was set to 20%. The script worked fully with little complaints. However, the script will call multiple tables if they are present in the  directory (even in a sub directory) therefore file management is highly important with this pipeline. 

python3 /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/DS1_outliers.py -i /Users/lewiswood/Desktop/project3datanew -coh1 non_salin_all_raw -coh2 salin_all_raw -o end -snps 100 -cut 0.5 -per 0.2 -suf _Aa -ploidy 2

The second scripts DS2 and GS2 – appear to fulfil the same purpose. However, staying on track with the pipeline (https://github.com/SailerChristian/Divergence_Scan ) DS2 was ran first. Although the script identifies candidate genes, it has issues presenting the results. This is an issue with pandas data frame merging command towards the very end of the script. because of the two data frames having columns located in differing places the tables cannot be merged without first performing extensive data management. Due to this being a probable issue arising from the data files itself. This issue is unsolvable without reading the method of how the author of the pipeline produced the .gff and .txt files. The work around for this is to print the worked results straight to the terminal from the script or to use pandas to read them to csv files so they can be analysed later. this was done through editing the script and where it reads its data frames to. 

python3 /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/DS2_genes.py -i /Users/lewiswood/Desktop/endnon_salin_all_rawsalin_all_raw -coh1 non_salin_all_raw -coh2 salin_all_raw -an /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/models_aed0.5.gff -gf /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/LyV2_TAIR10orth_des_20150927.txt -ovlp 0.00001 -suf _Aa


going forward the script GS2 was called. This script successfully writes candidate gene lists however it doesn’t write the gene functions. The issue with this is that while fulfilling the aims of the project, it doesn’t allow for a specific conclusion to be drawn about the two population groups. It is also worth noting that the script called is an edited version. The bedtools calling section of the script that searches for overlaps was changed to accommodate for differences in vocabulary in the github files and the given project files. Therefore, the word transcript and transcription were changed to ‘gene’.



python3 /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/G2x_genes.py -i /Users/lewiswood/Desktop/project3datanew/endnon_salin_all_rawsalin_all_raw -an /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/models_aed0.5.gff -gf /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/LyV2_TAIR10orth_des_20150927.txt -ovlp 0.00001 -suf _Aa worked

 Furthermore, bed tools can also be ran outside of the script to identify overlaps. 	

bedtools intersect -a /Users/lewiswood/Desktop/project3datanew/project3datanewnon_salin_all_rawsalin_all_raw/non_salin_all_rawsalin_all_raw_100SNPs_5000ppm_Fst_outliers_Aa.bed -b /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/models_aed0.5.gff -f 1.0000000000000001e-07 -wb | awk '{$1=$2=$3=""; print $4,$5,$6,$7,$8,$9,$10,$11,$12}' | grep gene

the final script from the pipeline that produced results was GS3__graphs_percentile.py. the script was successful in creating r files with the necessary information to create graphs however missing elements from the r scripts require it to be re-written. In addition to this the r scripts are written with differencing assigned variables from the python code. Therefore, the scripts need to be further edited, however little annotation of the r scripts on the GitHub pipeline page makes this challenging. An attempt was made to cannibalise the AFD calculation portions of the script however due to the amount of information needed to be read into computer memory this was not possible due to the majority of the ram and available disk space on the laptop being taken up by previously created files. 

python3 /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/G3_graphs_percentile.py -coh1 non_salin_all_raw -coh2 salin_all_raw -geneor /Users/lewiswood/Desktop/project3datanew/Divergence_Scan-master/Lyrata2genesOriented.out -snps 100 -cut 0.5 -win 50000 -i /non_salin_all_rawsalin_all_raw -ovlp 0.00001 -outlier ALL -suf _Aa -ploidy 2

finally, the last script D3_perSNP.py, was not ran successfully. While there is not an explanation of why the script couldn’t process the data (even after all returned errors were fixed) it is likely that the work arounds employed in the previous parts of the pipelines caused issues downstream. The second explanation is that the given files may be annotated differently to those given on the GitHub pipeline https://github.com/SailerChristian/Divergence_Scan . 

In summary while the G3 and D3 scripts failed to work as expected – candidate genes lists and outliers were still produced. By going into outlier files and sorting for fst and cleaving these regions from fasta sequences. It was then possible to use blast to identify functions. A second similar approach to obtain gene function was to extent the ovlp – in the second script DS2 to 95%. This reduced the candidate gene list significantly, using a reference file it was then possible to identify the position of the candidate genes and later, their function. 
