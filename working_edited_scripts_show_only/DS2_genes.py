## genome scan pipeline, written by Jeff DaCosta, streamlined by Christian Sailer
## September 2016, updated 30 November 2016
## updated 11 June 2018 by Christian Sailer

import os, sys, argparse, subprocess, statistics
from natsort import natsorted
import pandas as pd
import numpy as np


#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Part 2 of divergence scan pipeline. This script identifies the genes that lie '+
                                 'within the given interval (-win) using the annotation gff file. In the second step '+
                                 'it greps the gene functions from the ortholog gene function list.')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full or relative path to directory containing the contrast directory (as in step DS1)')
parser.add_argument('-coh1', type=str, metavar='cohort_1', required=True, help='REQUIRED: Name of first cohort')
parser.add_argument('-coh2', type=str, metavar='cohort_2', required=True, help='REQUIRED: Name of second cohort')
parser.add_argument('-an', type=str, metavar='gff_annotation_file', default='Data/annotations/LyV2.gff', help='Full path to annotation gff file [Data/annotations/LyV2.gff]')
parser.add_argument('-gf', type=str, metavar='gene_function_file', default='Data/annotations/AL_AT_GI_OG_Araport11.txt', help='Full path to lyrata - thaliana orthologs file [Data/annotations/AL_AT_GI_OG_Araport11.txt]')
parser.add_argument('-ovlp', type=float, metavar='percent_overlap', default='0.0001', help='Precentage of base pairs that overlap the search pattern as percentage [0.0001]')
parser.add_argument('-suf', type=str, metavar='suffix_species', default='', help='Suffix to append to the file name, we recommend a 2-letter species abbreviation, eg. _Aa')
args = parser.parse_args()

print('\nSearching directories...')

contrast = args.coh1+args.coh2
outputdir = str(args.i+contrast+'/genes/')
print('Place output into '+str(outputdir))
if os.path.exists(outputdir) == False:
    os.mkdir(outputdir)

###### STEP 1 ######
# obtain outlier gene annotation list using bedtools
# search for overlaps of most diverged windos with gene annotation
print('\nStep 1: Obtain outlier gene annotation list, '+str(int(args.ovlp))+' percent overlap\n')
a = []
#files = []
for dirName, subdirList, fileList in os.walk(args.i+contrast+'/'):
    for file in fileList:
        if file.endswith('_outliers'+args.suf+'.bed') == True:
            a.append(file)
a_sorted = natsorted(a)

print('\tFound '+str(len(a))+' bedfiles to select gene annotation:')
count = 0
for file in a_sorted:
    print('\t Processing '+str(file))
    basename = file.replace('_outliers'+args.suf+'.bed', '')
    gcmd = open(outputdir+'bedtools_gff.unix', 'w')
    gcmd.write('bedtools intersect -a '+args.i+contrast+'/'+str(file)+' -b '+str(args.an)+' -f '+str(args.ovlp/100)+' -wb | ')
    gcmd.write("""awk '{$1=$2=$3=""; print $4,$5,$6,$7,$8,$9,$10,$11,$12}' """)
    gcmd.write('| grep gene | sort -u | ')
    gcmd.write("""tr ' ' '\t' """)
    gcmd.write('> '+outputdir+basename+'_'+str(int(args.ovlp))+'ol'+args.suf+'.gff')
    gcmd.close()

    #run in unix
    cmd = (open(outputdir+'bedtools_gff.unix', 'r'))
    p = subprocess.Popen(cmd, shell=True)
    sts = os.waitpid(p.pid, 0)[1]

    count +=1

print('\n\tExtracted gene annotation lists from '+str(count)+' bedfiles')
#os.remove(outputdir+'bedtools_gff.unix')


###### STEP 2 ######
# obtain outlier gene list
# extract the AL gene ID from GFF file
print('\nStep 2: Obtain outlier gene list, '+str(int(args.ovlp))+' percent overlap\n')
gff = []
for dirName, subdirList, fileList in os.walk(outputdir):
    for file in fileList:
        if file.endswith(str(int(args.ovlp))+'ol'+args.suf+'.gff') == True:
            gff.append(file)
gff_sorted = natsorted(gff)

print('\tFound '+str(len(gff))+' gff candidate files to select genes:')
count = 0
for file in gff_sorted:
    print('\t Processing '+str(file))
    basename = file.replace('.gff', '')
    gcmd = open(outputdir+'sed.unix', 'w')
    gcmd.write("""sed -n -e 's/^.*;Name=//p' """)
    gcmd.write(outputdir+file+' | ')
    gcmd.write("""sed -n -e 's/;.*$//p' """)
    gcmd.write('| sort -u > '+outputdir+basename+'.txt')
    gcmd.close()

    #run in unix
    cmd = (open(outputdir+'sed.unix', 'r'))
    p = subprocess.Popen(cmd, shell=True)
    sts = os.waitpid(p.pid, 0)[1]

    count +=1

print('\n\tExtracted gene lists from '+str(count)+' bedfiles')
#os.remove(outputdir+'sed.unix')


###### STEP 3 ######
# create interval list file in bed format for candidate genes
print('\nStep 3: Create interval bedfiles for candidate genes\n')
inlist = []
for dirName, subdirList, fileList in os.walk(outputdir):
    for file in fileList:
        if file.endswith(str(int(args.ovlp))+'ol'+args.suf+'.gff'):
            inlist.append(file)
inlist_sorted = natsorted(inlist)
print('\tFound '+str(len(inlist))+' gff candidate files to obtain intervals:')

# extract gene interval including 2kb upstream of gene
# the 2kb upstream are important for the follow up EAA
for file in inlist_sorted:
    print('\t Processing '+str(file))
    basename = file.replace(args.suf+'.gff','')
    with open(outputdir+file,'r') as infile:
        outfile = open(outputdir+basename+'_intervals'+args.suf+'.bed', 'w')
        for line in infile:
            data = line.split()
            datgen = data[8].split(sep=';')
            genename = datgen[0].split(sep='=')
            if data[6] == str('+'):
                bedstart = int(data[3])
                bedend = int(data[4])
            else:
                bedstart = int(data[3]) - 1
                bedend = int(data[4])
            outfile.write(data[0]+'\t'+str(bedstart)+'\t'+str(bedend)+'\t'+genename[1]+'\n')
    infile.close()
    outfile.close()

print('\n\tCreated '+str(count)+' candidate interval bedfiles')


###### Step 4 ######
# obtain gene onthologies from A. thaliana ortholog look up table
print('\nSTEP 4: Obtain orthologs and GF terms for gene lists\n')

a = []
for dirName, subdirList, fileList in os.walk(outputdir):
    for file in fileList:
        if file.endswith(str(int(args.ovlp))+'ol'+args.suf+'.txt') == True:
            a.append(file)
a_sorted = natsorted(a)

print('\tFound '+str(len(a))+' genelists to select genes:')

count = 0
for file in a_sorted:
    print('\t Processing '+str(file))
    with open(outputdir+file, 'r') as original:
        data = original.read()
    with open(outputdir+file, 'w') as modified:
        modified.write('AL_ID\n'+data)
    basename = file.replace(args.suf+'.txt', '')
    query = pd.read_table(outputdir+file)
    test = pd.read_table(args.gf)
    # merge tables as left join (left=AL_ID gene list)
    inner = pd.merge(test, query, how='left')
    # add contrast column
    inner['contrast'] = contrast
    inner.to_csv(outputdir+basename+args.suf+'_GF.txt', sep='\t', index=False)
        # # flexible header, takes header from input file
        # theader = test.readline()
        # header = theader.replace('\n', '')
        # outfile.write(header+'\tcontrast\n') # add one extra column
    #     for tline in query:
    #         line = tline.replace('\n', '')
    #         test = open(args.gf, 'r')
    #         for testline in test:
    #             ptestline = testline.replace('\n', '')
    #             data = testline.split('\t')
    #             if line in data[0]:
    #                 outfile.write(ptestline+'\t'+contrast+'\n')
    # query.close()
    # outfile.close()
    count +=1

print('\nObtained Arabidopsis thaliana ortholog gene function for '+str(count)+' outlier gene lists\n')
print('\n\tDONE\n')
