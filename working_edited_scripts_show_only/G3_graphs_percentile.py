## genome scan pipeline, part 3, written by Christian Saile
## genome scan pipeline, part 3, written by Christian Sailer
## updated 5 January 2017

import os, sys, subprocess, statistics, argparse
from natsort import natsorted
import pandas as pd
import numpy as np

# create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(
    description='Part 3 of genome scan pipeline. This script runs through the candidate gene list ' +
                'and creates an overview plot (AFD, Dxy, Fst, DD, AFDwindow) for each candidate.')

parser.add_argument('-coh1', type=str, metavar='cohort_1', required=True, help='REQUIRED: Name of cohort 1')
parser.add_argument('-coh2', type=str, metavar='cohort_2', required=True, help='REQUIRED: Name of cohort 2')
parser.add_argument('-geneor', type=str, metavar='gene_orientation file', required=True,
                    help='REQUIRED: Relative or full path to genes orientation file [Lyrata2genesOriented.out]')
parser.add_argument('-snps', type=str, metavar='SNP_window_size', default='100',
                    help='Number of SNPs in SNP-based window [100]')
parser.add_argument('-cut', type=float, metavar='cutoff_ratio', default='0.5',
                    help='Top percentile, defines outliers [0.5]')
parser.add_argument('-win', type=int, metavar='window_size', default='50000',
                    help='Number of bps plotted around candidate gene [50000]')
parser.add_argument('-i', type=str, metavar='inputdir_path', default='na',
                    help='Optional input directory if script is run from outside the GS_1 results folder [na]')
parser.add_argument('-ovlp', type=float, metavar='percent_overlap', default='0.0001',
                    help='Precentage of base pairs that overlap the search pattern as percentage [0.0001]')
parser.add_argument('-outlier', type=str, metavar='outliers_to_plot', default='na',
                    help='Define for which outlier combinations to produce graphs, eg DxyDD [na]')
parser.add_argument('-suf', type=str, metavar='species_suffix', default='',
                    help='Species abbreviation as suffix, eg _Aa for Arabidopsis arenosa')
parser.add_argument('-ploidy', type=int, metavar='ploidy', default=2, help='Ploidy of population [2]')

args = parser.parse_args()

if args.i == 'na':
    infile = str(args.coh1 + args.coh2 + '/')
else:
   # infile = str(args.i + args.coh1 + args.coh2 +'/')
   infile = str("/Users/lewiswood/Desktop/non_salin_all_rawsalin_all_raw")
print('\nSearching ' + infile + '...')

cwd = str(os.getcwd() + '/')

###### STEP 1 #######
# specify file based on SNP window size and CI and append all files within the contrast folder with that pattern
# load gene orientation file
f = pd.read_table(args.geneor, header=None)
f.columns = ['locus','scaf','start','stop','orient']
print('\nLoaded gene orientation file')
f['geneor'] = np.where(f.orient == '+', 'last', 'first')
print(f.loc[0:10])
print(f[f.locus == 'AL1G44270'])
winwidth = int(args.win / 1000 * 2)

search = str('_' + args.snps + 'SNPs_')

inall = []
for dirName, subdirList, filelist in os.walk(infile):
    for file in filelist:
        if search in file:
            inall.append(file)
suball = []
for file in inall:
    if str(search + str(int(10000 * args.cut)) + 'ppm_') in file:
        suball.append(file)
suball.sort()

# append AFs.table
AF = []
for dirName, subdirList, filelist in os.walk(infile):
    for file in filelist:
        #if file.endswith('_AFs' + args.suf + '.table'):
        if file.endswith('.table'):
            AF.append(file)
AF.sort()
print('\nLoaded ' + str(len(AF)) + ' AF tables:')
for i in range(len(AF)):
    print('\t' + str(AF[i]))

# append specified allsites.csv file
allsites = []
for file in suball:
    if file.endswith('_allsites' + args.suf + '.csv'):
        allsites.append(file)
allsites.sort()
for i in range(len(allsites)):
    print('\nLoaded ' + str(allsites[i]))

# append specified descriptive stats
descr = []
for file in suball:
    if file.endswith('_stats' + args.suf + '.txt'):
        descr.append(file)
descr.sort()
for i in range(len(descr)):
    print('\nLoaded ' + str(descr[i]))

# append gene list
genelist = []
if args.outlier == 'na':
    for file in suball:
        if file.endswith(str(int(args.ovlp)) + 'ol' + args.suf + '.txt'):
            genelist.append(file)
else:
    for file in suball:
        if file.endswith(args.outlier + '_' + str(int(args.ovlp)) + 'ol' + args.suf + '.txt'):
            genelist.append(file)
genelist.sort()

print('\nLoaded ' + str(len(genelist)) + ' genelists:')
for i in range(len(genelist)):
    print('\t' + str(genelist[i]))

# make dir strings
#inputdir = str(infile + 'genes/')
inputdir = str('/Users/lewiswood/Desktop/non_salin_all_rawsalin_all_rawgenes/')
contrast = inputdir.replace('/genes/', '')

for file in genelist:
    tbasename = file.replace(contrast + '_', '')
    basename = tbasename.replace(args.suf + '.txt', '')
    if os.path.exists(inputdir + 'plots_' + basename) == False:
        os.mkdir(inputdir + 'plots_' + basename)
    outputdir = str(inputdir + 'plots_' + basename + '/')
    inputfile = open(inputdir + file , 'r')
    #with open('/Users/lewiswood/Desktop/non_salin_all_rawsalin_all_rawgenes/non_salin_all_rawsalin_all_raw_100SNPs_5000ppm_DD_0ol_Aa.txt', 'w') as f:
    #inputfile = open("/Users/lewiswood/Desktop/non_salin_all_rawsalin_all_rawgenes/non_salin_all_rawsalin_all_raw_100SNPs_5000ppm_DD_0ol_Aa.txt", 'r')
    #inputfile = open("/Users/lewiswood/Desktop/non_salin_all_rawsalin_all_rawgenes/non_salin_all_rawsalin_all_raw_100SNPs_5000ppm_DD_0ol_Aa.txt")
    for line in inputfile:
         count = 0
         ttline = line.replace('\n', '')
         tline = str('"' + ttline + '"')
         scaf = tline[3]
         for i in range(len(allsites)):
            print('\nCreating overview plot for candidate gene\n\n\t' + ttline + '\n\nfrom list ' + str(
                inputdir + file) + '\nin ' + outputdir)
            #rfile = open(outputdir + ttline + '_' + basename + args.suf + '.r', 'w')
            rfile = open(outputdir + ttline + '_' + basename + args.suf + '.r', 'w')
            rfile.write('# load libraries\n' +
                        'library(ggplot2)\n' +
                        'library(grid)\n' +
                        'library(grid)\n' +
                        'f <- read.table("' + args.geneor + '",header=F)\n' +
                        'locus=' + tline + '\n' +
                        'winsize = ' + str(args.win) + '\n' +
                        'graphwindow <- winsize/1000*2\n' +
                        'locstart <- subset(f[,3], as.vector(f[,1])==locus)\n' +
                        'locend <- subset(f[,4], as.vector(f[,1])==locus)\n' +
                        'middle=(locend+locstart)/2\n' +
                        'locstart <- subset(f[,3], as.vector(f[,1])==locus)\n' +
                        'locend <- subset(f[,4], as.vector(f[,1])==locus)\n' +
                        'middle=(locend+locstart)/2\n' +
                        'start=middle-winsize\n' +
                        'end=middle+winsize\n' +
                        'genes=subset(f, (f[,2]==' + scaf + ' & f[,3]>=start & f[,3]<=end) | (f[,2]==' + scaf + ' & f[,4]>=start & f[,4]<=end))\n' +
                        'genes$geneor <- ifelse(genes[,5]=="+","last","first")\n\n' +
                        '# the path will be read in by -i\n' +
                        'ts25 <- read.csv("' + cwd + infile + str(allsites[i]) + '", header=T, stringsAsFactor=F)\n' +
                        '# for the scaffolds, the beginning of the name has to be evaluated first, best by the python script\n' +
                        's25 <- subset(ts25[ts25$scaffold=="scaffold_' + scaf + '", ])\n' +
                        'scaf' + scaf + ' <- read.table("' + cwd + infile + args.coh1 + args.coh2 + '_scaf_' + scaf + '_AFs' + args.suf + '.table", header=T)\n' +
                        '# read in descriptive stats table\n' +
                        'descr <- read.table("' + cwd + infile + str(descr[i]) + '", header=T)\n' +
                        '# select subsets around locus\n' +
                        'sr25 <- s25[s25$start >= start & s25$end <=end,]\n' +
                        'scaf' + scaf + '$rawAFD_2 <- scaf' + scaf + '$AC1/scaf' + scaf + '$AN1 - scaf' + scaf + '$AC2/scaf' + scaf + '$AN2\n' +
                        'scaf' + scaf + '_sub <- scaf' + scaf + '[scaf' + scaf + '$POS >= start & scaf' + scaf + '$POS <=end,]\n\n' +
                        '# switch to R only here\n' +
                        'layout <- theme_bw(base_size=10, base_family="Helvetica") +\n' +
                        '   theme(axis.title.x = element_blank(),\n' +
                        '   axis.text.x = element_blank(),\n' +
                        '   axis.ticks.x = element_blank(),\n' +
                        '   panel.grid = element_blank())\n\n' +
                        'layouttix <- theme_bw(base_size=10, base_family="Helvetica") +\n' +
                        '   theme(axis.title.x = element_blank(),\n' +
                        '   panel.grid = element_blank())\n\n' +
                        'pafds <- ggplot(aes(POS, rawAFD_2), data=scaf' + scaf + '_sub) + ylab("AFDsnp") +\n' +
                        '   ggtitle(locus) +\n' +
                        '   scale_y_continuous(limits=c(-1,1.1)) +\n' +
                        '   geom_point(alpha=0.35, size=0.5) +\n' +
                        'layouttix\n' +
                        '   for(j in 1:nrow(genes)){\n' +
                        '   pafds <- pafds + geom_segment(x=genes[j,3], xend=genes[j,4], y=1.1, yend=1.1, colour="grey30",\n' +
                        '   arrow=arrow(length=unit(0.02,"npc"), ends=genes[j,"geneor"]))\n' +
                        '}\n' +
                        'p.AFDs <- pafds  + geom_segment(x=genes[which(genes[,1]==locus), 3],\n' +
                        '   xend=genes[which(genes[,1]==locus), 4], y=1.1, yend=1.1, color="red", arrow=arrow(length=unit(0.03, "npc"),\n' +
                        '   ends=genes[which(genes[,1]==locus),"geneor"]))\n\n' +
                        'pdxy <- qplot(start, Dxy, data=sr25, geom="line", ylim=c(0,1)) +\n' +
                        '   geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4], ymin=-Inf, ymax=Inf), fill="grey90") +\n' +
                        '   geom_line(color=I("green")) +\n' +
                        '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff"), "Dxy"], color="green", linetype=2) +\n' +
                        '   layout\n')
            if args.ploidy == 2:
                rfile.write('pfst <- qplot(start, Fst, data=sr25, geom="line", ylim=c(0, 1)) +\n' +
                            '   geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4],\n' +
                            '   ymin=-Inf, ymax=Inf), fill="grey90") +\n' +
                            '   geom_line(color=I("blue")) +\n' +
                            '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff"), "Fst"], color="blue", linetype=2) +\n' +
                            '   layout\n')
                rfile.write(
                    'p.Fst <- pfst  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n' +
                    '   y=1, yend=1, color="red", arrow=arrow(length=unit(0.03, "npc"), ends=genes[which(genes[,1]==locus),"geneor"]))\n\n' +
                    'p.Dxy <- pdxy  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n' +
                    '   y=1, yend=1, color="red", arrow=arrow(length=unit(0.03, "npc"), ends=genes[which(genes[,1]==locus), "geneor"]))\n\n' +
                    'pdd <- qplot(start, DD, data=sr25, geom="line", ylim=c(-1, 1 )) +\n' +
                    '   geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4],\n' +
                    '   ymin=-Inf, ymax=Inf), fill="grey90") +\n' +
                    '   geom_line(color=I("goldenrod")) +\n' +
                    '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff"), "DD"], color="goldenrod", linetype=2) +\n' +
                    '   layout\n' +
                    'p.DD <- pdd  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n' +
                    '   y=25, yend=25, color="red", arrow=arrow(length=unit(0.03, "npc"), \n' +
                    '   ends=genes[which(genes[,1]==locus),"geneor"]))\n\n' +
                    'pafdw <- ggplot(aes(start, ' + args.coh1 + '_' + args.coh2 + '), data=sr25) + ylab("AFDwin") +\n' +
                    '   scale_y_continuous(limits=c(-1,1.1)) +\n' +
                    '   geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4],\n' +
                    '   ymin=-Inf, ymax=Inf), fill="grey90") +\n' +
                    '   geom_point(color="red", size=0.5) +\n' +
                    '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff"), "AFD_1_2"], color="red", linetype=2) +\n' +
                    '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff2"), "AFD_1_2"], color="red", linetype=2) +\n' +
                    '   layout\n' +
                    'p.AFDw <- pafdw  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n' +
                    '   y=1.1, yend=1.1, color="red", arrow=arrow(length=unit(0.03, "npc"),\n' +
                    '   ends=genes[which(genes[,1]==locus),"geneor"]))\n\n' +
                    'pdf("' + outputdir + ttline + '_' + basename + '_' + str(
                        winwidth) + 'kb' + args.suf + '.pdf", width=4, height=8)\n' +
                    'grid.newpage()\n' +
                    'pushViewport(viewport(layout = grid.layout(6,1)))\n' +
                    '# define number of columns and rows for the graph #\n' +
                    'vplayout <- function(x,y)\n' +
                    '   viewport(layout.pos.row=x,layout.pos.col=y)\n' +
                    '# specifiy which plot goes where #\n' +
                    'print(p.AFDs, vp=vplayout(1:2,1))\n' +
                    'print(p.Dxy, vp=vplayout(3,1))\n' +
                    'print(p.Fst, vp=vplayout(4,1))\n' +
                    'print(p.DD, vp=vplayout(5,1))\n' +
                    'print(p.AFDw, vp=vplayout(6,1))\n' +
                    'dev.off()')

            else:
                rfile.write('pfst <- qplot(start, Fst_4C, data=sr25, geom="line", ylim=c(0, 1)) +\n' +
                            '   geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4],\n' +
                            '   ymin=-Inf, ymax=Inf), fill="grey90") +\n' +
                            '   geom_line(color=I("blue")) +\n' +
                            '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff"), "Fst_4C"], color="blue", linetype=2) +\n' +
                            '   layout\n')
                rfile.write(
                    'p.Fst <- pfst  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n' +
                    '   y=1, yend=1, color="red", arrow=arrow(length=unit(0.03, "npc"), ends=genes[which(genes[,1]==locus),"geneor"]))\n\n' +
                    'p.Dxy <- pdxy  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n' +
                    '   y=1, yend=1, color="red", arrow=arrow(length=unit(0.03, "npc"), ends=genes[which(genes[,1]==locus), "geneor"]))\n\n' +
                    'pdd <- qplot(start, DD, data=sr25, geom="line", ylim=c(-1, 1 )) +\n' +
                    '   geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4],\n' +
                    '   ymin=-Inf, ymax=Inf), fill="grey90") +\n' +
                    '   geom_line(color=I("goldenrod")) +\n' +
                    '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff"), "DD"], color="goldenrod", linetype=2) +\n' +
                    '   layout\n' +
                    'p.DD <- pdd  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n' +
                    '   y=25, yend=25, color="red", arrow=arrow(length=unit(0.03, "npc"), \n' +
                    '   ends=genes[which(genes[,1]==locus),"geneor"]))\n\n' +
                    'pafdw <- ggplot(aes(start, ' + args.coh1 + '_' + args.coh2 + '), data=sr25) + ylab("AFDwin") +\n' +
                    '   scale_y_continuous(limits=c(-1,1.1)) +\n' +
                    '   geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4],\n' +
                    '   ymin=-Inf, ymax=Inf), fill="grey90") +\n' +
                    '   geom_point(color="red", size=0.5) +\n' +
                    '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff"), "AFD_1_2"], color="red", linetype=2) +\n' +
                    '   geom_hline(yintercept=descr[which(descr$stat=="Cutoff2"), "AFD_1_2"], color="red", linetype=2) +\n' +
                    '   layout\n' +
                    'p.AFDw <- pafdw  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n' +
                    '   y=1.1, yend=1.1, color="red", arrow=arrow(length=unit(0.03, "npc"),\n' +
                    '   ends=genes[which(genes[,1]==locus),"geneor"]))\n\n' +
                    'pdf("' + outputdir + ttline + '_' + basename + '_' + str(
                        winwidth) + 'kb' + args.suf + '.pdf", width=4, height=8)\n' +
                    'grid.newpage()\n' +
                    'pushViewport(viewport(layout = grid.layout(6,1)))\n' +
                    '# define number of columns and rows for the graph #\n' +
                    'vplayout <- function(x,y)\n' +
                    '   viewport(layout.pos.row=x,layout.pos.col=y)\n' +
                    '# specifiy which plot goes where #\n' +
                    'print(p.AFDs, vp=vplayout(1:2,1))\n' +
                    'print(p.Dxy, vp=vplayout(3,1))\n' +
                    'print(p.Fst, vp=vplayout(4,1))\n' +
                    'print(p.DD, vp=vplayout(5,1))\n' +
                    'print(p.AFDw, vp=vplayout(6,1))\n' +
                    'dev.off()')

            rfile.close()

            cmd = ('Rscript ' + outputdir + ttline + '_' + basename + args.suf + '.r')
            p = subprocess.Popen(cmd, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

            count += 1
    inputfile.close()
    #print('\n\tCreated ' + str(count) + ' graphs to annotate for ' + str(file))
    print('\n\tCreated ' + str(round) + ' graphs to annotate for ' + str(file))