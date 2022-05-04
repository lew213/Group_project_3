## genome scan pipeline, written by Jeff DaCosta, streamlined by Christian Sailer
## September 2016, updated 13 February 2017

import os, sys, argparse
from DS_classes_1 import G1

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='Start of genome scan pipeline. This script combines the two allele counte tables of one contrast, '+
                                 'removes fixed non-variant sites and calculates Dxy, Fst, pi and allele frequency difference. Missing data is '+
                                 'allowed. Outliers are selected based on a percentile cutoff (-cut, ppm in the file name).')

parser.add_argument('-i', type=str, metavar='inputdir_path', required=True, help='REQUIRED: Full or relative path to input AC table directory')
parser.add_argument('-coh1', type=str, metavar='cohort_1', required=True, help='REQUIRED: Name of first cohort')
parser.add_argument('-coh2', type=str, metavar='cohort_2', required=True, help='REQUIRED: Name of second cohort')
parser.add_argument('-snps', type=int, metavar='snps_per_window', default='100', help='Number of SNPs per window [100]')
parser.add_argument('-o', type=str, metavar='outputdir_path', default='na', help='Optional output directory relative path, contrast directory will be a subdir of this one')
parser.add_argument('-cut', type=float, metavar='cutoff_ratio', default='0.5', help='Top percentile, defines outliers [0.5]')
parser.add_argument('-per', type=float, metavar='percent_missing_data', default='0.2', help='Ratio of missing data allowed, eg 0.2 allows 80percent missing data [0.20]')
parser.add_argument('-suf', type=str, metavar='suffix_species', default='', help='Suffix to append to the file name, 2-letter species abbreviation, eg. _Aa')
parser.add_argument('-ploidy', type=int, metavar='ploidy', default=2, help='Ploidy of population [2]')
args = parser.parse_args()

#cwd = os.getcwd()
if args.o == 'na':
    outputdir = str(args.coh1+args.coh2)
else:
    outputdir = str(args.o+args.coh1+args.coh2)
    if os.path.exists(args.o) == False:
        os.mkdir(args.o)

# check if output directory exists, create it if necessary
if os.path.exists(outputdir) == False:
    os.mkdir(outputdir)
    print('\n\tCreated directory '+outputdir)
contrast = args.coh1+args.coh2

###### EXECUTION ######
# before repeating filtering and SNP window calculation, test if the according output files have been generated and if yes, jump to according steps
if __name__ == '__main__':
    genome_scan = G1(args=args)
    if os.path.isfile(outputdir+'/'+contrast+'_WG_'+str(args.snps)+'SNPs_3metrics'+args.suf+'.txt') == False:
        if os.path.isfile(outputdir+'/'+contrast+'_scaf_1_AC'+args.suf+'.table') == False:
            genome_scan.step1to6()
            
        else:
            print('\n\t SNPs have already been filtered for '+str(args.coh1)+str(args.coh2)+':')
            print('-> jump to')
            genome_scan.step4to6()
     
    else:
        print('\n\tPopulation genetic summary statistics have already been calculated for '+str(args.coh1)+str(args.coh2)+':')
        print('-> jump to')
        genome_scan.step5to6()
        
    print('\n Succesfully found selective sweep candidates for '+args.coh1+'_'+args.coh2+' contrast and '+str(args.snps)+' SNP windows\n\n')
    print('\n\tDONE\n')


