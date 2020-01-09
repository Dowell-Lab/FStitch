#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import numpy as np
import subprocess
import io
import argparse
import os
import sys
import datetime
from pandas import DataFrame
from io import StringIO
from argparse import RawTextHelpFormatter

def main():

    parser = argparse.ArgumentParser(description='FStitch Annotation Module\n\n===============================================================================================================\n\nExpands gene regions based on transcriptional annotations derived from the segment module.', epilog='@Dowell Lab, Margaret Gruca, margaret.gruca@colorado.edu\nFor questions and issues, see https://github.com/Dowell-Lab/FStitch', usage='%(prog)s --cram sample.cram --bed fstitch.bed --genes gene_ref.bed --output /my/out/dir', formatter_class=RawTextHelpFormatter)
    
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    required.add_argument('-c', '--cram', dest='cram', metavar='SAMPLE.CRAM/BAM', \
                       help='Mapped CRAM/BAM file. Both file types supported, although CRAM is recommended. Will need to export path to reference genome if using CRAM.', required=True)
    
    required.add_argument('-g', '--genes', dest='gene_ref', metavar='GENE_REF.BED', \
                       help='Gene reference file in BED6(+) format. Only the first 6 columns will be used which must be in standard chr, start, end, name, score, strand format.', required=True)

    required.add_argument('-b', '--bed', dest='fstitch_seg_file', metavar='FSTITCH.BED', \
                       help='FStitch segment output (BED file), concatenated for both postive and negative strands.', required=True)
    
    required.add_argument('-o', '--output', dest='output', metavar='/path/to/out/dir', \
                        help='Directory where tmp files and output BED file will be saved. Full path must be specified.', required=True)
    
    optional.add_argument('-s', '--save', dest='save_raw_annotations', action='store_true', \
                        help='Save list of active genes and original annotated regions (prior to merging/expanding using FStitch). Default = False', default=False, required=False)
    
    optional.add_argument('-r', '--radius', dest='radius', metavar='<RADIUS>', \
                       help='Radius around the transcription start site (TSS) for which counts will be generated to determine activity. Default = 1500 (Recommended).', default=1500, required=False, type=int)
    
#    optional.add_argument('-d', '--distance', dest='merge_distance', metavar='<DISTANCE>', \
#                       help='Distance (in bp) in which segment data will be merged (see BEDTools merge for additional details, same as -d argument). Default = 100 (Recommended).', default=100, required=False, type=int)    
    
    optional.add_argument('-min', '--mincount', dest='count', metavar='<COUNT>', \
                       help='Minimum number of read counts for site to be considered "active". Depth/complexity of data may determine an adjustment. Default = 10.', default=10, required=False, type=int)    
    
    args = parser.parse_args()
    
##### Set any global variables
    RADIUS = args.radius
    COUNT = args.count
    rootname_path = os.path.splitext(os.path.basename(args.cram))[0]
    rootname = rootname_path.split('.', 1)[0]
    cram = args.cram
    fstitch = args.fstitch_seg_file
    outdir = args.output    
    
    pd.options.mode.chained_assignment = None    
    
    try:
        if os.stat(args.fstitch_seg_file).st_size > 0:
            print("Beginning Processing....." + str(datetime.datetime.now()))  
        else:
            raise ValueError("FStitch BED file is empty. Check your FStitch training file for errors occuring in FStitch segment module.")
    except OSError:
        raise ValueError("FStitch BED file missing ... exiting")      
    
##### Generate a dataframe for the full gene reference file and make a new one which expands TSS's for counting
    
    genes = pd.read_csv((args.gene_ref), sep='\t', usecols = range(6), \
                         names = ['chr', 'start', 'end', 'gene', 'score', 'strand'])
    
    genes_tss = pd.read_csv((args.gene_ref), sep='\t', usecols = range(6), \
                         names = ['chr', 'start', 'end', 'gene', 'score', 'strand'])
    
    genes_tss.loc[genes_tss['strand'] == '-', 'start'] = genes_tss['end'] - RADIUS
    genes_tss.loc[genes_tss['strand'] == '-', 'end'] = genes_tss['end'] + RADIUS
    
    genes_tss.loc[genes_tss['strand'] == '+', 'end'] = genes_tss['start'] + RADIUS
    genes_tss.loc[genes_tss['strand'] == '+', 'start'] = genes_tss['start'] - RADIUS
    
    # Make sure there are no start values less than 0 and if so reset them to 0
    genes_tss['start'] = np.where(genes_tss['start'] < 0, 0, genes_tss['start'])
    genes_tss['end'] = np.where(genes_tss['end'] < 0, 0, genes_tss['end'])
    genes_tss.to_csv('%s/.%s_tss_file.bed' % (outdir, rootname), header=None, sep='\t', index=False)    ### Need to figure out a temp directory/file solution for these potentially
    
#########################################################################################################    
    
##### Generate counts on the opposite strand to determine TSS
    # Need to check and make sure there won't be false negatives due to unidiretional TSS -- yet to see this but who knows. Can always set a flag to count on both strands
    
    print("Generating counts for TSS.....")    
    
    all_tss = '%s/.%s_tss_file.bed' % (outdir, rootname) ### Again... need to figure out if there's a different temp file solution here
    
    count = subprocess.Popen(['multiBamCov', '-bams', cram, '-bed', all_tss, '-S'],
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)
    
    countout,counterr = count.communicate()
    
    # Convert count output data to a dataframe, clean for CRAM related warnings, and drop anything with counts less than 10 (this can be a flag, too)
    print("mutliBamCov stderr output:")
    print(counterr)
    countdata = StringIO(str(countout,'utf-8'))
    count_df = pd.read_csv(countdata, sep='\t', header=None, \
                      names=['chr', 'start', 'end', 'gene', 'score', 'strand', 'count'])
    # This next step will take care of CRAM warning messages if they come up and prevent an errors in the subsequent steps
    clean_count_df = count_df.dropna()
    
    print("Counting done for TSS. Parsing data to determine active genes.....")       
    
#########################################################################################################     

##### Determine which genes are actively transcribed
    # This can be adjusted using an argument
    transcribed_genes_tss = clean_count_df[clean_count_df['count'] >= COUNT]
    # This is getting turned into a float... let's fix that for bedtools moving forward
    transcribed_genes_tss['start'] = transcribed_genes_tss['start'].astype(int)
    transcribed_genes_tss['end'] = transcribed_genes_tss['end'].astype(int)
    transcribed_genes_tss.to_csv('%s/.%s_transcribed_genes_tss.bed' % (outdir, rootname), header=None, sep='\t', index=False)  ### Need to figure out a temp directory/file solution for these potentially
    
    # Generate a list of "active" genes
    transcribed_genes_list = transcribed_genes_tss['gene'].tolist()
    # Grab gene annotations from original refseq file that will be merged with FStitch regions
    transcribed_gene_annotations = genes[genes['gene'].isin(transcribed_genes_list)]
    if (args.save_raw_annotations):
        transcribed_gene_annotations.to_csv('%s/%s_transcribed_gene_annotations.bed' % (outdir, rootname), header=None, sep='\t', index=False)
        
    transcribed_gene_annotations['original_length'] = transcribed_gene_annotations['end'] - transcribed_gene_annotations['start'] # Get original length for a stats file output & filter
    
##########################################################################################################     

##### Subtract any active TSS from FStitch regions. This will ensure the gene annotation ends before the next active gene begins if they run into one another
    tss = '%s/.%s_transcribed_genes_tss.bed' % (outdir, rootname)

    subtract = subprocess.Popen(['subtractBed', '-a', fstitch, '-b', tss],
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)
    
    subout,suberr = subtract.communicate()
    print("subtractBed stderr output:") # This is annoying to separate these but kept running into a subprocess error if I put them together if the suberr was empty
    print(suberr)
    
    subdata = StringIO(str(subout,'utf-8'))
    sub_df = pd.read_csv(subdata, sep='\t', header=None, index_col=False, \
                      names=['chr', 'start', 'end', 'gene', 'score', 'strand', 'count', 'igv'])
    
    fstitch_seg_file = sub_df[~sub_df.gene.str.contains("OFF")]
    fstitch_seg_file = fstitch_seg_file.drop(columns=['igv', 'count'])
    
    # Concatenate FStitch segments with genes for merging
    concat = pd.concat([fstitch_seg_file, transcribed_gene_annotations], sort=False)
    sorted_concat = concat.sort_values(by=['chr', 'start'])
    sorted_concat.to_csv('%s/.%s_all_regions.bed' % (outdir, rootname), header=None, sep='\t', index=False) ### Yet another temp file...
    
#########################################################################################################         

##### Merge genes with FStitch regions
    regions = '%s/.%s_all_regions.bed' % (outdir, rootname)
    #merge_d = args.merge_distance
    
    merge = subprocess.Popen(['mergeBed', '-i', regions, '-s', '-d', '100' , '-c', '4,6', '-o', 'collapse'],
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT) 
    
    mergeout,mergeerr = merge.communicate()
    print("mergeBed stderr output:")
    print(mergeerr)
    
    mergedata = StringIO(str(mergeout,'utf-8'))
    merge_df = pd.read_csv(mergedata, sep='\t', header=None, \
                names = ['chr', 'start', 'end', 'name', 'strand'])

#########################################################################################################        

##### Data clean-up -- get correct BED formatting with single unqiue entry per gene
    split_df = DataFrame(merge_df.name.str.split(',').tolist(), index=[merge_df.chr, merge_df.start, merge_df.end, merge_df.strand]).stack() # Create one unique entry per gene/region
    split_df = split_df.reset_index()[[0, 'chr', 'start', 'end', 'strand']] # name variable is currently labeled 0
    split_df.columns = ['gene', 'chr', 'start', 'end', 'strand'] # renaming "name" to gene
    unique_genes = split_df[~split_df.gene.str.contains('=')].sort_values(by=['chr', 'start']) # Drop any residual FStitch regions to just get back a list of genes

    # This next part is to parse the annotations accordingly as to whether the user gave accession_name2 or just accession
    row = unique_genes.iloc[1, unique_genes.columns.get_loc('gene')]

    if (row.count('_') <= 1) :
        unique_genes['score'] = 0
        final_gene_annotations = unique_genes[['chr', 'start', 'end', 'gene', 'score', 'strand']]        
    elif (row.count('_') > 1) :
        gene_accession_split = unique_genes['gene'].str.replace(r'([^_]+_[^_]+)_', r'\1|').str.split('|', expand=True).rename(lambda x: f'col{x + 1}', axis=1) # Split gene/accession number       
        final_gene_annotations = unique_genes.join(gene_accession_split).drop(columns=['gene']).rename(columns={'col1': 'accession', 'col2': 'gene'}) # Add gene/accession number back in as separate columns, drop combined accesssion_gene column
        final_gene_annotations['strand'] = final_gene_annotations['strand'].str.split(',').str[0] # Remove list of strands from all merged regions and drop to single identifier
        final_gene_annotations = final_gene_annotations[['chr', 'start', 'end', 'gene', 'accession', 'strand']] # Reorder columns to make a pseudo BED6
        
    final_gene_annotations.to_csv('%s/%s_expanded_gene_annotations.bed' % (outdir, rootname), header=None, sep='\t', index=False)            
        
    gene_lengths_stats = transcribed_gene_annotations[transcribed_gene_annotations['gene'].isin(unique_genes['gene'].tolist())].sort_values(by=['chr', 'start'])
    gene_lengths_stats['end_length'] = (unique_genes['end'] - unique_genes['start']).tolist()
    gene_lengths_stats.to_csv('%s/%s_annotated_gene_length_stats.bed' % (outdir, rootname), sep='\t', index=False)     
    
    clean = subprocess.Popen(['rm', regions, all_tss, tss],
           stdout=subprocess.PIPE, 
           stderr=subprocess.STDOUT)
    
    print('Annotation module complete.\n' + str(datetime.datetime.now()))    
    
if __name__=='__main__':
    main()    