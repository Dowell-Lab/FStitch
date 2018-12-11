#!/usr/bin/env python

from __future__ import print_function
import argparse
import numpy as np
import datetime
import pandas as pd
import chartify
import sys
import os
from pybedtools import BedTool
from argparse import RawTextHelpFormatter

def main():

    parser = argparse.ArgumentParser(description='FStitch Bidirectional Module\n\n===============================================================================================================\n\nAnnotate bidirectionals using FStitch segment output.', epilog='@Dowell Lab, Margaret Gruca, margaret.gruca@colorado.edu\nFor questions and issues, see https://github.com/Dowell-Lab/FStitch', usage='%(prog)s --bed fstitch.bed --genes gene_ref.bed --output /my/out/file.bed', formatter_class=RawTextHelpFormatter)
    
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    required.add_argument('-b', '--bed', dest='fstitch_seg_filename', metavar='FSTITCH.BED', \
                       help='FStitch segment output (BED file), concatenated for both postive and negative strands.', required=True)
    
    required.add_argument('-g', '--genes', dest='gene_ref', metavar='GENE_REF.BED', \
                       help='Gene reference file in BED format.', required=True)
    
    required.add_argument('-o', '--output', dest='output', metavar='SAMPLE.BED', \
                        help='Directory where output BED file, stats, and plots will be saved. Full path and file extension for the BED file must be specified.', required=True)
    
    optional.add_argument('-tss', '--removetss', dest='tss_remove', action='store_true', \
                        help='Remove bidirectionals that fall directly over a gene transcription start site. Default = False', default=False, required=False)
    
    optional.add_argument('-s', '--split', dest='split', action='store_true', \
                        help='Split files into short and long bidirectionals. Default = False', default=False, required=False)
    
    optional.add_argument('-f', '--footprint', dest='footprint', metavar='<FOOTPRINT>', \
                       help='The footprint is a gap between positive and negative reads. This function will add an integer value (in bp) to merge positive and negative segments that do not overlap. This value should likely be increased for lower complexity data and will have minimal effect for high complexity data. Default = 300.', default=300, required=False)
    
    optional.add_argument('-lg', '--mergelength', dest='merge_length', metavar='<MERGE_LENGTH>', \
                        help='Choose length (in bp) for short/long merge length. Short and long calls are segregated and merged separately to prevent short calls from being merged into long bidirectional regions (e.g. superenhancers, unanoated genes/lncRNAs). Default=12000', default=12000, required=False)
    
    optional.add_argument('-lm', '--maxlength', dest='bidir_length', metavar='<BIDIR_LENGTH>', \
                        help='Integer value (in bp) for max bidirectional length. Default=25000', default=25000, required=False)
    
    optional.add_argument('-ls', '--splitlength', dest='split_length', metavar='<SPLIT_LENGTH>', \
                        help='Choose length (in bp) for short/long bidirectional file split. Only an option if -s flag is specified. Default=8000', default=8000, required=False)
    
    optional.add_argument('-p', '--plotbidirs', dest='gen_plot', action='store_true', \
                        help='Generate a histogram plot for bidirectional lengths. Default = False', default=False, required=False)
    
    args = parser.parse_args()
    

    if (int(args.bidir_length) <= int(args.merge_length)):
            print ("Bidirectional max length cannot be less than the merge length.")
            sys.exit()
    
    elif (int(args.bidir_length) <= int(args.split_length)):
            print("Bidirectional max length cannot be less than the split length.")
            sys.exit()
    else:
        print('Max bidirectional length set to:', args.bidir_length,'bp')
        print('Footprint set to:', args.footprint,'bp')
        print('Split length set to:', args.split_length,'bp')
        print(str(datetime.datetime.now()) + '\nStarting bidirectional caller.....')
    
    fstitch_seg_file = pd.read_table((args.fstitch_seg_filename), header=None, skiprows=1, usecols=range(6))
    genes = open(args.gene_ref)
    rootname = os.path.splitext(args.output)[0]

### Format FStitch file

    print('Parsing strand information.....')
    
    fstitch_seg_file[3], fstitch_seg_file[6] = fstitch_seg_file[3].str.split('=', 1).str
    fstitch_seg_file.columns = ['chromosome', 'start', 'end', 'activity', 'igv_tag', 'strand', 'CI']
    fstitch_seg_file = fstitch_seg_file[~fstitch_seg_file['activity'].isin(['OFF'])]
    fstitch_seg_file = fstitch_seg_file.drop(columns= ['activity' , 'igv_tag' , 'CI'])
    
### Parse the strands and 300bp to end of negative strand -- this acts as a pseduo bedtools '-d' flag or a 'footprint' value 
    
# Getting a warning here... believe it is a false positive but will triple check later. For now, we'll supress the warning output
    pd.options.mode.chained_assignment = None
    
    fs_neg = fstitch_seg_file[fstitch_seg_file['strand'] != "+"]
    fs_neg['end'] = fs_neg.end + int(args.footprint)
    
    fs_pos = fstitch_seg_file[fstitch_seg_file['strand'] != "-"]
    
### Format the gene reference file
    
    print('Parsing gene reference information.....')
    
    genes = pd.read_table((args.gene_ref), header=None, usecols=range(6))
    genes.columns = ['chromosome', 'start', 'end', 'id' , 'identifier', 'strand']
    genes = genes.drop(columns= ['identifier'])
    genes_pos = genes[genes['strand'] != "-"]
    genes_neg = genes[genes['strand'] != "+"]
    
    if(args.tss_remove):
        print('Parsing transcription start sites.....')
        genes_pos_tss = genes_pos
        genes_pos_tss.end = genes_pos.start + 1
        genes_neg_tss = genes_neg
        genes_neg_tss.start = genes_neg.end - 1
        
        genes_tss = pd.concat([genes_neg_tss, genes_pos_tss])
                        
### Use pybedtools to filter intragenic bidirections using segregated fstitch and gene data
# First create 'BedTool' objects from dataframes
                        
    bt_genes = BedTool.from_dataframe(genes)
    bt_genes_pos = BedTool.from_dataframe(genes_pos)
    bt_genes_neg = BedTool.from_dataframe(genes_neg)
    bt_fs_pos = BedTool.from_dataframe(fs_pos)
    bt_fs_neg = BedTool.from_dataframe(fs_neg)
    bt_nogenes_fs_pos = bt_fs_pos.subtract(bt_genes)
    bt_nogenes_fs_neg = bt_fs_neg.subtract(bt_genes)
    
    print('Detecting intragenic bidirectionals.....')
     
# Intersect segements on opposite strand genes                    
                        
    pos_bidirs = bt_genes_pos.intersect(bt_fs_neg, wb=True)
    pos_bidirs = pos_bidirs.sort()
    pos_bidirs = pos_bidirs.merge()
                        
    neg_bidirs = bt_genes_neg.intersect(bt_fs_pos, wb=True)
    neg_bidirs = neg_bidirs.sort()
    neg_bidirs = neg_bidirs.merge()
    
# Get rid of duplicates over isoforms and expand regions
    
    pos_bidirs = BedTool.to_dataframe(pos_bidirs)
    pos_bidirs = pos_bidirs.drop_duplicates(subset=['start', 'end'])
    pos_bidirs['end'] = pos_bidirs['end'] + (pos_bidirs['end'] - pos_bidirs['start'])
    
    neg_bidirs = BedTool.to_dataframe(neg_bidirs)
    neg_bidirs = neg_bidirs.drop_duplicates(subset=['start', 'end'])
    neg_bidirs['start'] = neg_bidirs['start'] - (neg_bidirs['end'] - neg_bidirs['start'])
    
    intragenic_bidirs = pd.concat([pos_bidirs, neg_bidirs])
                          
    print('Detecting intergenic bidirectionals.....')
                          
### Get all other bidirectionals (intergenic)
                          
    bd1 = bt_fs_neg.intersect(bt_fs_pos, wo=True)
    bd1 = BedTool.to_dataframe(bd1)
    bd1 = bd1.drop(columns = ['name', 'score', 'strand', 'thickEnd', 'itemRgb', 'end'])
    bd1.columns = ['chrom', 'start', 'end']
    
    bd2 = bt_fs_neg.intersect(bt_fs_pos)
    bd2 = BedTool.to_dataframe(bd2)
    bd2 = bd2.drop(columns = ['name'])
    
    bd3 = bt_nogenes_fs_neg.intersect(bt_nogenes_fs_pos, wo=True)
    bd3 = BedTool.to_dataframe(bd3)
    bd3 = bd3.drop(columns = ['name', 'score', 'strand', 'thickEnd', 'itemRgb', 'end'])
    bd3.columns = ['chrom', 'start', 'end']
    
    concat_intergenic_bidirs = pd.concat([bd1, bd2, bd3])
                          
    bt_concat_intergenic_bidirs = BedTool.from_dataframe(concat_intergenic_bidirs)
    bt_intergenic_bidirs = bt_concat_intergenic_bidirs.subtract(bt_genes, A=True)
    intergenic_bidirs = BedTool.to_dataframe(bt_intergenic_bidirs)
    
### Now that we have both intragenic and intergenic bidirectionals, it's time to filter and merge based on size                      
                          
    df = pd.concat([intergenic_bidirs, intragenic_bidirs])
    
# Optionally remove bidirectionals that fall directly over a tss
    
    if(args.tss_remove):
        bt_genes_tss = BedTool.from_dataframe(genes_tss)
        bt_df = BedTool.from_dataframe(df)
        bt_df = bt_df.subtract(bt_genes_tss, A=True)
        df = BedTool.to_dataframe(bt_df)
    
    df['diff'] = df.end - df.start
    
    print('Filtering bidirectionals by length.....')                      
    
# This is controled by parameter 'l' currently. Should the following filters be based on this parameter, as well...?
    
    dropped_bidirs_short = df[(df['diff'] <= 100)]
    dropped_bidirs_long = df[(df['diff'] >= int(args.bidir_length))]
    df = df[df['diff'] <= int(args.bidir_length)]
    
# These steps merge based on size. Do not want to merge large things with small things... this ends up mering all discrete calls with gene/lcnRNA/superenhancer calls
    
    df1 = df[(df['diff'] <= int(args.merge_length)) & (df['diff'] >= 100)]
    bt_df1 = BedTool.from_dataframe(df1)
    bt_df1 = bt_df1.sort().merge()
    df1 = BedTool.to_dataframe(bt_df1)
    
    df2 = df[df['diff'] > int(args.merge_length)]
    bt_df2 = BedTool.from_dataframe(df2)
    bt_df2 = bt_df2.sort().merge()
    df2 = BedTool.to_dataframe(bt_df2)
                          
    dff = pd.concat([df1, df2])
    bt_dff = BedTool.from_dataframe(dff)
    bt_dff = bt_dff.sort()
    dff = BedTool.to_dataframe(bt_dff)
    dff['id'] = dff.index + 1
    dff['id'] = 'bidir_' + dff['id'].astype(str)
    dff['length'] = dff.end - dff.start
    dff.to_csv((args.output), sep="\t", header=None, index=False)
    
# Save options : split file based on a length (-ls)
# Print stats
    
    if args.split:
        print('Splitting output into short and long bidirectionals...')
        
        dff_short = dff[dff['length'] <= int(args.split_length)]
        dff_short.to_csv((rootname + '.short.bed'), sep="\t", header=None, index=False)
        dff_long = dff[dff['length'] > int(args.split_length)]
        dff_long.to_csv((rootname + '.long.bed'), sep="\t", header=None, index=False)
    
        stat_name = ['footprint', 'max_length', 'split_length' 'fstitch_pos_segs', 'fstitch_neg_segs', 'dropped_short_bidirs', 'dropped_long_bidirs', 'intragenic_bidirs', 'intergenic_bidirs', 'mean_length', 'median_length', 'total_bidirs_short', 'total_bidirs_long', 'total_bidirs']
        stat_value = [(args.footprint), (args.bidir_length), (args.split_length), len(fs_pos.index), len(fs_neg.index), len(dropped_bidirs_short.index), len(dropped_bidirs_long.index), len(intragenic_bidirs.index), len(intergenic_bidirs.index),  round(dff['length'].mean()), round(dff['length'].median()), len(dff_short.index), len(dff_long.index), len(dff.index)]
            
        stats = pd.DataFrame([stat_name, stat_value])
        stats.to_csv((rootname + '.stats.txt'), sep='\t', header=None, index=False)
    
    else:
    
        stat_name = ['footprint', 'max_length', 'fstitch_pos_segs', 'fstitch_neg_segs', 'dropped_short_bidirs', 'dropped_long_bidirs', 'intragenic_bidirs', 'intergenic_bidirs', 'mean_length', 'median_length', 'total_bidirs']
        stat_value = [(args.footprint), (args.bidir_length), len(fs_pos.index), len(fs_neg.index), len(dropped_bidirs_short.index), len(dropped_bidirs_long.index), len(intragenic_bidirs.index), len(intergenic_bidirs.index),  round(dff['length'].mean()), round(dff['length'].median()), len(dff.index)]
            
        stats = pd.DataFrame([stat_name, stat_value])
        stats.to_csv((rootname + '.stats.txt'), sep='\t', header=None, index=False)
                          
    if args.gen_plot:
        print('Generating a histogram for bidirectional lengths...')
        print('Mean length: ', dff['length'].mean())
        print('Median length: ', dff['length'].median())
        
        (chartify.Chart(blank_labels=True, y_axis_type='density')
        .plot.histogram(
           data_frame=dff,
           values_column='length',
           bins=50)
        .set_title('Bidirectional Lengths')
        .set_subtitle('All Bidirectionals')
        .axes.set_xaxis_label('Length')
        .axes.set_yaxis_label('Number')
        .save(rootname + '.length_hist.html'))
    
        dff_short = dff[dff['length'] <= int(args.split_length)]
        print('Adjusted mean length: ', dff_short['length'].mean())
        print('Adjusted median length: ', dff_short['length'].median())
        
        (chartify.Chart(blank_labels=True, y_axis_type='density')
        .plot.histogram(
           data_frame=dff_short,
           values_column='length',
           bins=50)
        .set_title('Bidirectional Lengths')
        .set_subtitle('Length < ' + str(args.split_length) + 'bp')
        .axes.set_xaxis_label('Length')
        .axes.set_yaxis_label('Number')
        .save(rootname + '.short_length_hist.html'))
                          
    print('Bidirectional module complete.\n' + str(datetime.datetime.now()))
    sys.exit(0)


if __name__=='__main__':
    main()