#!/usr/bin/env python3
import pandas as pd
import numpy as np
import pysam
import sys
import re
import time
#from multiprocessing import Pool, Process, Lock
import os
import collections
import datetime
import argparse
import csv

bam_filename = sys.argv[1]
gtf_filename = sys.argv[2]
vcf_filename = sys.argv[3]
output_prefix = sys.argv[4]

# Parse commandline arguments
# parser = argparse.ArgumentParser(description='This is python implementation of TimeLapse mutation calling')
# requiredNamed = parser.add_argument_group('required named arguments')
# requiredNamed.add_argument('-b', '--bam', type=str, required=True, metavar = 'in_file.bam',
#                     help='Bam file to process')
# requiredNamed.add_argument("-g", "--gtf", type=str, required=True, metavar='reference.gtf', help="GTF with exon annotation")
# parser.add_argument('--mutType', default='TC', type=str, choices=['TC', 'GA', 'TC,GA'],
#                     help='Type of mutation to record (default: TC)')
# parser.add_argument('--reads', default='PE', type=str, choices=['PE', 'SE'],
#                     help='Type of mutation to record (default: PE)')
# parser.add_argument('--minDist', default=5, type=int, metavar = '<int>',
#                     help='Base distance from read-end filter (default: 5)')
# parser.add_argument('--minQual', default=40, type=int, metavar = '<int>',
#                     help='Base minimal quality filter (default: 40)')
# parser.add_argument('--tracks', action='store_true', # Automatically stored as default = FALSE
#                     help='Generate files necessary for creating browser tracks')
# parser.add_argument('--mutsRDS', action='store_true',
#                     help='Generate _muts.rds like file with mutation frequency output')
# parser.add_argument('--mutPos', action='store_true',
#                     help='Generate _cU.rds like file and mutation position bedGraph files')
# parser.add_argument('--SNPs', help='Path to snp.txt file')
# parser.add_argument('--strandedness', default='F', type=str, choices=['F', 'R'],
#                     help='Is first read forward or reverse orientation? F = forward is the default.')
# args = parser.parse_args()

# args.mutType = args.mutType.split(',')
# args.base = [x[0] for x in args.mutType]        # base nucleotide: TC => T
# inputName = args.bam.split('.bam')[0]           # name without .bam suffix

# if args.strandedness == 'F':
#     strand_check = True
# else:
#     strand_check = False




# # Load SNPs for filtering
# snp = {}
# snpFile = open(args.SNPs, 'r')
# for line in snpFile:
#     line = line.strip().split(':')
#     snp[line[2] + ':' + line[3]] = line[0] + ':' + line[1]


class GTF:
    def __init__(self, gtf_filename):
        self.gene_exon_dict = self.populate_gene_info(gtf_filename)
        self.gene_region_dict = self.merge_gene_intervals()
        self.filename = gtf_filename

    def populate_gene_info(self, gtf_filename):
        gene_exon_dict = collections.defaultdict(list)
        with open(gtf_filename, "r") as in_gtf:
            j = 0
            for i in in_gtf.readlines():
                if "#!" in i:
                    continue
                chrom, _, feature, start, end, _, strand, _, info = i.strip("\n").split("\t")
                start = int(start)
                end = int(end)
                if feature == "exon":
                    gene_id = info.split(";")["gene_id" in info.split(";")].split()[1].strip("\"")
                    gene_exon_dict[gene_id].append((chrom, start, end, strand))
                j += 1
                if j % 100000 == 0:
                    print("Processed {} lines".format(j))
        return gene_exon_dict

    def merge_gene_intervals(self):
        gene_region_dict = collections.defaultdict(list)
        for gene in self.gene_exon_dict.keys():
            gene_region_dict[gene] = self.merge_intervals(gene)

        return gene_region_dict

    def merge_intervals(self, geneid):

        # get the intervals corresponding to a geneid
        # order by starting position
        all_exons = list(set(self.gene_exon_dict[geneid]))
        all_exons.sort(key=lambda i:i[1])

        merged_intervals = []
        chrom, curr_start, curr_end, strand = all_exons[0]
        for (_, new_start, new_end, _) in all_exons[1:]:
            if new_start < curr_end: 
                curr_end = max(new_end, curr_end)
            else:
                merged_intervals.append((chrom, curr_start, curr_end, strand))
                curr_start = new_start
                curr_end = new_end

        merged_intervals.append((chrom, curr_start, curr_end, strand))

        return merged_intervals
    

    def get_gene_regions(self, geneid):

        return self.gene_region_dict[geneid]

    def get_genes(self):
        return list(self.gene_region_dict.keys())


class MutationCounter:
    def __init__(self) :
        self.muts={'TA': 0, 'CA': 0, 'GA': 0, 'NA': 0, 
                   'AT': 0, 'CT': 0, 'GT': 0, 'NT': 0, 
                   'AC': 0, 'TC': 0, 'GC': 0, 'NC': 0, 
                   'AG': 0, 'TG': 0, 'CG': 0, 'NG': 0, 
                   'AN': 0, 'TN': 0, 'CN': 0, 'GN': 0,
                   'TT': 0, 'AA': 0, 'CC': 0, 'GG': 0
        }
        
    def add(self, read_loc_tuple):
        
       # read_loc_tuple: (ref_base, query_base, align_quality, orientation)
       base_info = [read_loc_tuple.ref_base, read_loc_tuple.query_base]
       if read_loc_tuple.orientation == "F":
           mut = "".join(base_info)
       else:
           mut = "".join([get_revcomp(x) for x in base_info])
           
       if mut in self.muts.keys():
           self.muts[mut] += 1
    
       else:
           raise ValueError("Invalid mutation type")
           
    def get_muts(self, muts_lst):
        """
        

        Parameters
        ----------
        muts_lst : list of mutations.

        Returns
        -------
        sum of all mutations in list

        """
        
        return sum([self.muts[x] for x in muts_lst])
    
# define the read alignment named tuple
readAlignmentInfo = collections.namedtuple('readAlignmentInfo', 
                                           ['ref_base', 'query_base', 
                                            'align_quality', 'orientation'])
          

def get_revcomp(base):
    """
    

    Parameters
    ----------
    base : string DNA nucleotide base (upper or lower).

    Returns
    -------
    The reverse complement of that base. Preserves upper/lower case.

    """

    return {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N', 
            'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n'}[base] 

def get_read_alignment_info(r, strandedness):
    """
    

    Parameters
    ----------
    r : pysam.AlignedSegment
    strandedness: string "F" or "R", which determines whether read 1 of library is F or R

    Returns
    -------
    Dictionary with { genomic_position: (ref_base, query_base, align_quality)}

    """

    
    read_alignment_info = r.get_aligned_pairs(matches_only = True, with_seq=True)
    ref_pos = [x[1] for x in read_alignment_info]
    ref_base = [x[2] for x in read_alignment_info]
    
    read_info = zip(ref_pos, ref_base,
                                           list(r.query_alignment_sequence), 
                                           r.query_alignment_qualities)
                                           
                                           

    # figure out gene orientation using libtype and read fwd/rev info
    # if an odd number of libtype == "F", r.is_forward, and r.is_read1, 
    # then gene is in forward orientation
    # else gene is in reverse orientation
    if sum([strandedness == "F", r.is_read1, not r.is_reverse]) % 2 == 1:
        orientation = "F"
    else:
        orientation = "R"
  
    return { x[0] : readAlignmentInfo(x[1].upper(), x[2].upper(), x[3], orientation) for x in read_info }
    

def get_fragment_conversions(readlst, snps, strandedness, minqual=30):
    """
    

    Parameters
    ----------
    readlst : list of length 2 [read_1, read_2] if PE and [read_1] if SE 
    snps: list of genomic positions with SNPs 
    strandedness : string "F" or "R" designating which strand is read 1
    minqual: minimum base quality to call a mutation

    Returns
    -------
    MutationCounter with mutations for that region with SNP masking
    use the read with the higher quality base call
    
    Previous (confusing) specs
    1) If there is no mutation in 1st and in 2nd read => replace with higher quality 2nd read
    2) If there is mutation in 1st and in 2nd read => replace with higher quality 2nd read
    3) If there is mutation in 1nd but not in 2st read => replace only if 1st read quality is less than threshold and less than 2nd read
    4) If there is mutation in 2nd read but not in 1st read => replace even with lower quality 2nd read as long as 2nd read quality is higher than threshold

    """
    
    
    r1 = readlst[0]
    mcounter = MutationCounter()
    
    # get location and read info
    #  { genomic_position: (ref_base, query_base, align_quality, gene orientation)}
    r1_align = get_read_alignment_info(r1, strandedness)
    
    if len(readlst) == 1:
        
        for loc in r1_align.keys():
            if not loc in snps and r1_align[loc].align_quality >= minqual:
                mcounter.add(r1_align[loc])
                
    # otherwise its paired end and both pairs exist
    else:
        r2 = readlst[1]
        r2_align = get_read_alignment_info(r2, strandedness)
    
        # get set of overlapping genomic loci
        dovetail = set(r1_align.keys()) & set(r2_align.keys())
        
        # if there's overlap, add mutations based on highest quality base
        for loc in dovetail:
            # make sure not SNP location
            if loc in snps:
                continue
            
            # if quality score as good or better for read1 use that
            elif r1_align[loc].align_quality >= r2_align[loc].align_quality: 
                if r1_align[loc].align_quality >= minqual:
                    mcounter.add(r1_align[loc])
            # else use read 2
            else:
                if r2_align[loc].align_quality >= minqual:
                    mcounter.add(r2_align[loc])
                     
        # add the remaining genomic positions:
        for loc in set(r1_align.keys())-set(dovetail):
            if r1_align[loc].align_quality >= minqual:
                mcounter.add(r1_align[loc])
                
        for loc in set(r2_align.keys())-set(dovetail):
            if r2_align[loc].align_quality >= minqual:
                mcounter.add(r2_align[loc])
                
                
    return mcounter
    
               
# TODO: probably want to turn this into a generator
def get_conversions_region(region, bamfile, vcf, libtype, strandedness, muttype):
    """
    

    Parameters
    ----------
    region : TYPE
        DESCRIPTION.
    bamfile : TYPE
        DESCRIPTION.
    vcf : TYPE
        DESCRIPTION.
    libtype : TYPE
        DESCRIPTION.
    muttype: string with first element reference base and second query base of mutation to quantify

    Returns
    -------
    None.

    """
    
    snp_pos = [variant.pos for variant in vcf.fetch(region[0],
                                                    region[1], 
                                                    region[2])]
    
    ref_base = muttype[0]
    mut_col = muttype
    muts_to_sum = ["{}{}".format(ref_base, y) for y in ["A", "T", "G", "C", "N"]]
    
    all_convs = pd.DataFrame(columns=[muttype, 'nref', 'ai', 'io', 'ei', 'sj', 'nfrag'])

    
    for readlst in get_reads_region(region, bamfile, libtype):
        
        gf = readlst[0].get_tag('GF')
        xf = readlst[0].get_tag('XF')
        
        # check if intron anywhere (0 if false, 1 if true for 1 of pair, 2 if true for both)
        ai = 1 if sum([r.get_tag('XF') == '__no_feature' and r.get_tag('GF') != '__no_feature' for r in readlst]) > 0 else 0
        
        # num of intron only
        io = 1 if sum([r.get_tag('EF') == '__no_feature' and r.get_tag('GF') != '__no_feature' for r in readlst]) == len(readlst) else 0
        
        # spans exon intron
        ei = 1 if sum([r.get_tag('XF') == '__no_feature' and r.get_tag('EF') != '__no_feature' for r in readlst]) > 0 else 0
        
        # spans exon-exon junction
        sj = 1 if sum([1 if re.match('N', r.cigarstring) is not None else 0 for r in readlst]) > 0 else 0
        

        # filter out reads with insertions/deletions
        readlst_filt = [r for r in filter(lambda r: not bool(re.match('I|D', r.cigarstring)),
                              readlst)]
        if len(readlst_filt) > 0:
                                   
            convs = get_fragment_conversions(readlst_filt, snp_pos, strandedness)
            all_convs = pd.concat([all_convs, 
                                  pd.Series(data={muttype: convs.get_muts([muttype]),
                                    'nref' : convs.get_muts(muts_to_sum),
                                    'ai'   : ai,
                                    'io'   : io,
                                    'ei'   : ei,
                                    'sj'   : sj,
                                    'nfrag': 1 }).to_frame().T],
                                    axis=0).groupby([muttype, 'nref', 'ai', 'io', 'ei', 'sj']).agg('sum').reset_index()

    return all_convs
        
        
        
    

def get_reads_region(region, bamfile, libtype):
    
    if libtype == "SE":
        for r in bamfile.fetch(region[0], region[1], region[2]):
            yield [r]
            
    elif libtype == "PE":
        read_dict = {}
        for r in bamfile.fetch(region[0], region[1], region[2]):
            
            # operate under the assumption that reads are uniquely mapped
            # so once a read name is observed twice, return both mates and delete from dict
            # otherwise more likely to result in segmentation fault 
            if not r.query_name in read_dict.keys(): # might be slow, could opt for trie structure
               read_dict[r.query_name] = r
            else:
               mate = read_dict[r.query_name]
               del read_dict[r.query_name]
               if mate.is_read1:
                  yield [mate, r]
               else:
                  yield [r, mate]
        # for reads with mates that are unmapped, still report
        for rname in read_dict.keys():
            
            yield [read_dict[rname]]  
    


def get_conversions_gene(bam, gtf, geneid, vcf, libtype, strandedness, muttype):
    """
    bam: bamfile
    gtf: GTF object defining genes
    geneid: gene identifier (or 4th column of bedfile)
    vcf: VCF object 
    libtype: "PE" or "SE"
    strandedness: "F" or "R" designating which orientation is read_1 on txpt
    muttype: string of length 2 with first char as ref_base and second as mut_base (assuming forward orientation on transcript)
    """
    regions = gtf.get_gene_regions(geneid)
    
    gene_convs = []
    for region in regions:
        gene_convs.append(get_conversions_region(region, bam, vcf, libtype, strandedness, muttype))
        
    # remove any empty conversions
    #gene_convs = [g for g in gene_convs if g is not None]

    agg_gene_convs = pd.concat(gene_convs, axis=0).groupby([muttype, 'nref', 'ai', 'io', 'ei', 'sj']).agg('sum').reset_index()
    agg_gene_convs['GF'] = geneid
    return agg_gene_convs

    
    




    
#### Main ####
if __name__ == "__main__": 
    start = time.time()

    sys.stdout.write("Bam filename: {}\n".format(bam_filename))
    sys.stdout.write("GTF filename: {}\n".format(gtf_filename))
    sys.stdout.write("Outprefix: {}\n".format(output_prefix))

    gtf = GTF(gtf_filename)
    in_bam = pysam.AlignmentFile(bam_filename, "rb")
    vcf = pysam.VariantFile(vcf_filename)

    strandedness = "F"
    libtype = "PE"
    muttype = "TC"
    
    all_genes = gtf.get_genes()
    sys.stdout.write("Number of genes to process: {}".format(str(len(all_genes))))

                             
    colnames = ["GF", muttype, "nref", "ai", "io", "ei", "sj", "nfrag"]                         
  
    out_filename="{}_counts.csv".format(output_prefix)
    with open(out_filename, "w") as outfile:
        outfile.write("{}\n".format(",".join(colnames)))

    g = 0
    for gene in all_genes:
            
        conversion_tbl = get_conversions_gene(in_bam, gtf, gene, vcf, libtype, strandedness, muttype)
        if len(conversion_tbl.index) > 0:
            conversion_tbl = conversion_tbl[colnames]
            conversion_tbl.to_csv(path_or_buf=out_filename, header=False, index=False, mode='a')
        g += 1
        if g % 100 == 0:
            sys.stdout.write("{} genes processed".format(str(g)))

    end = time.time()
    sys.stdout.write("Done. Elapsed time: {} seconds".format(int(end-start)))





#print(get_conversions_gene(in_bam, gtf, "ENSMUSG00000021025", max_num_conversions=5))
#with open('{}_region_gtf_compressed.pickle'.format(output_prefix), 'wb') as handle:
#    pickle.dump(gtf, handle, protocol=pickle.HIGHEST_PROTOCOL)

#with open('{}_region_gtf_compressed.pickle'.format(output_prefix), 'rb') as handle:
#    b = pickle.load(handle)

#
#chrom  = "4" #"12"
#start  = 119135178 #55536195
#end    = 119151801 #55539432
#strand = "-"

#tcount_convs_df, all_convs_df, tcount_total, tc_count_total, bp_count_total, bp_conv_count_total, invalid_reads = get_conversions_region(in_bam, chrom, start, end, strand)
#print(tcount_convs_df)
#print(all_convs_df)
#print("Tcount total: {}".format(tcount_total))
#print("TC count total: {}".format(tc_count_total))
#print("bp count total: {}".format(bp_count_total))
#print("bp conversion count total: {}".format(bp_conv_count_total))
#print("number of invalid reads: {}".format(invalid_reads))
