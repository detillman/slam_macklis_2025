#!/usr/bin/env python3
import pandas as pd
import numpy as np
import pysam
import sys
from itertools import product
import re
import time



bam_filename = sys.argv[1]
vcf_filename = sys.argv[2]
output_file = sys.argv[3]


def annotate_read_mismatch_indel(read, vcf):

    letter_lst = ['A', 'T', 'G', 'C', 'N']
    letter_dict = dict(zip(letter_lst, range(len(letter_lst))))
    conversions = np.zeros((5, 5))
    
    # populate conversions by iterating over read
    # get tuples corresponding to the (location query, location reference, reference base)
    align_tuples = read.get_aligned_pairs(matches_only=True, with_seq=True)

    # get the read query sequence
    query_seq = read.query_sequence #read.query_alignment_sequence

    snp_pos = [variant.pos for variant in vcf.fetch(read.reference_name, read.reference_start, read.reference_end)]

    if len(query_seq) < len(align_tuples):
            print("Fail")
            print(align_tuples)
            print(len(align_tuples))
            print(query_seq)
            print(len(query_seq))

    for tup in align_tuples:
       ref = tup[2].upper()
       query = query_seq[tup[0]].upper()

       # check VCF
       if tup[1] not in snp_pos:
          ref_idx =  letter_dict[ref]
          query_idx = letter_dict[query]
          conversions[ref_idx, query_idx] += 1

          global master_conversions
          master_conversions[ref_idx, query_idx] += 1
       else:
          global snps_pos_masked
          snps_pos_masked += 1

    for (ref_idx, query_idx) in product(range(len(letter_lst)), repeat=2):
        read.set_tag(
            "{}{}".format(letter_lst[ref_idx], letter_lst[query_idx]).lower(), 
            conversions[ref_idx, query_idx], 
            value_type="i")

    # also annotate indel
    if re.search("I|D", read.cigarstring):
       read.set_tag("di", 1, value_type="i")
    else:
        read.set_tag("di", 0, value_type="i")

    return read


start = time.time()
in_bam = pysam.AlignmentFile(bam_filename, "rb")
out_bam = pysam.AlignmentFile(output_file, "wb", template=in_bam)

vcf = pysam.VariantFile(vcf_filename)

master_conversions = np.zeros((5, 5))

snps_pos_masked = 0


for chrom in [str(x+1) for x in range(19)] + ["X", "Y"]: #in_bam.header.references:
    sys.stdout.write(chrom)

    for read in in_bam.fetch(chrom):
       out_bam.write(annotate_read_mismatch_indel(read, vcf))

    
    

out_bam.close()
in_bam.close()

sys.stdout.write("Global conversion rate: rows=reference, columns=read sequence")
sys.stdout.write(pd.DataFrame(master_conversions, 
    columns=['A', 'T', 'G', 'C', 'N'], 
    index=['A', 'T', 'G', 'C', 'N']))
sys.stdout.write("Number of masked SNP positions: {}".format(snps_pos_masked))

end = time.time()
sys.stdout.write("Elapsed time: {} seconds".format(int(end-start)))
