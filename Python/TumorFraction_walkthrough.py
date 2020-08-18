# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 11:32:38 2020

@author: amurtha
"""

import pandas as pd
import scipy.stats as stats
import numpy as np

pd.set_option("display.max_columns", None)

def get_adj_maf(maf, depth):
    alt = maf*depth
    for p in np.arange(0,1,0.005):
        dist = stats.binom(depth, p)
        if dist.cdf(alt) < 0.95:
            maf = p; break;
    return maf;

# =============================================================================
# Create mutation table
# =============================================================================

sample_id = ['17-095-3rd-cfDNA'] * 5
chrom = ['chr5', 'chr8', 'chr11', 'chr17', 'chr17']
gene = ['TERT', 'FGFR1', 'ATM', 'TP53', 'TP53']
af = [0.032, 0.074, 0.053, 0.07, 0.072] #allele frequency is equal to mutant reads / total reads
read_depth = [312, 1519, 1351, 1262, 1283]
effect = ['Upstream',
 'Missense',
 'Missense',
 'Synonymous',
 'Missense']

muts = pd.DataFrame({'Sample_ID':sample_id,'CHROM':chrom,'GENE':gene,'EFFECT':effect,'Allele_frequency':af, 'Read_depth':read_depth})

del sample_id, chrom,gene,af,read_depth,effect

# =============================================================================
# Create copy number table
# =============================================================================

sample_id = ['17-095-3rd-cfDNA'] * 4
gene = ['TERT', 'FGFR1', 'ATM', 'TP53']
copy_num = [0, 0, 0, 0]
lr = [0.01, 0.04, -0.06, -0.03]

cn = pd.DataFrame({'Sample_ID':sample_id,'GENE':gene,'Copy_num':copy_num,'Log_ratio':lr,})

del sample_id,gene,lr,copy_num

# =============================================================================
# Tumor fraction estimation
# =============================================================================

print('Sample tumor fraction (TF) calculation: 17-095-3rd-cfDNA')

print('\nSample mutations:')
print(muts.head())

print('\nCheck log-ratio of mutated genes. Do not use mutations on genes with log-ratio > 0.3. In this case, no genes amplified.')
print(cn)

print('\nCheck depth. Do not use mutation with depth < 30 reads. In this case, all depths are above 30.')
print(muts[['GENE','Read_depth']])

print('\nCheck for mutations on allosomes. Mutations on allosomes will not be considered when calculating tumor fraction. In this case, all mutations are on autosomes.')
print(muts[['GENE','CHROM']])

print('\nAdjust allele frequencies to conservatively calculate the VAF if the observed VAF is a 95% quantile outlier.')
muts['Adj_allele_frequency'] = muts.apply(lambda row: get_adj_maf(row['Allele_frequency'], row['Read_depth']), axis = 1)
print(muts[['GENE','Allele_frequency','Adj_allele_frequency']])

print('\nCalculate tumor fraction (TF) for every mutation using TF = 2 / (1 + 1 / AF), where AF is the adjusted allele frequency.')
muts['Tumor fraction'] = 2 / (1 + 1/muts['Adj_allele_frequency'])
print(muts[['GENE','Allele_frequency','Adj_allele_frequency','Tumor fraction']])

print('\nThe tumor fraction of the sample is the highest calculated. Mutation must not have low depth, an amplification, or be located on an allosome.')
print('\nIn this case, the tumor fraction of 17-095-3rd-cfDNA is 12.2%, based on the FGFR1 and TP53 mutations.')