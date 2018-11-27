#/usr/bin/python

"""Common definition of colors and plotting utilities."""

import pandas as pd
import seaborn as sns
import matplotlib
import numpy as np


from Bio import SeqIO
from os import path

from scipy.stats import ttest_ind

from matplotlib import pyplot as plt

import matplotlib.font_manager

matplotlib.font_manager.findSystemFonts(fontpaths='/System/Library/Fonts/*', fontext='ttf')

greyish = sns.xkcd_rgb['greyish']

qual_palette = sns.color_palette('Paired', n_colors=10, desat=0.9)

font_dict = {'family':'sans-serif','sans-serif':['Helvetica','DejaVu Sans']}
plt.rc('font', **font_dict)
plt.rc('axes', labelsize=26) 


# Genome and pool data filenames
genome_fname = "../data/berkeleylab-feba-7de9bcaf0960/bin/g/Halo/genome.fna"
pool_fname1 = "../data/HnTnSeqAnalysis/DSJD02HTn1_R1/DSJDo2Htn1_S1_l006_R1_001.pool"
pool_fname2 = "../data/HnTnSeqAnalysis/DSJD03Htn2_R1/DSJD03Htn2_S2_l006_R1_001.pool"

gene_stats_fname1 = '../data/HnTnSeqAnalysis/geneStatsRep1.csv'
gene_stats_fname2 = '../data/HnTnSeqAnalysis/geneStatsRep2.csv'

# DataFrames for pool and gene data
pool_df1 = pd.read_csv(pool_fname1, sep='\t')
pool_df2 = pd.read_csv(pool_fname2, sep='\t')
rep_df1 = pd.read_csv(gene_stats_fname1)
rep_df2 = pd.read_csv(gene_stats_fname2)

# DataFrame for essentiality calls 
essentiality_fname = '../data/essentiality_calls.csv'
essentiality_df = pd.read_csv(essentiality_fname)

# Grab genome length
genome_size_bp = None
for req in SeqIO.parse(genome_fname, "fasta"):
    genome_size_bp = len(req)

# When we are plotting the library we want to include insertions mapped in either replicate.
# Therefore we merge the two dataframes. 
merge_on = pool_df1.columns.tolist()
total_pool_df = pool_df1.merge(pool_df2, on=merge_on, how='outer')
total_pool_df['n_total'] = total_pool_df.n + total_pool_df.n2


# loading more data from for plotting
pathToData = "../data/2017_11_03_miseq_barSeq1/"
logratios_df = pd.read_csv(path.join(pathToData, 'fit_logratios.tab'), sep='\t')

# smaller_df
small_log_df = pd.DataFrame()
small_log_df["locusId"] = logratios_df["locusId"]
small_log_df["sysName"] = logratios_df["sysName"]
small_log_df["desc"] = logratios_df["desc"]
small_log_df["HighCO2Log"] = logratios_df["setAS2 HCO"]
small_log_df["LowCO2Log"] = logratios_df["setAS3 LCO"]

pathToData = "../data/2017_11_03_miseq_barSeq2/"
logratios_df = pd.read_csv(path.join(pathToData, 'fit_logratios.tab'), sep='\t')

# smaller_df
small_log_df2 = pd.DataFrame()
small_log_df2["locusId"] = logratios_df["locusId"]
small_log_df2["sysName"] = logratios_df["sysName"]
small_log_df2["desc"] = logratios_df["desc"]
small_log_df2["HighCO2Log"] = logratios_df["setAS2 HCO"]
small_log_df2["LowCO2Log"] = logratios_df["setAS3 LCO"]

merge_on = ['locusId','sysName','desc']

small_log_df = small_log_df.merge(small_log_df2, on = merge_on, how = 'outer', suffixes=('_1','_2'))

# Absolute difference of 1.0 in the log fitness indicates twofold difference in linear scale 
# since we are using log base2 here. 
diff1 = small_log_df.LowCO2Log_1 - small_log_df.HighCO2Log_1
diff2 = small_log_df.LowCO2Log_2 - small_log_df.HighCO2Log_2

small_log_df['hcr'] = (diff1 <= -1.0) & (diff2 <= -1.0)

big_diff_mask =  (np.abs(diff1) >= 1.0) & (np.abs(diff2) >= 1.0)
small_diff_mask = np.logical_not(big_diff_mask)
big_diff_df = small_log_df[big_diff_mask]
small_diff_df = small_log_df[small_diff_mask]