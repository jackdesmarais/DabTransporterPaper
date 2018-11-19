#/usr/bin/python

"""Common definition of colors and plotting utilities."""

import seaborn as sns

import matplotlib

qual_palette = sns.color_palette('Paired', n_colors=10, desat=0.9)

matplotlib.rcParams['font.sans-serif'] = ['Helvetica']+matplotlib.rcParams['font.sans-serif']
matplotlib.rcParams['font.sans-serif']
