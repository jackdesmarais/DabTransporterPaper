#/usr/bin/python

"""Common definition of colors and plotting utilities."""

import seaborn as sns

import matplotlib
from matplotlib import pyplot as plt

import matplotlib.font_manager

matplotlib.font_manager.findSystemFonts(fontpaths='/System/Library/Fonts/*', fontext='ttf')

qual_palette = sns.color_palette('Paired', n_colors=10, desat=0.9)

font_dict = {'family':'sans-serif','sans-serif':['Helvetica']}
plt.rc('font', **font_dict)
plt.rc('axes', labelsize=26) 
