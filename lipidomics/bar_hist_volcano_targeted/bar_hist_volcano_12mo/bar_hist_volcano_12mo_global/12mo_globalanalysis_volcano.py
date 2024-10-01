__author__ = "James Colter"
__copyright__ = "Open Source"
__credits__ = ["James Colter"]
__license__ = "GPL"
__version__ = "1"
__maintainer__ = "James Colter"
__email__ = "jdcolter@ucalgary.ca"
__status__ = "Research"

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import numpy as np
import pandas as pd
from bioinfokit.analys import stat
import scipy.stats as stats
from scipy.stats import norm
import random

sns.set_theme(style='ticks', font_scale=2)

pval = pd.read_csv(r'pvalue.csv')
mko = pd.read_csv(r'meanko.csv')
mctrl = pd.read_csv(r'meanCTRL.csv')

fold = pd.concat([mko, mctrl.iloc[:,1:]], axis=1)
imputed_fold = fold.fillna(fold.iloc[:,1:].median())
log2_fold = np.log2(imputed_fold.iloc[:,1:])
print(log2_fold)

# Mean Differentials
diff = mko.iloc[:,1:].div(mctrl.iloc[:,1:])
log2FC = log2_fold.iloc[:,:6] - log2_fold.iloc[:,6:12]
print(log2FC)

# Plots

# # Ventricle
# fig, ax = plt.subplots(nrows=1, ncols=4)
# fig.suptitle('Ventricular Downregulation by Lipid Group')
# sns.barplot(-0.01/pval.iloc[0:4,1], ax=ax[0], color='darkred')
# ax[0].set_xlabel('Fatty Acyls')
# sns.barplot(-0.01/pval.iloc[4:5,1], ax=ax[1], color='maroon')
# ax[1].set_xlabel('Sterols')
# sns.barplot(-0.01/pval.iloc[5:9,1], ax=ax[2], color='dimgray')
# ax[2].set_xlabel('Sphingolipids')
# sns.barplot(-0.01/pval.iloc[9:13,1], ax=ax[3], color='lightgray')
# ax[3].set_xlabel('Glycerophospholipids')
# plt.setp(ax, ylim=[-1,0], ylabel=None)
# ax[0].set_ylabel('Relative Downregulation')

# # Cortex
# i=2
# fig, ax = plt.subplots(nrows=1, ncols=4)
# fig.suptitle('Cortex Downregulation by Lipid Group')
# sns.barplot(-0.005/pval.iloc[0:4,i], ax=ax[0], color='darkred')
# ax[0].set_xlabel('Fatty Acyls')
# sns.barplot(-0.01/pval.iloc[4:5,i], ax=ax[1], color='maroon')
# ax[1].set_xlabel('Sterols')
# sns.barplot(-0.01/pval.iloc[5:9,i], ax=ax[2], color='dimgray')
# ax[2].set_xlabel('Sphingolipids')
# sns.barplot(-0.01/pval.iloc[9:13,i], ax=ax[3], color='lightgray')
# ax[3].set_xlabel('Glycerophospholipids')
# plt.setp(ax, ylim=[-1,0], ylabel=None)
# ax[0].set_ylabel('Relative Downregulation')


# # Caudate
# i=3
# fig, ax = plt.subplots(nrows=1, ncols=4)
# fig.suptitle('Caudate Downregulation by Lipid Group')
# sns.barplot(-0.01/pval.iloc[0:4,i], ax=ax[0], color='darkred')
# ax[0].set_xlabel('Fatty Acyls')
# sns.barplot(-0.01/pval.iloc[4:5,i], ax=ax[1], color='maroon')
# ax[1].set_xlabel('Sterols')
# sns.barplot(-0.01/pval.iloc[5:9,i], ax=ax[2], color='dimgray')
# ax[2].set_xlabel('Sphingolipids')
# sns.barplot(-0.01/pval.iloc[9:13,i], ax=ax[3], color='lightgray')
# ax[3].set_xlabel('Glycerophospholipids')
# plt.setp(ax, ylim=[-1,0], ylabel=None)
# ax[0].set_ylabel('Relative Downregulation')

# # Cortex + Caudate
# i=4
# fig, ax = plt.subplots(nrows=1, ncols=4)
# fig.suptitle('Cortex+Caudate Downregulation by Lipid Group')
# sns.barplot(-0.01/pval.iloc[0:4,i], ax=ax[0], color='darkred')
# ax[0].set_xlabel('Fatty Acyls')
# sns.barplot(-0.01/pval.iloc[4:5,i], ax=ax[1], color='maroon')
# ax[1].set_xlabel('Sterols')
# sns.barplot(-0.01/pval.iloc[5:9,i], ax=ax[2], color='dimgray')
# ax[2].set_xlabel('Sphingolipids')
# sns.barplot(-0.01/pval.iloc[9:13,i], ax=ax[3], color='lightgray')
# ax[3].set_xlabel('Glycerophospholipids')
# plt.setp(ax, ylim=[-1,0], ylabel=None)
# ax[0].set_ylabel('Relative Downregulation')

# # Corpus Callosum
# i=5
# fig, ax = plt.subplots(nrows=1, ncols=4)
# fig.suptitle('Corpus Callosum Downregulation by Lipid Group')
# sns.barplot(-0.01/pval.iloc[0:4,i], ax=ax[0], color='darkred')
# ax[0].set_xlabel('Fatty Acyls')
# sns.barplot(-0.01/pval.iloc[4:5,i], ax=ax[1], color='maroon')
# ax[1].set_xlabel('Sterols')
# sns.barplot(-0.01/pval.iloc[5:9,i], ax=ax[2], color='dimgray')
# ax[2].set_xlabel('Sphingolipids')
# sns.barplot(-0.01/pval.iloc[9:13,i], ax=ax[3], color='lightgray')
# ax[3].set_xlabel('Glycerophospholipids')
# plt.setp(ax, ylim=[-1,0], ylabel=None)
# ax[0].set_ylabel('Relative Downregulation')

# # Ventral
# i=6
# fig, ax = plt.subplots(nrows=1, ncols=4)
# fig.suptitle('Ventral Downregulation by Lipid Group')
# sns.barplot(-0.01/pval.iloc[0:4,i], ax=ax[0], color='darkred')
# ax[0].set_xlabel('Fatty Acyls')
# sns.barplot(-0.01/pval.iloc[4:5,i], ax=ax[1], color='maroon')
# ax[1].set_xlabel('Sterols')
# sns.barplot(-0.01/pval.iloc[5:9,i], ax=ax[2], color='dimgray')
# ax[2].set_xlabel('Sphingolipids')
# sns.barplot(-0.01/pval.iloc[9:13,i], ax=ax[3], color='lightgray')
# ax[3].set_xlabel('Glycerophospholipids')
# plt.setp(ax, ylim=[-1,0], ylabel=None)
# ax[0].set_ylabel('Relative Downregulation')


# # Volcano Plot (Diff)
# for i in range(0,6):
#     fig, ax = plt.subplots()
#     plt.scatter(x=-1/diff.iloc[0:4,i], y=-np.log10(pval.iloc[0:4,i+1]), color='darkred')
#     plt.scatter(x=-1/diff.iloc[4:5,i], y=-np.log10(pval.iloc[4:5,i+1]), color='salmon')
#     plt.scatter(x=-1/diff.iloc[5:9,i], y=-np.log10(pval.iloc[5:9,i+1]), color='dimgray')
#     plt.scatter(x=-1/diff.iloc[9:13,i], y=-np.log10(pval.iloc[9:13,i+1]), color='lightgray')
#     fig.suptitle('Year-one Mice')
#     ax.set_title(list(diff)[i])
#     ax.set_xlabel('Fold Signal Change')
#     ax.set_ylabel('-log10(p-value)')
#     plt.axvline(-0.5,color="grey",linestyle="--")
#     plt.axhline(1.3,color="grey",linestyle="--")
#     ax.set_xlim([-4,0.1])
#     ax.set_ylim([0,3])
#     plt.tight_layout()

# Volcano Plot (Log2FC)
for i in range(0,6):
    fig, ax = plt.subplots()
    plt.scatter(x=log2FC.iloc[0:4,i], y=-np.log10(pval.iloc[0:4,i+1]), color='darkred')
    plt.scatter(x=log2FC.iloc[4:5,i], y=-np.log10(pval.iloc[4:5,i+1]), color='salmon')
    plt.scatter(x=log2FC.iloc[5:9,i], y=-np.log10(pval.iloc[5:9,i+1]), color='dimgray')
    plt.scatter(x=log2FC.iloc[9:13,i], y=-np.log10(pval.iloc[9:13,i+1]), color='lightgray')
    fig.suptitle('Month-One Mice')
    ax.set_title(list(diff)[i])
    ax.set_xlabel('Fold Signal Change')
    ax.set_ylabel('-log10(p-value)')
    plt.axvline(-0.25,color="grey",linestyle="--")
    plt.axhline(1.3,color="grey",linestyle="--")
    ax.set_xlim([-4,0.1])
    ax.set_ylim([0,3])
    plt.tight_layout()
plt.show()