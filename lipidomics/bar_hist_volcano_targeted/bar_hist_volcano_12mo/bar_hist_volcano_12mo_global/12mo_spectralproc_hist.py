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

# Data Split (Statistics)
def split(data_in, coeff):
    array_out = np.zeros(coeff)
    for j in range(coeff):
        array_out[j] = np.mean(random.sample(sorted(data_in), int(len(data_in)/coeff)))
        #array_out[j] = np.mean(data_in[int(j/coeff*len(data_in)):int(((j+1)/coeff)*len(data_in))])
    return array_out


# load ventricle data
ventricle_162 = pd.read_csv('162-ventricle.csv', delimiter=';')
ventricle_163 = pd.read_csv('163-ventricle.csv', delimiter=';')
ventricle_164 = pd.read_csv('164-ventricle.csv', delimiter=';')
ventricle_565 = pd.read_csv('565-ventricle.csv', delimiter=';')
ventricle_271 = pd.read_csv('271-ventricle.csv', delimiter=';')
ventricle_272 = pd.read_csv('272-ventricle.csv', delimiter=';')
ventricle_292 = pd.read_csv('292-ventricle.csv', delimiter=';')
ventricle_373 = pd.read_csv('373-ventricle.csv', delimiter=';')
# load cortex data
cortex_162 = pd.read_csv('162-cortex.csv', delimiter=';')
cortex_163 = pd.read_csv('163-cortex.csv', delimiter=';')
cortex_164 = pd.read_csv('164-cortex.csv', delimiter=';')
cortex_565 = pd.read_csv('565-cortex.csv', delimiter=';')
cortex_271 = pd.read_csv('271-cortex.csv', delimiter=';')
cortex_272 = pd.read_csv('272-cortex.csv', delimiter=';')
cortex_292 = pd.read_csv('292-cortex.csv', delimiter=';')
cortex_373 = pd.read_csv('373-cortex.csv', delimiter=';')
# load caudate data
caudate_162 = pd.read_csv('162-caudate.csv', delimiter=';')
caudate_163 = pd.read_csv('163-caudate.csv', delimiter=';')
caudate_164 = pd.read_csv('164-caudate.csv', delimiter=';')
caudate_565 = pd.read_csv('565-caudate.csv', delimiter=';')
caudate_271 = pd.read_csv('271-caudate.csv', delimiter=';')
caudate_272 = pd.read_csv('272-caudate.csv', delimiter=';')
caudate_292 = pd.read_csv('292-caudate.csv', delimiter=';')
caudate_373 = pd.read_csv('373-caudate.csv', delimiter=';')
# load ventral data
ventral_162 = pd.read_csv('162-ventral.csv', delimiter=';')
ventral_163 = pd.read_csv('163-ventral.csv', delimiter=';')
ventral_164 = pd.read_csv('164-ventral.csv', delimiter=';')
ventral_565 = pd.read_csv('565-ventral.csv', delimiter=';')
ventral_271 = pd.read_csv('271-ventral.csv', delimiter=';')
ventral_272 = pd.read_csv('272-ventral.csv', delimiter=';')
ventral_292 = pd.read_csv('292-ventral.csv', delimiter=';')
ventral_373 = pd.read_csv('373-ventral.csv', delimiter=';')
# load corpus data
corpus_162 = pd.read_csv('162-corpus.csv', delimiter=';')
corpus_163 = pd.read_csv('163-corpus.csv', delimiter=';')
corpus_164 = pd.read_csv('164-corpus.csv', delimiter=';')
corpus_565 = pd.read_csv('565-corpus.csv', delimiter=';')
corpus_271 = pd.read_csv('271-corpus.csv', delimiter=';')
corpus_272 = pd.read_csv('272-corpus.csv', delimiter=';')
corpus_292 = pd.read_csv('292-corpus.csv', delimiter=';')
corpus_373 = pd.read_csv('373-corpus.csv', delimiter=';')
# load caudate+cortex data
cc_162 = pd.read_csv('162-cc.csv', delimiter=';')
cc_163 = pd.read_csv('163-cc.csv', delimiter=';')
cc_164 = pd.read_csv('164-cc.csv', delimiter=';')
cc_565 = pd.read_csv('565-cc.csv', delimiter=';')
cc_271 = pd.read_csv('271-cc.csv', delimiter=';')
cc_272 = pd.read_csv('272-cc.csv', delimiter=';')
cc_292 = pd.read_csv('292-cc.csv', delimiter=';')
cc_373 = pd.read_csv('373-cc.csv', delimiter=';')

# -------------------------------------- FA 16:0 - iloc[0,:] ----------------------------------------------- #
i=0

# # Ventricle
# # Region Data
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 16:0 Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# k=2
# mean_162 = (split(ventricle_162.iloc[i,:], k))
# mean_163 = (split(ventricle_163.iloc[i,:], k))
# mean_164 = (split(ventricle_164.iloc[i,:], k))
# mean_565 = (split(ventricle_565.iloc[i,:], k))
# mean_271 = (split(ventricle_271.iloc[i,:], k))
# mean_272 = (split(ventricle_272.iloc[i,:], k))
# mean_292 = (split(ventricle_292.iloc[i,:], k))
# mean_373 = (split(ventricle_373.iloc[i,:], k))

# # mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# # mean_ko = np.concatenate([mean_272, mean_292])

# # mean_ctrl = np.concatenate([mean_565])
# # mean_ko = np.concatenate([mean_271, mean_373])

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
# mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

# def statistic(x,y,axis):
#     return np.mean(x, axis=axis)-np.mean(y, axis=axis)

# res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
# print(res.statistic, res.pvalue)
# u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

# #sns.color_palette("vlag", as_cmap=True)
# fig, ax = plt.subplots()
# label=['ko', 'ctrl']
# sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
# fig.suptitle("FA 16:0 - Ventricle")
# ax.set_xlabel('Groups')
# ax.set_ylabel('Normalized Ion Intensity')


plt.show()

# Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[0,:], cortex_163.iloc[0,:], cortex_164.iloc[0,:], cortex_565.iloc[0,:]])
# ko_concat = pd.concat([cortex_271.iloc[0,:], cortex_272.iloc[0,:], cortex_292.iloc[0,:], cortex_373.iloc[0,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 16:0 Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# k=2
# mean_162 = (split(cortex_162.iloc[i,:], k))
# mean_163 = (split(cortex_163.iloc[i,:], k))
# mean_164 = (split(cortex_164.iloc[i,:], k))
# mean_565 = (split(cortex_565.iloc[i,:], k))
# mean_271 = (split(cortex_271.iloc[i,:], k))
# mean_272 = (split(cortex_272.iloc[i,:], k))
# mean_292 = (split(cortex_292.iloc[i,:], k))
# mean_373 = (split(cortex_373.iloc[i,:], k))

# # mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# # mean_ko = np.concatenate([mean_272, mean_292])

# # mean_ctrl = np.concatenate([mean_565])
# # mean_ko = np.concatenate([mean_271, mean_373])

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
# mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

# def statistic(x,y,axis):
#     return np.mean(x, axis=axis)-np.mean(y, axis=axis)

# res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
# print(res.statistic, res.pvalue)
# u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

# #sns.color_palette("vlag", as_cmap=True)
# fig, ax = plt.subplots()
# label=['ko', 'ctrl']
# sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
# fig.suptitle("FA 16:0 - Cortex")
# ax.set_xlabel('Groups')
# ax.set_ylabel('Normalized Ion Intensity')
# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[0,:], caudate_163.iloc[0,:], caudate_164.iloc[0,:], caudate_565.iloc[0,:]])
# ko_concat = pd.concat([caudate_271.iloc[0,:], caudate_272.iloc[0,:], caudate_292.iloc[0,:], caudate_373.iloc[0,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 16:0 Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[0,:], cc_163.iloc[0,:], cc_164.iloc[0,:], cc_565.iloc[0,:]])
# ko_concat = pd.concat([cc_271.iloc[0,:], cc_272.iloc[0,:], cc_292.iloc[0,:], cc_373.iloc[0,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 16:0 Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[0,:], ventral_163.iloc[0,:], ventral_164.iloc[0,:], ventral_565.iloc[0,:]])
# ko_concat = pd.concat([ventral_271.iloc[0,:], ventral_272.iloc[0,:], ventral_292.iloc[0,:], ventral_373.iloc[0,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 16:0 Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[0,:], corpus_163.iloc[0,:], corpus_164.iloc[0,:], corpus_565.iloc[0,:]])
# ko_concat = pd.concat([corpus_271.iloc[0,:], corpus_272.iloc[0,:], corpus_292.iloc[0,:], corpus_373.iloc[0,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 16:0 Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()





# # ----------------------------------- FA 18:1 (iloc[1,:]) --------------------------------------- #
i=1
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 18:1 Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

k=2
mean_162 = (split(ventricle_162.iloc[i,:], k))
mean_163 = (split(ventricle_163.iloc[i,:], k))
mean_164 = (split(ventricle_164.iloc[i,:], k))
mean_565 = (split(ventricle_565.iloc[i,:], k))
mean_271 = (split(ventricle_271.iloc[i,:], k))
mean_272 = (split(ventricle_272.iloc[i,:], k))
mean_292 = (split(ventricle_292.iloc[i,:], k))
mean_373 = (split(ventricle_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
# mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

# def statistic(x,y,axis):
#     return np.mean(x, axis=axis)-np.mean(y, axis=axis)

# res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
# print(res.statistic, res.pvalue)
# u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

# #sns.color_palette("vlag", as_cmap=True)
# fig, ax = plt.subplots()
# label=['ko', 'ctrl']
# sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
# fig.suptitle("FA 18:1 - Ventricle")
# ax.set_xlabel('Groups')
# ax.set_ylabel('Normalized Ion Intensity')

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 18:1 Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 18:1 Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 18:1 Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 18:1 Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 18:1 Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()




# # ----------------------------------- FA 20:4, ST 20:1;O2 (iloc[1,:]) --------------------------------------- #
# i=4
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 20:4, ST 20:1;O2 Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 20:4, ST 20:1;O2 Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 20:4, ST 20:1;O2 Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 20:4, ST 20:1;O2 Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 20:4, ST 20:1;O2 Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('FA 20:4, ST 20:1;O2 Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()






# # ----------------------------------- NAE 15:1, SPB 17:2;O2 (iloc[1,:]) --------------------------------------- #
# i=2
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('NAE 15:1, SPB 17:2;O2 Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('NAE 15:1, SPB 17:2;O2 Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('NAE 15:1, SPB 17:2;O2 Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('NAE 15:1, SPB 17:2;O2 Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('NAE 15:1, SPB 17:2,O2 Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('NAE 15:1, SPB 17:2,O2 Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()





# # ----------------------------------- ST 21:5;O8;GlcA (iloc[1,:]) --------------------------------------- #
# i=3
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('ST 21:5;O8;GlcA Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('ST 21:5;O8;GlcA Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('ST 21:5;O8;GlcA Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('ST 21:5;O8;GlcA Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('ST 21:5;O8;GlcA Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('ST 21:5;O8;GlcA Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()





# # ----------------------------------- HexCer 43:6;O6 (iloc[1,:]) --------------------------------------- #
# i=7
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 43:6;O6 Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 43:6;O6 Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 43:6;O6 Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 43:6;O6 Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 43:6;O6 Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 43:6;O6 Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()





# # ----------------------------------- HexCer 41:2;O3 (iloc[1,:]) --------------------------------------- #
# i=6
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 41:2;O3 Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 41:2;O3 Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 41:2;O3 Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 41:2;O3 Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 41:2;O3 Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('HexCer 41:2;O3 Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()




# # ----------------------------------- CerPE 42:2;O6, SM 39:2;O6 (iloc[1,:]) --------------------------------------- #
# i=8
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 42:2;O6, SM 39:2;O6 Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 42:2;O6, SM 39:2;O6 Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 42:2;O6, SM 39:2;O6 Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 42:2;O6, SM 39:2;O6 Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 42:2;O6, SM 39:2;O6 Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 42:2;O6, SM 39:2;O6 Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()





# # ----------------------------------- CerPE 48:3,O5, ... (iloc[1,:]) --------------------------------------- #
# i=12
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 48:3,O5, ...  Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 48:3,O5, ... Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 48:3,O5, ... Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 48:3,O5, ... Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 48:3,O5, ... Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('CerPE 48:3,O5, ... Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()



# # ----------------------------------- PI 36:4, PG 38:5;O (iloc[1,:]) --------------------------------------- #
# i=5
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 36:4, PG 38:5;O Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 36:4, PG 38:5;O Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 36:4, PG 38:5;O Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 36:4, PG 38:5;O Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 36:4, PG 38:5;O Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 36:4, PG 38:5;O Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()






# # ----------------------------------- PI 38:4, PG 40:5;O (iloc[1,:]) --------------------------------------- #
# i=9
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 38:4, PG 40:5;O Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 38:4, PG 40:5;O Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 38:4, PG 40:5;O Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 38:4, PG 40:5;O Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 36:4, PG 38:5;O Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PI 38:4, PG 40:5;O Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()





# # ----------------------------------- PC O-43:11, ... (iloc[1,:]) --------------------------------------- #
# i=10
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PC O-43:11, ... Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PC O-43:11, ... Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PC O-43:11, ... Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PC O-43:11, ... Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PC O-43:11, ... Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PC O-43:11, ... Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()





# # ----------------------------------- PG 42:5, PG O-42:6;O (iloc[1,:]) --------------------------------------- #
# i=11
# # Ventricle
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventricle_162.iloc[i,:], ventricle_163.iloc[i,:], ventricle_164.iloc[i,:], ventricle_565.iloc[i,:]])
# ko_concat = pd.concat([ventricle_271.iloc[i,:], ventricle_272.iloc[i,:], ventricle_292.iloc[i,:], ventricle_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PG 42:5, PG O-42:6;O Ventricular Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cortex_162.iloc[i,:], cortex_163.iloc[i,:], cortex_164.iloc[i,:], cortex_565.iloc[i,:]])
# ko_concat = pd.concat([cortex_271.iloc[i,:], cortex_272.iloc[i,:], cortex_292.iloc[i,:], cortex_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PG 42:5, PG O-42:6;O Cortex Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([caudate_162.iloc[i,:], caudate_163.iloc[i,:], caudate_164.iloc[i,:], caudate_565.iloc[i,:]])
# ko_concat = pd.concat([caudate_271.iloc[i,:], caudate_272.iloc[i,:], caudate_292.iloc[i,:], caudate_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PG 42:5, PG O-42:6;O Caudate Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Cortex+Caudate
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([cc_162.iloc[i,:], cc_163.iloc[i,:], cc_164.iloc[i,:], cc_565.iloc[i,:]])
# ko_concat = pd.concat([cc_271.iloc[i,:], cc_272.iloc[i,:], cc_292.iloc[i,:], cc_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PG 42:5, PG O-42:6;O Cortex+Caudate Regions')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Ventral
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([ventral_162.iloc[i,:], ventral_163.iloc[i,:], ventral_164.iloc[i,:], ventral_565.iloc[i,:]])
# ko_concat = pd.concat([ventral_271.iloc[i,:], ventral_272.iloc[i,:], ventral_292.iloc[0,:], ventral_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PG 42:5, PG O-42:6;O Ventral Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

# plt.show()

# # Corpus Callosum
# # concatenate ctrl, ko data
# ctrl_concat = pd.concat([corpus_162.iloc[i,:], corpus_163.iloc[i,:], corpus_164.iloc[i,:], corpus_565.iloc[i,:]])
# ko_concat = pd.concat([corpus_271.iloc[i,:], corpus_272.iloc[i,:], corpus_292.iloc[i,:], corpus_373.iloc[i,:]])
# ctrl_concat = ctrl_concat.reset_index(drop=True)
# ko_concat = ko_concat.reset_index(drop=True)
# data = { 'ctrl': ctrl_concat.iloc[1:], 'ko': ko_concat.iloc[1:]}

# coupled_data = pd.DataFrame(data)
# ctrl = pd.DataFrame(coupled_data['ctrl'])
# ctrl = ctrl.loc[(ctrl!=0).any(axis=1)]
# ko = pd.DataFrame(coupled_data['ko'])
# ko = ko.loc[(ko!=0).any(axis=1)]

# # Plot Distributions
# sns.set_theme(style="ticks", font_scale=2)
# fig, (ax_box1, ax_box2, ax_hist) = plt.subplots(3, sharex=True, gridspec_kw={"height_ratios": (0.1, 0.1, 0.8)})
# sns.boxplot(ctrl, ax=ax_box1, orient="h", color='slategrey')
# sns.boxplot(ko, ax=ax_box2, orient="h", color='darkred')
# sns.kdeplot(ctrl['ctrl'], ax=ax_hist, color='slategrey', bw_adjust=0.04, fill=True, label='control')
# sns.kdeplot(ko['ko'], ax=ax_hist, color='darkred', bw_adjust=0.04, fill=True, label='knockout')

# ax_box1.set(yticks=[])
# ax_box2.set(yticks=[])
# ax_hist.set(xlabel='Ion Intensity (Total Ion Count Normalization)')
# ax_hist.set(ylabel='Spatial Density')
# ax_hist.legend()
# sns.despine(ax=ax_hist)
# sns.despine(ax=ax_box1, left=True)
# sns.despine(ax=ax_box2, left=True)
# fig.suptitle('PG 42:5, PG O-42:6;O Corpus Callosum Region')

# # Statistics
# print('Std. Deviation of x (ctrl):', np.std(data['ctrl']))
# print('Std. Deviation of y (ko):', np.std(data['ko']))
# res = stat()
# res.ztest(df=coupled_data, x='ctrl', y='ko', x_std=np.std(data['ctrl']), y_std=np.std(data['ko']), test_type=2)
# print(res.summary)
# u_stat, u_p = stats.mannwhitneyu(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])
# w_stat, w_p = stats.ranksums(x=data['ctrl'], y=data['ko'], alternative='greater')
# print('Non-parametric Wilcoxon Rank-Sum Test', [w_stat, w_p])

plt.show()