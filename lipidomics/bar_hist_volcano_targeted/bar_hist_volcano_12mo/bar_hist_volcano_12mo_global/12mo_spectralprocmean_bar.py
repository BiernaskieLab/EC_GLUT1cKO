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

plt.figure()
plt.imshow(ventricle_272.iloc[1:,1:])

# Random Sampling (k = # sub-population groups)
# sample_size =
k=2

# Variables (Hold data)
all_means_ctrl = np.zeros((13,6))
all_means_ko = np.zeros((13,6))
all_pvalue = np.zeros((13,6))

# Data Split (Statistics)
def split(data_in, coeff):
    array_out = np.zeros(coeff)
    for j in range(coeff):
        array_out[j] = np.mean(random.sample(population=sorted(data_in), k=int(len(data_in)/coeff)))
        #array_out[j] = np.mean(data_in[int(j/coeff*len(data_in)):int(((j+1)/coeff)*len(data_in))])
    return array_out

# -------------------------------------- FA 16:0 - iloc[0,:] ----------------------------------------------- #
i=0

# ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

print('\n')
all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 16:0 - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 16:0 - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 16:0 - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 16:0 - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 16:0 - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 16:0 - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)

plt.tight_layout()
plt.show()


# ----------------------------------- FA 18:1 (iloc[1,:]) --------------------------------------- #
i=1
# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 18:1 - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 18:1 - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 18:1 - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 18:1 - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

print('\n')
all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 18:1 - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')
print('\n')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 18:1 - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()




# # ----------------------------------- FA 20:4, ST 20:1;O2 (iloc[1,:]) --------------------------------------- #
i=4
# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 20:4 - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 20:4 - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 20:4 - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 20:4 - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 20:4 - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("FA 20:4 - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()

# # ----------------------------------- NAE 15:1, SPB 17:2;O2 (iloc[1,:]) --------------------------------------- #
i=2
# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("NAE 15:1, SPB 17:2;O2 - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("NAE 15:1, SPB 17:2;O2 - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("NAE 15:1, SPB 17:2;O2 - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("NAE 15:1, SPB 17:2;O2 - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("NAE 15:1, SPB 17:2;O2 - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("NAE 15:1, SPB 17:2;O2 - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()





# # ----------------------------------- ST 21:5;O8;GlcA (iloc[1,:]) --------------------------------------- #
i=3

# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("ST 21:5;O8;GlcA - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("ST 21:5;O8;GlcA - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("ST 21:5;O8;GlcA - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("ST 21:5;O8;GlcA - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("ST 21:5;O8;GlcA - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("ST 21:5;O8;GlcA - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()



# # ----------------------------------- HexCer 43:6;O6 (iloc[1,:]) --------------------------------------- #
i=7

# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 43:6;O6 - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 43:6;O6 - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 43:6;O6 - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 43:6;O6 - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 43:6;O6 - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 43:6;O6 - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()





# # ----------------------------------- HexCer 41:2;O3 (iloc[1,:]) --------------------------------------- #
i=6
# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 41:2;O3 - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 41:2;O3 - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 41:2;O3 - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 41:2;O3 - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 41:2;O3 - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("HexCer 41:2;O3 - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()



# # ----------------------------------- CerPE 42:2;O6, SM 39:2;O6 (iloc[1,:]) --------------------------------------- #
i=8

# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 42:2;O6, SM 39:2;O6 - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 42:2;O6, SM 39:2;O6 - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 42:2;O6, SM 39:2;O6 - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 42:2;O6, SM 39:2;O6 - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 42:2;O6, SM 39:2;O6 - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 42:2;O6, SM 39:2;O6 - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()




# # ----------------------------------- CerPE 48:3,O5, ... (iloc[1,:]) --------------------------------------- #
i=12

# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 48:3,O5, ... - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 48:3,O5, ...- Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 48:3,O5, ... - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 48:3,O5, ... - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 48:3,O5, ... - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("CerPE 48:3,O5, ... - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()




# # ----------------------------------- PI 36:4, PG 38:5;O (iloc[1,:]) --------------------------------------- #
i=5

# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 36:4, PG 38:5;O - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 36:4, PG 38:5;O - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 36:4, PG 38:5;O - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 36:4, PG 38:5;O - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 36:4, PG 38:5;O - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 36:4, PG 38:5;O - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()




# # ----------------------------------- PI 38:4, PG 40:5;O (iloc[1,:]) --------------------------------------- #
i=9

# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 38:4, PG 40:5;O - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 38:4, PG 40:5;O - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 38:4, PG 40:5;O - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 38:4, PG 40:5;O - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 38:4, PG 40:5;O - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PI 38:4, PG 40:5;O - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()




# # ----------------------------------- PC O-43:11, ... (iloc[1,:]) --------------------------------------- #
i=10

# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PC O-43:11, ...- Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PC O-43:11, ... - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PC O-43:11, ... - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PC O-43:11, ... - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PC O-43:11, ... - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PC O-43:11, ... - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()




# # ----------------------------------- PG 42:5, PG O-42:6;O (iloc[1,:]) --------------------------------------- #
i=11

# # Ventricle
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

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PG 42:5, PG O-42:6;O - Ventricle")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,0] = res.pvalue
all_means_ctrl[i,0] = np.mean(mean_ctrl)
all_means_ko[i,0] = np.mean(mean_ko)
plt.tight_layout()

# Cortex
mean_162 = (split(cortex_162.iloc[i,:], k))
mean_163 = (split(cortex_163.iloc[i,:], k))
mean_164 = (split(cortex_164.iloc[i,:], k))
mean_565 = (split(cortex_565.iloc[i,:], k))
mean_271 = (split(cortex_271.iloc[i,:], k))
mean_272 = (split(cortex_272.iloc[i,:], k))
mean_292 = (split(cortex_292.iloc[i,:], k))
mean_373 = (split(cortex_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PG 42:5, PG O-42:6;O - Cortex")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,1] = res.pvalue
all_means_ctrl[i,1] = np.mean(mean_ctrl)
all_means_ko[i,1] = np.mean(mean_ko)
plt.tight_layout()

# # Caudate
mean_162 = (split(caudate_162.iloc[i,:], k))
mean_163 = (split(caudate_163.iloc[i,:], k))
mean_164 = (split(caudate_164.iloc[i,:], k))
mean_565 = (split(caudate_565.iloc[i,:], k))
mean_271 = (split(caudate_271.iloc[i,:], k))
mean_272 = (split(caudate_272.iloc[i,:], k))
mean_292 = (split(caudate_292.iloc[i,:], k))
mean_373 = (split(caudate_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PG 42:5, PG O-42:6;O - Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,2] = res.pvalue
all_means_ctrl[i,2] = np.mean(mean_ctrl)
all_means_ko[i,2] = np.mean(mean_ko)
plt.tight_layout()

# # Cortex+Caudate
mean_162 = (split(cc_162.iloc[i,:], k))
mean_163 = (split(cc_163.iloc[i,:], k))
mean_164 = (split(cc_164.iloc[i,:], k))
mean_565 = (split(cc_565.iloc[i,:], k))
mean_271 = (split(cc_271.iloc[i,:], k))
mean_272 = (split(cc_272.iloc[i,:], k))
mean_292 = (split(cc_292.iloc[i,:], k))
mean_373 = (split(cc_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PG 42:5, PG O-42:6;O - Cortex+Caudate")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,3] = res.pvalue
all_means_ctrl[i,3] = np.mean(mean_ctrl)
all_means_ko[i,3] = np.mean(mean_ko)
plt.tight_layout()

# # Corpus Callosum
mean_162 = (split(corpus_162.iloc[i,:], k))
mean_163 = (split(corpus_163.iloc[i,:], k))
mean_164 = (split(corpus_164.iloc[i,:], k))
mean_565 = (split(corpus_565.iloc[i,:], k))
mean_271 = (split(corpus_271.iloc[i,:], k))
mean_272 = (split(corpus_272.iloc[i,:], k))
mean_292 = (split(corpus_292.iloc[i,:], k))
mean_373 = (split(corpus_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PG 42:5, PG O-42:6;O - Corpus Callosum")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,4] = res.pvalue
all_means_ctrl[i,4] = np.mean(mean_ctrl)
all_means_ko[i,4] = np.mean(mean_ko)
plt.tight_layout()

#ventral
mean_162 = (split(ventral_162.iloc[i,:], k))
mean_163 = (split(ventral_163.iloc[i,:], k))
mean_164 = (split(ventral_164.iloc[i,:], k))
mean_565 = (split(ventral_565.iloc[i,:], k))
mean_271 = (split(ventral_271.iloc[i,:], k))
mean_272 = (split(ventral_272.iloc[i,:], k))
mean_292 = (split(ventral_292.iloc[i,:], k))
mean_373 = (split(ventral_373.iloc[i,:], k))

# mean_ctrl = np.concatenate([mean_162, mean_163, mean_164])
# mean_ko = np.concatenate([mean_272, mean_292])

# mean_ctrl = np.concatenate([mean_565])
# mean_ko = np.concatenate([mean_271, mean_373])

mean_ctrl = np.concatenate([mean_162, mean_163, mean_164, mean_565])
mean_ko = np.concatenate([mean_272, mean_292, mean_271, mean_373])

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)

res = stats.permutation_test((mean_ko, mean_ctrl), statistic, batch=len(mean_ctrl)-1, vectorized=True, permutation_type='independent', n_resamples=1e4, alternative='less')
print(res.statistic, res.pvalue)
u_stat, u_p = stats.mannwhitneyu(x=mean_ctrl, y=mean_ko, alternative='greater')
print('Non-parametric Mann-Whitney U Test', [u_stat, u_p])

#sns.color_palette("vlag", as_cmap=True)
fig, ax = plt.subplots()
label=['ko', 'ctrl']
sns.boxplot(data = {'ctrl': mean_ctrl, 'ko': mean_ko}, ax=ax, palette='coolwarm')
fig.suptitle("PG 42:5, PG O-42:6;O - Ventral")
ax.set_xlabel('Groups')
ax.set_ylabel('Normalized Ion Intensity')

all_pvalue[i,5] = res.pvalue
all_means_ctrl[i,5] = np.mean(mean_ctrl)
all_means_ko[i,5] = np.mean(mean_ko)
plt.tight_layout()

plt.show()

print(all_means_ctrl)
print(all_means_ko)
print(all_pvalue)