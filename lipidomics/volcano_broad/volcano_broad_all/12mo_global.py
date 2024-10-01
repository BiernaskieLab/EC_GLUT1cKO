# # Libraries
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import numpy as np
import pandas as pd
from bioinfokit.analys import stat
import scipy.stats as stats
from scipy.stats import norm
import random

# # Formatting
# Seaborn (matplotlib)
sns.set_theme(style='ticks', font_scale=2)

# # File IO
    # Note - Easier to not change variable names so just look at files to see which dataset is being used!!

# load ventricle data
ventricle_c1 = pd.read_csv('162-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_c2 = pd.read_csv('163-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_c3 = pd.read_csv('164-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_c4 = pd.read_csv('565-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_ko1 = pd.read_csv('ventricle-271-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_ko2 = pd.read_csv('272-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_ko3 = pd.read_csv('292-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_ko4 = pd.read_csv('373-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
# load cortex data
cortex_c1 = pd.read_csv('162-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_c2 = pd.read_csv('163-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_c3 = pd.read_csv('164-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_c4 = pd.read_csv('565-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_ko1 = pd.read_csv('cortex-271-Total Ion Count.csv', delimiter=';', comment='#')
cortex_ko2 = pd.read_csv('272-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_ko3 = pd.read_csv('292-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_ko4 = pd.read_csv('373-cortex-Total Ion Count.csv', delimiter=';', comment='#')
# load caudate data
caudate_c1 = pd.read_csv('162-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_c2 = pd.read_csv('163-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_c3 = pd.read_csv('164-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_c4 = pd.read_csv('565-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_ko1 = pd.read_csv('caudate-271-Total Ion Count.csv', delimiter=';', comment='#')
caudate_ko2 = pd.read_csv('272-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_ko3 = pd.read_csv('292-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_ko4 = pd.read_csv('373-caudate-Total Ion Count.csv', delimiter=';', comment='#')

# load ventral data
ventral_c1 = pd.read_csv('162-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_c2 = pd.read_csv('163-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_c3 = pd.read_csv('164-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_c4 = pd.read_csv('565-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_ko1 = pd.read_csv('ventral-271-Total Ion Count.csv', delimiter=';', comment='#')
ventral_ko2 = pd.read_csv('272-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_ko3 = pd.read_csv('292-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_ko4 = pd.read_csv('373-ventral-Total Ion Count.csv', delimiter=';', comment='#')
# load corpus data
corpus_c1 = pd.read_csv('162-corpus-callosum-Total Ion Count.csv', delimiter=';', comment='#')
corpus_c2 = pd.read_csv('163-corpus-callosum-Total Ion Count.csv', delimiter=';', comment='#')
corpus_c3 = pd.read_csv('164-corpus-callosum-Total Ion Count.csv', delimiter=';', comment='#')
corpus_c4 = pd.read_csv('565-corpus-Total Ion Count.csv', delimiter=';', comment='#')
corpus_ko1 = pd.read_csv('corpus-271-Total Ion Count.csv', delimiter=';', comment='#')
corpus_ko2 = pd.read_csv('272-corpus-callosum-Total Ion Count.csv', delimiter=';', comment='#')
corpus_ko3 = pd.read_csv('292-corpus-callosum-Total Ion Count.csv', delimiter=';', comment='#')
corpus_ko4 = pd.read_csv('373-corpus-Total Ion Count.csv', delimiter=';', comment='#')
# load caudate+cortex data
cc_c1 = pd.read_csv('162-cortex+caudate-Total Ion Count.csv', delimiter=';', comment='#')
cc_c2 = pd.read_csv('163-cortex+caudate-Total Ion Count.csv', delimiter=';', comment='#')
cc_c3 = pd.read_csv('164-cortex+caudate-Total Ion Count.csv', delimiter=';', comment='#')
cc_c4 = pd.read_csv('565-cc-Total Ion Count.csv', delimiter=';', comment='#')
cc_ko1 = pd.read_csv('cc-271-Total Ion Count.csv', delimiter=';', comment='#')
cc_ko2 = pd.read_csv('272-cortex+caudate-Total Ion Count.csv', delimiter=';', comment='#')
cc_ko3 = pd.read_csv('292-cortex+caudate-Total Ion Count.csv', delimiter=';', comment='#')
cc_ko4 = pd.read_csv('373-cc-Total Ion Count.csv', delimiter=';', comment='#')

# Reference
ref = pd.read_csv('lipids-correlation-fullset-youngmice.csv', delimiter=',')

# # Global variables
# Set variables
ion_count = len(ventricle_c1.iloc[:,0])
size_set = 8
ident_set = ['c1', 'c2', 'c3', 'c4', 'ko1', 'ko2', 'ko3', 'ko4']
grp_annot = ref.iloc[:ion_count,8]
name_annot = ref.iloc[:ion_count,7]

# Ion means by brain
dist_ventricle = np.zeros([ion_count*2, size_set])
dist_ventricle_log2 = np.zeros([ion_count*2, size_set])

dist_cortex = np.zeros([ion_count*2, size_set])
dist_caudate = np.zeros([ion_count*2, size_set])
dist_cc = np.zeros([ion_count*2, size_set])
dist_corpus = np.zeros([ion_count*2, size_set])
dist_ventral = np.zeros([ion_count*2, size_set])
# Test statistic
stat_ventricle = np.zeros(100)
stat_cortex = np.zeros(100)
stat_caudate = np.zeros(100)
stat_cc = np.zeros(100)
stat_corpus = np.zeros(100)
stat_ventral = np.zeros(100)

# # Functions
def ret_dist(input):
# Arguments:    
    #input - dataframe

# Returns:
    # np.append(col_mean, col_std)
        # col_mean - np.ndarray containing mean for each ion intensity distribution
        # col_std - np.ndarray containing std. dev. for each ion intensity distribution

    col_mean = np.zeros(ion_count)
    col_std = np.zeros(ion_count)
    for i in range(ion_count):
        col_mean[i] = np.mean(input.iloc[i,:])
        col_std[i] = np.mean(input.iloc[i,:])
    return np.append(col_mean, col_std)


index_ref =  [0, ventricle_c1.shape[1], ventricle_c2.shape[1], ventricle_c3.shape[1], ventricle_c4.shape[1], \
              ventricle_ko1.shape[1], ventricle_ko2.shape[1], ventricle_ko3.shape[1], ventricle_ko4.shape[1]]

def ret_dist_log2(input, index):
# Arguments:    
    #input - dataframe
    #index - position in original dataframe

# Returns:
    # np.append(col_mean, col_std)
        # col_mean - np.ndarray containing mean for each ion intensity distribution
        # col_std - np.ndarray containing std. dev. for each ion intensity distribution

    init = 0
    for k in range(index+1):
        if k < index:
            init+=index_ref[k]
        if k == index:
            fin = init + index_ref[k]

    col_mean = np.zeros(ion_count)
    col_std = np.zeros(ion_count)
    for i in range(ion_count):
        col_mean[i] = np.mean(input.iloc[i,init:fin])
        col_std[i] = np.mean(input.iloc[i,init:fin])
    return np.append(col_mean, col_std)

def isolate_lipid_group(input, vector):
# Arguments:
    # input - dataframe
    # vector - group to sort

# Returns:
    # iso - np.ndarray, sorted means by lipid group
    # anno - np.ndarray, sorted annotations by lipid name within group

    j=[]
    for i in range(ion_count):
        if grp_annot[i] == vector:
            j.append(i)
    
    iso = np.zeros([len(j), size_set])
    anno = []

    for i in range(len(iso)):
        iso[i] = input[j[i],:]
        anno.append(name_annot[j[i]])

    return iso, anno

def diff_compare_local(reference, input):
# Arguments:
    # reference - Condition to compare against
    # input - Condition to compare with reference

# Returns:
    # Fold fractional difference relative to reference 
    
    output = []
    for i in range(np.shape(input)[0]):
        output.append((np.mean(input[i,:])-np.mean(reference[i,:])) / np.mean(reference[i,:]))
    return output

def log2FC(reference, input):
# Arguments:
    # reference - Condition to compare against
    # input - Condition to compare with reference

    # Note that this function requires log2scale input vectors

# Returns:
    # Log2FC

    output = []
    for i in range(np.shape(input)[0]):
        output.append((np.mean(input[i,:])-np.mean(reference[i,:])))
    return output


# # Main Program

# Ventricle
ventricle_all = [ventricle_c1, ventricle_c2, ventricle_c3, ventricle_c4, \
                 ventricle_ko1, ventricle_ko2, ventricle_ko3, ventricle_ko4]
vent_all = pd.concat([ventricle_c1, ventricle_c2.iloc[:,1:], ventricle_c3.iloc[:,1:], ventricle_c4.iloc[:,1:], \
                       ventricle_ko1.iloc[:,1:], ventricle_ko2.iloc[:,1:], ventricle_ko3.iloc[:,1:], \
                        ventricle_ko4.iloc[:,1:]], axis=1)

vent_all = vent_all.replace(0, np.nan)
imputed_all = vent_all.fillna(vent_all.iloc[:,1:].median())
log2_all = np.log2(imputed_all.iloc[:,1:])

for j in range(size_set):
    dist_ventricle[:,j] = ret_dist(ventricle_all[j])
    dist_ventricle_log2[:,j] = ret_dist_log2(log2_all, j+1)

# sort = dist_ventricle
sort = dist_ventricle_log2

ventricle_FA, annot_FA = isolate_lipid_group(sort, 'FA')
ventricle_ST, annot_ST = isolate_lipid_group(sort, 'ST')
ventricle_SP, annot_SP = isolate_lipid_group(sort, 'SP')
ventricle_GP, annot_GP = isolate_lipid_group(sort, 'GP')
ventricle_GL, annot_GL = isolate_lipid_group(sort, 'GL')

fatty_acyl = ventricle_FA
sterol = ventricle_ST
sphingolipid = ventricle_SP
glycerophospholipid = ventricle_GP
glycerolipid = ventricle_GL

# fa_diff = diff_compare_local(fatty_acyl[:,0:4], fatty_acyl[:,4:8])
# st_diff = diff_compare_local(sterol[:,0:4], sterol[:,4:8])
# sp_diff = diff_compare_local(sphingolipid[:,0:4], sphingolipid[:,4:8])
# gp_diff = diff_compare_local(glycerophospholipid[:,0:4], glycerophospholipid[:,4:8 ])
# gl_diff = diff_compare_local(glycerolipid[:,0:4], glycerolipid[:,4:8])

fa_diff = log2FC(fatty_acyl[:,0:4], fatty_acyl[:,4:8])
st_diff = log2FC(sterol[:,0:4], sterol[:,4:8])
sp_diff = log2FC(sphingolipid[:,0:4], sphingolipid[:,4:8])
gp_diff = log2FC(glycerophospholipid[:,0:4], glycerophospholipid[:,4:8 ])
gl_diff = log2FC(glycerolipid[:,0:4], glycerolipid[:,4:8])

fig, ax = plt.subplots(nrows=1, ncols=5)
ax[0].bar(x=annot_FA, height=fa_diff, label='fatty acyls', color='darkred')
ax[0].set(ylabel='Mean Fold Signal Differential')
ax[0].set(xlabel='fatty acyls')
ax[1].bar(x=annot_ST, height=st_diff, label='sterols', color='darkred')
ax[1].set(xlabel='sterols')
ax[2].bar(x=annot_SP, height=sp_diff, label='sphingolipids', color='darkred')
ax[2].set(xlabel='sphingolipids')
ax[3].bar(x=annot_GP, height=gp_diff, label='glycerophospholipids', color='darkred')
ax[3].set(xlabel='glycerophospholipids')
ax[4].bar(x=annot_GL, height=gl_diff, label='glycerolipids', color='darkred')
ax[4].set(xlabel='glycerolipids')


for k in range(len(ax)):
    ax[k].tick_params(axis='x', labelsize=2, labelrotation=80)
    ax[k].set(ylim=[-1.5,1])
    if k > 0 :
        ax[k].set(yticks=[])
fig.align_labels()

plt.show()

