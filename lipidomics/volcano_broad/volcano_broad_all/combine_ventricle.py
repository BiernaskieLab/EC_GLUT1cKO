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

# load ventricle data
ventricle_c1y = pd.read_csv('286-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_c2y = pd.read_csv('287-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_c3y = pd.read_csv('288-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_ko1y = pd.read_csv('509-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_ko2y = pd.read_csv('521-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
ventricle_ko3y = pd.read_csv('522-ventricle-Total Ion Count.csv', delimiter=';', comment='#')
# load cortex data
cortex_c1y = pd.read_csv('286-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_c2y = pd.read_csv('287-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_c3y = pd.read_csv('288-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_ko1y = pd.read_csv('509-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_ko2y = pd.read_csv('521-cortex-Total Ion Count.csv', delimiter=';', comment='#')
cortex_ko3y = pd.read_csv('522-cortex-Total Ion Count.csv', delimiter=';', comment='#')
# load caudate data
caudate_c1y = pd.read_csv('286-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_c2y = pd.read_csv('287-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_c3y = pd.read_csv('288-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_ko1y = pd.read_csv('509-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_ko2y = pd.read_csv('521-caudate-Total Ion Count.csv', delimiter=';', comment='#')
caudate_ko3y = pd.read_csv('522-caudate-Total Ion Count.csv', delimiter=';', comment='#')
# load ventral data
ventral_c1y = pd.read_csv('286-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_c2y = pd.read_csv('287-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_c3y = pd.read_csv('288-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_ko1y = pd.read_csv('509-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_ko2y = pd.read_csv('521-ventral-Total Ion Count.csv', delimiter=';', comment='#')
ventral_ko3y = pd.read_csv('522-ventral-Total Ion Count.csv', delimiter=';', comment='#')
# load corpus data
corpus_c1y = pd.read_csv('286-corpus-Total Ion Count.csv', delimiter=';', comment='#')
corpus_c2y = pd.read_csv('287-corpus-Total Ion Count.csv', delimiter=';', comment='#')
corpus_c3y = pd.read_csv('288-corpus-Total Ion Count.csv', delimiter=';', comment='#')
corpus_ko1y = pd.read_csv('509-corpus-Total Ion Count.csv', delimiter=';', comment='#')
corpus_ko2y = pd.read_csv('521-corpus-Total Ion Count.csv', delimiter=';', comment='#')
corpus_ko3y = pd.read_csv('522-corpus-Total Ion Count.csv', delimiter=';', comment='#')
# load caudate+cortex data
cc_c1y = pd.read_csv('286-cc-Total Ion Count.csv', delimiter=';', comment='#')
cc_c2y = pd.read_csv('287-cc-Total Ion Count.csv', delimiter=';', comment='#')
cc_c3y = pd.read_csv('288-cc-Total Ion Count.csv', delimiter=';', comment='#')
cc_ko1y = pd.read_csv('509-cc-Total Ion Count.csv', delimiter=';', comment='#')
cc_ko2y = pd.read_csv('521-cc-Total Ion Count.csv', delimiter=';', comment='#')
cc_ko3y = pd.read_csv('522-cc-Total Ion Count.csv', delimiter=';', comment='#')

# Reference
ref = pd.read_csv('lipids-correlation-fullset-youngmice.csv', delimiter=',')

# # Global variables
# Set variables
drift = 0.6
ion_count = 100
size_set = 8
size_sety = 6
grp_annot = ref.iloc[:ion_count,8]
name_annot = ref.iloc[:ion_count,7]
k=2

index = [0,2,3,6,13,14,15,38,48,51,56,10,16,17,20,29,37,41,60,36,69,80,87,88,93,99,31,32,34,39,46,54,59,66,67,72,76,77,78,79,81,82,83,84,85,89,90,92,95,97,35,44,45,52,58]

# mean
dist_ventricle = np.zeros([ion_count, size_set])
dist_ventricle_log2 = np.zeros([ion_count, size_set])
dist_cortex = np.zeros([ion_count, size_set])
dist_caudate = np.zeros([ion_count, size_set])
dist_cc = np.zeros([ion_count, size_set])
dist_corpus = np.zeros([ion_count, size_set])
dist_ventral = np.zeros([ion_count, size_set])

dist_ventricley = np.zeros([ion_count, size_sety])
dist_ventricley_log2 = np.zeros([ion_count, size_sety])
dist_cortexy = np.zeros([ion_count, size_sety])
dist_caudatey = np.zeros([ion_count, size_sety])
dist_ccy = np.zeros([ion_count, size_sety])
dist_corpusy = np.zeros([ion_count, size_sety])
dist_ventraly = np.zeros([ion_count, size_sety])

# resample mean
dist_ventricle_re = np.zeros([ion_count, size_set*k])
dist_ventricle_re_l2 = np.zeros([ion_count, size_set*k])
dist_cortex_re = np.zeros([ion_count, size_set*k])
dist_caudate_re = np.zeros([ion_count, size_set*k])
dist_cc_re = np.zeros([ion_count, size_set*k])
dist_corpus_re = np.zeros([ion_count, size_set*k])
dist_ventral_re = np.zeros([ion_count, size_set*k])

dist_ventricley_re = np.zeros([ion_count, size_sety*k])
dist_ventricley_re_l2 = np.zeros([ion_count, size_sety*k])
dist_cortexy_re = np.zeros([ion_count, size_sety*k])
dist_caudatey_re = np.zeros([ion_count, size_sety*k])
dist_ccy_re = np.zeros([ion_count, size_sety*k])
dist_corpusy_re = np.zeros([ion_count, size_sety*k])
dist_ventraly_re = np.zeros([ion_count, size_sety*k])

#std
std_ventricle = np.zeros([ion_count, size_set])
std_ventricle_log2 = np.zeros([ion_count, size_set])
std_cortex = np.zeros([ion_count, size_set])
std_caudate = np.zeros([ion_count, size_set])
std_cc = np.zeros([ion_count, size_set])
std_corpus = np.zeros([ion_count, size_set])
std_ventral = np.zeros([ion_count, size_set])

std_ventricley = np.zeros([ion_count, size_set])
std_ventricley_log2 = np.zeros([ion_count, size_set])
std_cortexy = np.zeros([ion_count, size_set])
std_caudatey = np.zeros([ion_count, size_set])
std_ccy = np.zeros([ion_count, size_set])
std_corpusy = np.zeros([ion_count, size_set])
std_ventraly = np.zeros([ion_count, size_set])

# Test statistic
stat_ventricle = np.zeros(100)
stat_cortex = np.zeros(100)
stat_caudate = np.zeros(100)
stat_cc = np.zeros(100)
stat_corpus = np.zeros(100)
stat_ventral = np.zeros(100)
stat_ventricley = np.zeros(100)
stat_cortexy = np.zeros(100)
stat_caudatey = np.zeros(100)
stat_ccy = np.zeros(100)
stat_corpusy = np.zeros(100)
stat_ventraly = np.zeros(100)

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
        col_std[i] = np.std(input.iloc[i,:])/col_mean[i]
    return col_mean, col_std

def ret_dist_log2(input, ref, index):
# Arguments:    
    #input - dataframe

# Returns:
    # np.append(col_mean, col_std)
        # col_mean - np.ndarray containing mean for each ion intensity distribution
        # col_std - np.ndarray containing std. dev. for each ion intensity distribution

    init = 0
    for k in range(index+1):
        if k < index:
            init+=ref[k]
        if k == index:
            fin = init + ref[k]

    col_mean = np.zeros(ion_count)
    col_std = np.zeros(ion_count)
    for i in range(ion_count):
        col_mean[i] = np.mean(input.iloc[i,init:fin])
        col_std[i] = np.std(input.iloc[i,init:fin])
    return col_mean, col_std

def ret_dist_resample(input, resample):
# Arguments:    
    #input - dataframe

# Returns:
    # np.append(col_mean, col_std)
        # col_means - np.ndarray containing mean for each ion intensity distribution

    col_means = np.zeros([ion_count, resample])
    for i in range(ion_count):
        for j in range(resample):
            col_means[i,j] = np.mean(random.choices(population=sorted(input.iloc[i,1:]), k=int(len(input.iloc[i,1:])/resample)))
    return col_means

def ret_dist_resample_l2(input, ref, index_r, resample):
# Arguments:    
    #input - dataframe
    # ref - reference indexer (identifies range of data)
    # index_r - pointer to data location in ref
    # resample - number of random means to augment
    
# Returns:
    # np.append(col_mean, col_std)
        # col_means - np.ndarray containing mean for each ion intensity distribution
    
    init = 0
    for q in range(index_r+1):
        if q < index_r:
            init+=ref[q]
        if q == index_r:
            fin = init + ref[q]

    col_means = np.zeros([ion_count, resample])
    for i in range(ion_count):
        for j in range(resample):
            #col_means[i,j] = np.mean(random.choices(population=sorted(input.iloc[i,1:]), k=int(len(input.iloc[i,1:])/resample)))
            col_means[i,j] = np.mean(random.sample(population=sorted(input.iloc[i,init:fin]), \
                                                   k=int(len(input.iloc[i,init:fin])/resample)))
    return col_means

def isolate_lipid_group(input, vector, size):
# Arguments:
    # input - dataframe
    # vector - group to sort

# Returns:
    # iso - np.ndarray, sorted means by lipid group
    # anno - np.ndarray, sorted annotations by lipid name within group

    j=[]
    for i in range(ion_count):
        if grp_annot[i] == vector:
            #print(i)
            j.append(i)

    iso = np.zeros([len(j), size])
    anno = []

    for i in range(len(j)):
        iso[i] = (input[j[i],:])
        anno.append(name_annot[j[i]])

    return iso, anno

def diff_compare_local(reference, input):
    output = np.zeros(np.shape(input)[0])
    for i in range(np.shape(input)[0]):
        output[i] = ((np.mean(input[i,:])-np.mean(reference[i,:])) / np.mean(reference[i,:]))
    return output

def diff_compare_local_l2(reference, input):
    output = np.zeros(np.shape(input)[0])
    for i in range(np.shape(input)[0]):
        output[i] = (np.mean(input[i,:])-np.mean(reference[i,:]))
    return output

def diff_compare_split(reference, input):
    output = np.zeros(np.shape(input)[0])
    for i in range(np.shape(input)[0]):
        a = ((np.mean(input[i,1:3])-np.mean(reference[i,0:3])) / np.mean(reference[i,0:3]))
        b = ((np.mean([input[i,0], input[i,3]*drift])-np.mean(reference[i,3])) / np.mean(reference[i,3]))
        output[i] = np.mean([a,b])
    return output

def diff_compare_split_l2(reference, input):
    output = np.zeros(np.shape(input)[0])
    for i in range(np.shape(input)[0]):
        a = (np.mean(input[i,1:3])-np.mean(reference[i,0:3]))
        b = (np.mean([input[i,0], input[i,3]])-np.mean(reference[i,3]))
        output[i] = np.mean([a,b])
    return output

def log2FC(reference, input):
    output = []
    for i in range(np.shape(input)[0]):
        output.append((np.mean(input[i,:])-np.mean(reference[i,:])))
    return output

def stat_compare_split(input):
    output = np.zeros(np.shape(input))
    for i in range(ion_count):
        output[i,14:16] = (input[i,14:16]*drift)-np.mean(input[i,6:8])
        output[i,8:10] = input[i, 8:10] - np.mean(input[i,6:8])
        output[i,0:6] = input[i,0:6] - np.mean(input[i,0:6])
        output[i,6:8] = input[i,6:8] - np.mean(input[i,6:8])
        output[i,10:14] = input[i,10:14] - np.mean(input[i,0:6])
    return output

def statistic(x,y,axis):
    return np.mean(x, axis=axis)-np.mean(y, axis=axis)


def stat_test(input, reference):
    s = stats.permutation_test((reference, input), statistic, batch=len(reference)-1, vectorized=True, permutation_type='independent', n_resamples=5e3, alternative='two-sided')
    return s.pvalue

# # ---- Main Program ------------------------------------------------------------------------------------------------- # #

index_ref =  [0, ventricle_c1.shape[1], ventricle_c2.shape[1], ventricle_c3.shape[1], ventricle_c4.shape[1], \
              ventricle_ko1.shape[1], ventricle_ko2.shape[1], ventricle_ko3.shape[1], ventricle_ko4.shape[1]]

index_refy =  [0, ventricle_c1y.shape[1], ventricle_c2y.shape[1], ventricle_c3y.shape[1], \
              ventricle_ko1y.shape[1], ventricle_ko2y.shape[1], ventricle_ko3y.shape[1]]

# Ventricle
# Data Org (Tuple)
ventricle_all = [ventricle_c1, ventricle_c2, ventricle_c3, ventricle_c4, ventricle_ko1, ventricle_ko2, ventricle_ko3, ventricle_ko4]
ventricle_ally = [ventricle_c1y, ventricle_c2y, ventricle_c3y, ventricle_ko1y, ventricle_ko2y, ventricle_ko3y]

# As DataFrame
vent_all = pd.concat([ventricle_c1, ventricle_c2.iloc[:,1:], ventricle_c3.iloc[:,1:], ventricle_c4.iloc[:,1:], \
                       ventricle_ko1.iloc[:,1:], ventricle_ko2.iloc[:,1:], ventricle_ko3.iloc[:,1:], \
                        ventricle_ko4.iloc[:,1:]], axis=1)
vent_ally = pd.concat([ventricle_c1y, ventricle_c2y.iloc[:,1:], ventricle_c3y.iloc[:,1:], \
                       ventricle_ko1y.iloc[:,1:], ventricle_ko2y.iloc[:,1:], ventricle_ko3y.iloc[:,1:]], axis=1)

vent_all = vent_all.replace(0, np.nan)
imputed_all = vent_all.fillna(vent_all.iloc[:,1:].median())
log2_all = np.log2(imputed_all.iloc[:,1:])

vent_ally = vent_ally.replace(0, np.nan)
imputed_ally = vent_ally.fillna(vent_ally.iloc[:,1:].median())
log2_ally = np.log2(imputed_ally.iloc[:,1:])

# Compute mean and resampling vectors
for l in range(size_set):
    dist_ventricle[:,l], std_ventricle[:,l] = ret_dist(ventricle_all[l].iloc[:,1:])
    dist_ventricle_log2[:,l], std_ventricle_log2[:,l]  = ret_dist_log2(log2_all, index_ref, l+1)
for l in range(size_sety):
    dist_ventricley[:,l], std_ventricley[:,l] = ret_dist(ventricle_ally[l].iloc[:,1:])
    dist_ventricley_log2[:,l], std_ventricley_log2[:,l] = ret_dist_log2(log2_ally, index_refy, l+1)
m=0
for l in range(size_set):
    dist_ventricle_re[:,m:m+k] = ret_dist_resample(ventricle_all[l], k)
    dist_ventricle_re_l2[:,m:m+k] = ret_dist_resample_l2(log2_all, index_ref, l+1, k)
    m+=2
m=0
for l in range(size_sety):
    dist_ventricley_re[:,m:m+k] = ret_dist_resample(ventricle_ally[l], k)
    dist_ventricley_re_l2[:,m:m+k] = ret_dist_resample_l2(log2_ally, index_refy, l+1, k)

    m+=2

dist_ventricle = dist_ventricle_log2
dist_ventricley = dist_ventricley_log2
dist_ventricle_re = dist_ventricle_re_l2
dist_ventricley_re = dist_ventricley_re_l2

# Global mean differential computation (Barplots)
ventricle_FA, annot_FA = isolate_lipid_group(dist_ventricle, 'FA', size_set)
ventricle_ST, annot_ST = isolate_lipid_group(dist_ventricle, 'ST', size_set)
ventricle_SP, annot_SP = isolate_lipid_group(dist_ventricle, 'SP', size_set)
ventricle_GP, annot_GP = isolate_lipid_group(dist_ventricle, 'GP', size_set)
ventricle_GL, annot_GL = isolate_lipid_group(dist_ventricle, 'GL', size_set)

ventricle_FAy, annot_FA = isolate_lipid_group(dist_ventricley, 'FA', size_sety)
ventricle_STy, annot_ST = isolate_lipid_group(dist_ventricley, 'ST', size_sety)
ventricle_SPy, annot_SP = isolate_lipid_group(dist_ventricley, 'SP', size_sety)
ventricle_GPy, annot_GP = isolate_lipid_group(dist_ventricley, 'GP', size_sety)
ventricle_GLy, annot_GL = isolate_lipid_group(dist_ventricley, 'GL', size_sety)

# fa_ventricle_diff = diff_compare_split(ventricle_FA[:,0:4], ventricle_FA[:,4:8])
# st_ventricle_diff = diff_compare_split(ventricle_ST[:,0:4], ventricle_ST[:,4:8])
# sp_ventricle_diff = diff_compare_split(ventricle_SP[:,0:4], ventricle_SP[:,4:8])
# gp_ventricle_diff = diff_compare_split(ventricle_GP[:,0:4], ventricle_GP[:,4:8])
# gl_ventricle_diff = diff_compare_split(ventricle_GL[:,0:4], ventricle_GL[:,4:8])

fa_ventricle_diff = diff_compare_split_l2(ventricle_FA[:,0:4], ventricle_FA[:,4:8])
st_ventricle_diff = diff_compare_split_l2(ventricle_ST[:,0:4], ventricle_ST[:,4:8])
sp_ventricle_diff = diff_compare_split_l2(ventricle_SP[:,0:4], ventricle_SP[:,4:8])
gp_ventricle_diff = diff_compare_split_l2(ventricle_GP[:,0:4], ventricle_GP[:,4:8])
gl_ventricle_diff = diff_compare_split_l2(ventricle_GL[:,0:4], ventricle_GL[:,4:8])

# fa_ventricle_diffy = diff_compare_local(ventricle_FAy[:,0:3], ventricle_FAy[:,3:6])
# st_ventricle_diffy = diff_compare_local(ventricle_STy[:,0:3], ventricle_STy[:,3:6])
# sp_ventricle_diffy = diff_compare_local(ventricle_SPy[:,0:3], ventricle_SPy[:,3:6])
# gp_ventricle_diffy = diff_compare_local(ventricle_GPy[:,0:3], ventricle_GPy[:,3:6])
# gl_ventricle_diffy = diff_compare_local(ventricle_GLy[:,0:3], ventricle_GLy[:,3:6])

fa_ventricle_diffy = diff_compare_local_l2(ventricle_FAy[:,0:3], ventricle_FAy[:,3:6])
st_ventricle_diffy = diff_compare_local_l2(ventricle_STy[:,0:3], ventricle_STy[:,3:6])
sp_ventricle_diffy = diff_compare_local_l2(ventricle_SPy[:,0:3], ventricle_SPy[:,3:6])
gp_ventricle_diffy = diff_compare_local_l2(ventricle_GPy[:,0:3], ventricle_GPy[:,3:6])
gl_ventricle_diffy = diff_compare_local_l2(ventricle_GLy[:,0:3], ventricle_GLy[:,3:6])

#ventricle_FA_std, annot_FA = isolate_lipid_group(std_ventricle, 'FA', size_set)
#ventricle_FA_stdy, annot_FA = isolate_lipid_group(std_ventricley, 'FA', size_set)
#fa_ventricle_diff_std = diff_compare_split(ventricle_FA_std[:,0:4], ventricle_FA_std[:,4:8])
#fa_ventricle_diff_stdy = diff_compare_local(ventricle_FA_stdy[:,0:3], ventricle_FA_stdy[:,3:6])

# Test statistic computation
data_ventricle_stat = stat_compare_split(dist_ventricle_re)

for i in range(ion_count):
    stat_ventricley[i] = stat_test(dist_ventricley_re[i,0:6], dist_ventricley_re[i,6:12])
    stat_ventricle[i] = stat_test(data_ventricle_stat[i, 0:8], data_ventricle_stat[i,8:16])

diff_ventricle = diff_compare_split_l2(dist_ventricle[:,0:4], dist_ventricle[:,4:8])
diff_ventricley = diff_compare_local_l2(dist_ventricley[:,0:3], dist_ventricley[:,3:6])


index = [0,2,3,6,13,14,15,38,48,51,56,10,16,17,20,29,37,41,60,36,69,80,87,88,93,99,31,32,34,39,46,54,59,66,67,72,76,77,78,79,81,82,83,84,85,89,90,92,95,97,35,44,45,52,58]
FA_index = [0,2,3,6,13,14,15,38,48,51,56]
ST_index = [10,16,17,20,29,37,41,60]
SP_index = [36,69,80,87,88,93,99]
GP_index = [31,32,34,39,46,54,59,66,67,72,76,77,78,79,81,82,83,84,85,89,90,92,95,97]
GL_index = [35,44,45,52,58]

alpha = 0.05

# # null outliers
# st_ventricle_diffy[2] = 0.1
# st_ventricle_diffy[3] = 0.1

# weight_FA = np.zeros(len(annot_FA))
# weight_FAy = np.zeros(len(annot_FA))
# for i in range(len(annot_FA)):
#     weight_FA[i] = fa_ventricle_diff[i]*(np.log(stat_ventricle[FA_index[i]])/np.log(alpha))
#     weight_FAy[i] = fa_ventricle_diffy[i]*(np.log(stat_ventricley[FA_index[i]])/np.log(alpha))

# weight_ST = np.zeros(len(annot_ST))
# weight_STy = np.zeros(len(annot_ST))
# for i in range(len(annot_ST)):
#     weight_ST[i] = st_ventricle_diff[i]*(np.log(stat_ventricle[ST_index[i]])/np.log(alpha))
#     weight_STy[i] = st_ventricle_diffy[i]*(np.log(stat_ventricley[ST_index[i]])/np.log(alpha))

# weight_SP = np.zeros(len(annot_SP))
# weight_SPy = np.zeros(len(annot_SP))
# for i in range(len(annot_SP)):
#     weight_SP[i] = sp_ventricle_diff[i]*(np.log(stat_ventricle[SP_index[i]])/np.log(alpha))
#     weight_SPy[i] = sp_ventricle_diffy[i]*(np.log(stat_ventricley[SP_index[i]])/np.log(alpha))

# weight_GL = np.zeros(len(annot_GL))
# weight_GLy = np.zeros(len(annot_GL))
# for i in range(len(annot_GL)):
#     weight_GL[i] = gl_ventricle_diff[i]*(np.log(stat_ventricle[GL_index[i]])/np.log(alpha))
#     weight_GLy[i] = gl_ventricle_diffy[i]*(np.log(stat_ventricley[GL_index[i]])/np.log(alpha))

# weight_GP = np.zeros(len(annot_GP))
# weight_GPy = np.zeros(len(annot_GP))
# for i in range(len(annot_GP)):
#     weight_GP[i] = gp_ventricle_diff[i]*(np.log(stat_ventricle[GP_index[i]])/np.log(alpha))
#     weight_GPy[i] = gp_ventricle_diffy[i]*(np.log(stat_ventricley[GP_index[i]])/np.log(alpha))

# # Statistically weighted plots
# fig, ax = plt.subplots(nrows=1, ncols=5)
# fig.suptitle('12mo. Glut1KO')
# ax[0].bar(x=annot_FA, height=weight_FA, label='fatty acyls', color='darkred')
# ax[0].set(ylabel='Statistically Weighted Differential Signal')
# ax[0].set(xlabel='fatty acyls')
# ax[1].bar(x=annot_ST, height=weight_ST, label='sterols', color='darkmagenta')
# ax[1].set(xlabel='sterols')
# ax[2].bar(x=annot_SP, height=weight_SP, label='sphingolipids', color='slategrey')
# ax[2].set(xlabel='sphingolipids')
# ax[3].bar(x=annot_GP, height=weight_GP, label='glycerophospholipids', color='black')
# ax[3].set(xlabel='glycerophospholipids')
# ax[4].bar(x=annot_GL, height=weight_GL, label='glycerolipids', color='darkgrey')
# ax[4].set(xlabel='glycerolipids')
# for k in range(len(ax)):
#     ax[k].tick_params(axis='x', labelsize=2, labelrotation=80)
#     ax[k].set(ylim=[-1,1])
#     if k > 0 :
#         ax[k].set(yticks=[])
#     # ax[k].axhline(0.075,color="grey",linestyle="--")
#     # ax[k].axhline(-0.075,color="grey",linestyle="--")
# fig.align_labels()

# # Statistically weighted plots
# fig, ax = plt.subplots(nrows=1, ncols=5)
# fig.suptitle('1mo. Glut1KO')
# ax[0].bar(x=annot_FA, height=weight_FAy, label='fatty acyls', color='darkred')
# ax[0].set(ylabel='Statistically Weighted Differential Signal')
# ax[0].set(xlabel='fatty acyls')
# ax[1].bar(x=annot_ST, height=weight_STy, label='sterols', color='darkmagenta')
# ax[1].set(xlabel='sterols')
# ax[2].bar(x=annot_SP, height=weight_SPy, label='sphingolipids', color='slategrey')
# ax[2].set(xlabel='sphingolipids')
# ax[3].bar(x=annot_GP, height=weight_GPy, label='glycerophospholipids', color='black')
# ax[3].set(xlabel='glycerophospholipids')
# ax[4].bar(x=annot_GL, height=weight_GLy, label='glycerolipids', color='darkgrey')
# ax[4].set(xlabel='glycerolipids')
# for k in range(len(ax)):
#     ax[k].tick_params(axis='x', labelsize=2, labelrotation=80)
#     ax[k].set(ylim=[-1,1])
#     if k > 0 :
#         ax[k].set(yticks=[])
#     # ax[k].axhline(0.075,color="grey",linestyle="--")
#     # ax[k].axhline(-0.075,color="grey",linestyle="--")

# fig.align_labels()

# fig, ax = plt.subplots(nrows=1, ncols=5)
# fig.suptitle('Difference (12mo. from 1mo.)')
# ax[0].bar(x=annot_FA, height=weight_FA-weight_FAy, label='fatty acyls', color='darkred')
# ax[0].set(ylabel='Age-Relative Weighted Differential Signal')
# ax[0].set(xlabel='fatty acyls')
# ax[1].bar(x=annot_ST, height=weight_ST-weight_STy, label='sterols', color='darkmagenta')
# ax[1].set(xlabel='sterols')
# ax[2].bar(x=annot_SP, height=weight_SP-weight_SPy, label='sphingolipids', color='slategrey')
# ax[2].set(xlabel='sphingolipids')
# ax[3].bar(x=annot_GP, height=weight_GP-weight_GPy, label='glycerophospholipids', color='black')
# ax[3].set(xlabel='glycerophospholipids')
# ax[4].bar(x=annot_GL, height=weight_GL-weight_GLy, label='glycerolipids', color='darkgrey')
# ax[4].set(xlabel='glycerolipids')
# for k in range(len(ax)):
#     ax[k].tick_params(axis='x', labelsize=2, labelrotation=80)
#     ax[k].set(ylim=[-1,1])
#     if k > 0 :
#         ax[k].set(yticks=[])
# fig.align_labels()

# fig, ax = plt.subplots()
# plt.plot(dist_ventricle[:,3], dist_ventricle[:,7]/dist_ventricle[:,3])

# fig, ax = plt.subplots(nrows=1, ncols=5)
# fig.suptitle('12mo. Glut1KO')
# ax[0].bar(x=annot_FA, height=fa_ventricle_diff, label='fatty acyls', color='darkred')
# ax[0].set(ylabel='Mean Relative Differential Signal')
# ax[0].set(xlabel='fatty acyls')
# ax[1].bar(x=annot_ST, height=st_ventricle_diff, label='sterols', color='darkmagenta')
# ax[1].set(xlabel='sterols')
# ax[2].bar(x=annot_SP, height=sp_ventricle_diff, label='sphingolipids', color='slategrey')
# ax[2].set(xlabel='sphingolipids')
# ax[3].bar(x=annot_GP, height=gp_ventricle_diff, label='glycerophospholipids', color='black')
# ax[3].set(xlabel='glycerophospholipids')
# ax[4].bar(x=annot_GL, height=gl_ventricle_diff, label='glycerolipids', color='darkgrey')
# ax[4].set(xlabel='glycerolipids')

# for k in range(len(ax)):
#     ax[k].tick_params(axis='x', labelsize=2, labelrotation=80)
#     ax[k].set(ylim=[-1,1])
#     if k > 0 :
#         ax[k].set(yticks=[])
# fig.align_labels()

# fig, ax = plt.subplots(nrows=1, ncols=5)
# fig.suptitle('1mo. Glut1KO')
# ax[0].bar(x=annot_FA, height=fa_ventricle_diffy, label='fatty acyls', color='darkred')
# ax[0].set(ylabel='Mean Relative Differential Signal')
# ax[0].set(xlabel='fatty acyls')
# ax[1].bar(x=annot_ST, height=st_ventricle_diffy, label='sterols', color='darkmagenta')
# ax[1].set(xlabel='sterols')
# ax[2].bar(x=annot_SP, height=sp_ventricle_diffy, label='sphingolipids', color='slategrey')
# ax[2].set(xlabel='sphingolipids')
# ax[3].bar(x=annot_GP, height=gp_ventricle_diffy, label='glycerophospholipids', color='black')
# ax[3].set(xlabel='glycerophospholipids')
# ax[4].bar(x=annot_GL, height=gl_ventricle_diffy, label='glycerolipids', color='darkgrey')
# ax[4].set(xlabel='glycerolipids')
# for k in range(len(ax)):
#     ax[k].tick_params(axis='x', labelsize=2, labelrotation=80)
#     ax[k].set(ylim=[-1,1])
#     if k > 0 :
#         ax[k].set(yticks=[])
# fig.align_labels()

# fig, ax = plt.subplots(nrows=1, ncols=5)
# fig.suptitle('Difference (12mo. from 1mo.)')
# ax[0].bar(x=annot_FA, height=fa_ventricle_diff-fa_ventricle_diffy, label='fatty acyls', color='darkred')
# ax[0].set(ylabel='Mean Relative Differential Signal')
# ax[0].set(xlabel='fatty acyls')
# ax[1].bar(x=annot_ST, height=st_ventricle_diff-st_ventricle_diffy, label='sterols', color='darkmagenta')
# ax[1].set(xlabel='sterols')
# ax[2].bar(x=annot_SP, height=sp_ventricle_diff-sp_ventricle_diffy, label='sphingolipids', color='slategrey')
# ax[2].set(xlabel='sphingolipids')
# ax[3].bar(x=annot_GP, height=gp_ventricle_diff-gp_ventricle_diffy, label='glycerophospholipids', color='black')
# ax[3].set(xlabel='glycerophospholipids')
# ax[4].bar(x=annot_GL, height=gl_ventricle_diff-gl_ventricle_diffy, label='glycerolipids', color='darkgrey')
# ax[4].set(xlabel='glycerolipids')
# for k in range(len(ax)):
#     ax[k].tick_params(axis='x', labelsize=2, labelrotation=80)
#     ax[k].set(ylim=[-1,1])
#     if k > 0 :
#         ax[k].set(yticks=[])
# fig.align_labels()

colors = []
for i in range(ion_count):
    if grp_annot[i] == 'FA':
        colors.append('darkred')
    elif grp_annot[i] == 'ST':
        colors.append('darkmagenta')
    elif grp_annot[i] == 'SP':
        colors.append('slategrey')
    elif grp_annot[i] == 'GP':
        colors.append('black')
    elif grp_annot[i] == 'GL':
        colors.append('darkgrey')
    else:
        colors.append('whitesmoke')

# Volcano Plots
# log_diff_ventricle = np.log2(1+diff_ventricle)
# log_diff_ventricley = np.log2(1+diff_ventricley)

fig, ax = plt.subplots()
ax.scatter(x=diff_ventricle, y=-np.log10(stat_ventricle), s=16, color=colors, label=None)
ax.scatter(x=0, y=0, color='darkred', s=16, label='Fatty Acyls')
ax.scatter(x=0, y=0, color='darkmagenta', s=16, label='Sterols')
ax.scatter(x=0, y=0, color='slategray', s=16, label='Sphingolipids')
ax.scatter(x=0, y=0, color='black', s=16, label='Glycerophospholipids')
ax.scatter(x=0, y=0, color='darkgrey', s=16, label='Glycerolipids')
plt.axvline(-0.1,color="grey",linestyle="--")
plt.axvline(0.1,color="grey",linestyle="--")
plt.axhline(1.3,color="grey",linestyle="--")
ax.legend(fontsize=13)
ax.set(xlim=[-2.75,1],ylim=[-0.1, 4], ylabel='-log10(pvalue)', xlabel='Estimation of Fold-Change Lipid Expression', title='Observed Lipid Dysregulation at 12mo.')
for i in index:
    if -np.log10(stat_ventricle[i]) > 0.75:
        plt.text(x=1-(1/np.power(10,diff_ventricle[i])), y=-np.log10(stat_ventricle[i]), s=name_annot[i], fontsize=9)


fig, ax = plt.subplots()
ax.scatter(x=diff_ventricley, y=-np.log10(stat_ventricley), s=16, color=colors, label=None)
ax.scatter(x=0, y=0, color='darkred', s=16, label='Fatty Acyls')
ax.scatter(x=0, y=0, color='darkmagenta', s=16, label='Sterols')
ax.scatter(x=0, y=0, color='slategray', s=16, label='Sphingolipids')
ax.scatter(x=0, y=0, color='black', s=16, label='Glycerophospholipids')
ax.scatter(x=0, y=0, color='darkgrey', s=16, label='Glycerolipids')
plt.axvline(-0.1,color="grey",linestyle="--")
plt.axvline(0.1,color="grey",linestyle="--")
plt.axhline(1.3,color="grey",linestyle="--")
ax.legend(fontsize=13)
ax.set(xlim=[-2.75,1],ylim=[-0.1, 4], ylabel='-log10(pvalue)', xlabel='Estimation of Fold-Change Lipid Expression', title='Observed Lipid Dysregulation at 1mo.')
for i in index:
    if -np.log10(stat_ventricley[i]) > 0.75:
        plt.text(x=1-(1/np.power(10,diff_ventricley[i])), y=-np.log10(stat_ventricley[i]), s=name_annot[i], fontsize=9)


# fig, ax = plt.subplots()
# ax.scatter(x=1-(1/np.power(10,diff_ventricle)), y=-np.log10(stat_ventricle), s=16, color='darkred', label='12mo.')
# ax.scatter(x=1-(1/np.power(10,diff_ventricley)), y=-np.log10(stat_ventricley), s=16, color='lightsteelblue', label='1mo')
# plt.axvline(-0.2,color="grey",linestyle="--")
# plt.axvline(0.2,color="grey",linestyle="--")
# plt.axhline(1.3,color="grey",linestyle="--")
# ax.set(xlim=[-2.75,1],ylim=[-0.1, 4], ylabel='-log10(pvalue)', xlabel='Estimation of Fold-Change Lipid Expression', title='Overlay')
# ax.legend()

plt.show()


# for i in range(len(fa_ventricle_diff)):
#     print(fa_ventricle_diff[i])
# print(' ')
# for i in range(len(fa_ventricle_diffy)):
#     print(fa_ventricle_diffy[i])
# print(' ')
# print(' ')
# for i in range(len(st_ventricle_diff)):
#     print(st_ventricle_diff[i])
# print(' ')
# for i in range(len(st_ventricle_diffy)):
#     print(st_ventricle_diffy[i])
# print(' ')
# print(' ')
# for i in range(len(sp_ventricle_diff)):
#     print(sp_ventricle_diff[i])
# print(' ')
# for i in range(len(sp_ventricle_diffy)):
#     print(sp_ventricle_diffy[i])
# print(' ')
# print(' ')
# for i in range(len(gp_ventricle_diff)):
#     print(gp_ventricle_diff[i])
# print(' ')
# for i in range(len(gp_ventricle_diffy)):
#     print(gp_ventricle_diffy[i])
# print(' ')
# print(' ')
# for i in range(len(gl_ventricle_diff)):
#     print(gl_ventricle_diff[i])
# print(' ')
# for i in range(len(gl_ventricle_diffy)):
#     print(gl_ventricle_diffy[i])
# print(' ')
# print(' ')

# for i in index:
#     print(1-(1/np.power(10,diff_ventricle[i])))
# print(' ')
# for i in index:
#     print(1-(1/np.power(10,diff_ventricley[i])))
# print(' ')
# print(' ')

# for i in index:
#     print(stat_ventricle[i])
# print(' ')
# for i in index:
#     print(stat_ventricley[i])
# print(' ')
# print(' ')

