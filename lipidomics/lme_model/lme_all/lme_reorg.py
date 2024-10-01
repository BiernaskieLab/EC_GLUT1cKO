# # Libraries
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import numpy as np
import pandas as pd

import statsmodels.api as sm
import statsmodels.formula.api as sm

from scipy.stats import mannwhitneyu

from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler

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
size_set_old = 8
size_set_young = 6
grp_annot = ref.iloc[:ion_count,8]
name_annot = ref.iloc[:ion_count,7]
k=2

index = [0,2,3,6,13,14,15,38,48,51,56,10,16,17,20,29,37,41,60,36,69,80,87,88,93,99,31,32,34,39,46,54,59,66,67,72,76,77,78,79,81,82,83,84,85,89,90,92,95,97,35,44,45,52,58]

scalerR = RobustScaler()
scalerM = MinMaxScaler()


def scale_data(input):
    k=0
    output = input
    for i in range(input.shape[0]):
        output.iloc[i] = (input.iloc[i] - np.min(input.iloc[i])) / ( np.max(input.iloc[i])-np.min(input.iloc[i]))
    return pd.DataFrame(output)

def differential_compare_means_by_condition(input_means_by_condition):
     output = np.zeros(input_means_by_condition.shape[0])
     for i in range(input_means_by_condition.shape[0]):
        output[i] = (input_means_by_condition.iloc[i,1] - input_means_by_condition.iloc[i,0]) / input_means_by_condition.iloc[i,0]
     return pd.DataFrame(output)
     

def compute_means(input, indices, size):
    output_mouse = np.zeros([input.shape[0], size])
    output_condition = np.zeros([input.shape[0], 2])

    for i in range(input.shape[0]):
        ind = 0
        for j in range(size):
            output_mouse[i,j] = np.nanmean(input.iloc[i,ind:ind+indices[j+1]])
            ind+=indices[j+1]
    for i in range(input.shape[0]):            
        output_condition[i,0] = np.mean(output_mouse[i,0:int(size/2)])
        output_condition[i,1] = np.mean(output_mouse[i,int(size/2):size])
    return pd.DataFrame(output_mouse), pd.DataFrame(output_condition)

def nonparametric_test(input, size):
    output = np.zeros(input.shape[0])
    for i in range(input.shape[0]):
        stat, pval = mannwhitneyu(input.iloc[i,0:int(size/2)], input.iloc[i,int(size/2):size], method='auto')
        output[i] = pval
    return pd.DataFrame(output)
    
          
# #
# # ---- Main Program ------------------------------------------------------------------------------------------------- # #
# #

# # Ventricle
# vent_index_ref_old =  [0, ventricle_c1.shape[1]-1, ventricle_c2.shape[1]-1, ventricle_c3.shape[1]-1, ventricle_c4.shape[1]-1, \
#               ventricle_ko1.shape[1]-1, ventricle_ko2.shape[1]-1, ventricle_ko3.shape[1]-1, ventricle_ko4.shape[1]-1]

# vent_index_ref_young =  [0, ventricle_c1y.shape[1]-1, ventricle_c2y.shape[1]-1, ventricle_c3y.shape[1]-1, \
#               ventricle_ko1y.shape[1]-1, ventricle_ko2y.shape[1]-1, ventricle_ko3y.shape[1]-1]

# # Data Org (Tuple)
# vent_all_tuple = [ventricle_c1, ventricle_c2, ventricle_c3, ventricle_c4, ventricle_ko1, ventricle_ko2, ventricle_ko3, ventricle_ko4]
# vent_all_young_tuple = [ventricle_c1y, ventricle_c2y, ventricle_c3y, ventricle_ko1y, ventricle_ko2y, ventricle_ko3y]

# # As DataFrames
# data_index = 1

# vent_all_old = pd.concat([ventricle_c1.iloc[:,data_index:], ventricle_c2.iloc[:,data_index:], ventricle_c3.iloc[:,data_index:], \
#                         ventricle_c4.iloc[:,data_index:], ventricle_ko1.iloc[:,data_index:], ventricle_ko2.iloc[:,data_index:], \
#                         ventricle_ko3.iloc[:,data_index:], ventricle_ko4.iloc[:,data_index:]], axis=1)

# vent_all_young = pd.concat([ventricle_c1y.iloc[:,data_index:], ventricle_c2y.iloc[:,data_index:], ventricle_c3y.iloc[:,data_index:], \
#                         ventricle_ko1y.iloc[:,data_index:], ventricle_ko2y.iloc[:,data_index:], \
#                         ventricle_ko3y.iloc[:,data_index:]], axis=1)

# vent_all_old = vent_all_old.replace(0, np.nan)
# vent_imputed_all_old = vent_all_old.fillna(vent_all_old.median())
# vent_log_all_old = np.log2(vent_imputed_all_old)
# vent_all_scale_old = scale_data(vent_imputed_all_old)

# data_proc_old = vent_log_all_old

# output_old = pd.DataFrame(np.copy(data_proc_old))
# old_cutoff = np.sum(vent_index_ref_old[0:5])

# ko_wt_old = []
# for i in range(np.sum(vent_index_ref_old)+1):
#     if i < old_cutoff:
#         ko_wt_old.append(str(0)) #str('control'))
#     if i > old_cutoff:
#         ko_wt_old.append(str(1)) #str('knockout'))

# mouse_id_old = []
# track = 0
# for i in range(size_set_old):
#     for j in range(vent_index_ref_old[i+1]):
#         mouse_id_old.append(str(i+1))
#         track+=1

# data_annot = name_annot.to_list()
# data_annot.insert(0, 'condition')
# data_annot.insert(0, 'mouse')
# output_old.loc[-1] = ko_wt_old
# output_old.loc[-2] = mouse_id_old
# output_old.index = output_old.index + 2
# output_old = output_old.sort_index()
# output_old.index = data_annot

# data_out_old = output_old.transpose()
# data_out_old.to_csv('ventricle_data_old.csv')


# vent_all_young = vent_all_young.replace(0, np.nan)
# vent_imputed_all_young = vent_all_young.fillna(vent_all_young.median())
# vent_log_all_young = np.log2(vent_imputed_all_young)
# vent_all_scale_young = scale_data(vent_imputed_all_young)

# data_proc_young = vent_log_all_young

# output_young = pd.DataFrame(np.copy(data_proc_young))
# young_cutoff = np.sum(vent_index_ref_young[0:4])

# ko_wt_young = []
# for i in range(np.sum(vent_index_ref_young)+1):
#     if i < young_cutoff:
#         ko_wt_young.append(str(0)) #str('control'))
#     if i > young_cutoff:
#         ko_wt_young.append(str(1)) #str('knockout'))

# mouse_id_young = []
# track = 0
# for i in range(size_set_young):
#     for j in range(vent_index_ref_young[i+1]):
#         mouse_id_young.append(str(i+1))
#         track+=1

# data_annot = name_annot.to_list()
# data_annot.insert(0, 'condition')
# data_annot.insert(0, 'mouse')
# output_young.loc[-1] = ko_wt_young
# output_young.loc[-2] = mouse_id_young
# output_young.index = output_young.index + 2
# output_young = output_young.sort_index()
# output_young.index = data_annot

# data_out_young = output_young.transpose()
# data_out_young.to_csv('ventricle_data_young.csv')


# # Cortex
# cortex_index_ref_old =  [0, cortex_c1.shape[1], cortex_c2.shape[1]-1, cortex_c3.shape[1]-1, cortex_c4.shape[1]-1, \
#               cortex_ko1.shape[1]-1, cortex_ko2.shape[1]-1, cortex_ko3.shape[1]-1, cortex_ko4.shape[1]-1]

# cortex_index_ref_young =  [0, cortex_c1y.shape[1], cortex_c2y.shape[1]-1, cortex_c3y.shape[1]-1, \
#               cortex_ko1y.shape[1]-1, cortex_ko2y.shape[1]-1, cortex_ko3y.shape[1]-1]

# # # Data Org (Tuple)
# # cortex_all_tuple_old = [cortex_c1, cortex_c2, cortex_c3, cortex_c4, cortex_ko1, cortex_ko2, cortex_ko3, cortex_ko4]
# # cortex_all_young_tuple = [cortex_c1y, cortex_c2y, cortex_c3y, cortex_ko1y, cortex_ko2y, cortex_ko3y]

# # As DataFrame
# cortex_all_old = pd.concat([cortex_c1, cortex_c2.iloc[:,1:], cortex_c3.iloc[:,1:], cortex_c4.iloc[:,1:], \
#                        cortex_ko1.iloc[:,1:], cortex_ko2.iloc[:,1:], cortex_ko3.iloc[:,1:], cortex_ko4.iloc[:,1:]], axis=1)
# cortex_all_young = pd.concat([cortex_c1y, cortex_c2y.iloc[:,1:], cortex_c3y.iloc[:,1:], \
#                        cortex_ko1y.iloc[:,1:], cortex_ko2y.iloc[:,1:], cortex_ko3y.iloc[:,1:]], axis=1)

# cortex_all_old = cortex_all_old.replace(0, np.nan)
# cortex_imputed_all_old = cortex_all_old.fillna(cortex_all_old.median())
# cortex_log_all_old = np.log2(cortex_imputed_all_old)
# cortex_all_scale_old = scale_data(cortex_imputed_all_old)

# data_proc_old = cortex_log_all_old

# output_old = pd.DataFrame(np.copy(data_proc_old))
# old_cutoff = np.sum(cortex_index_ref_old[0:5])

# ko_wt_old = []
# for i in range(np.sum(cortex_index_ref_old)):
#     if i < old_cutoff:
#         ko_wt_old.append(str(0)) #str('control'))
#     if i >= old_cutoff:
#         ko_wt_old.append(str(1)) #str('knockout'))

# mouse_id_old = []
# track = 0
# for i in range(size_set_old):
#     for j in range(cortex_index_ref_old[i+1]):
#         mouse_id_old.append(str(i+1))
#         track+=1

# data_annot = name_annot.to_list()
# data_annot.insert(0, 'condition')
# data_annot.insert(0, 'mouse')
# output_old.loc[-1] = ko_wt_old
# output_old.loc[-2] = mouse_id_old
# output_old.index = output_old.index + 2
# output_old = output_old.sort_index()
# output_old.index = data_annot

# data_out_old = output_old.transpose()
# data_out_old.to_csv('cortex_data_old.csv')


# cortex_all_young = cortex_all_young.replace(0, np.nan)
# cortex_imputed_all_young = cortex_all_young.fillna(cortex_all_young.median())
# cortex_log_all_young = np.log2(cortex_imputed_all_young)
# cortex_all_scale_young = scale_data(cortex_imputed_all_young)

# data_proc_young = cortex_log_all_young

# output_young = pd.DataFrame(np.copy(data_proc_young))
# young_cutoff = np.sum(cortex_index_ref_young[0:4])

# ko_wt_young = []
# for i in range(np.sum(cortex_index_ref_young)+1):
#     if i < young_cutoff:
#         ko_wt_young.append(str(0)) #str('control'))
#     if i > young_cutoff:
#         ko_wt_young.append(str(1)) #str('knockout'))

# print(len(ko_wt_young))
# mouse_id_young = []
# track = 0
# for i in range(size_set_young):
#     for j in range(cortex_index_ref_young[i+1]):
#         mouse_id_young.append(str(i+1))
#         track+=1
# print(len(mouse_id_young))
# print(output_young.shape)
# data_annot = name_annot.to_list()
# data_annot.insert(0, 'condition')
# data_annot.insert(0, 'mouse')
# output_young.loc[-1] = ko_wt_young
# output_young.loc[-2] = mouse_id_young
# output_young.index = output_young.index + 2
# output_young = output_young.sort_index()
# output_young.index = data_annot

# data_out_young = output_young.transpose()
# data_out_young.to_csv('cortex_data_young.csv')


# Cortex
corpus_index_ref_old =  [0, corpus_c1.shape[1], corpus_c2.shape[1]-1, corpus_c3.shape[1]-1, corpus_c4.shape[1]-1, \
              corpus_ko1.shape[1]-1, corpus_ko2.shape[1]-1, corpus_ko3.shape[1]-1, corpus_ko4.shape[1]-1]

corpus_index_ref_young =  [0, corpus_c1y.shape[1], corpus_c2y.shape[1]-1, corpus_c3y.shape[1]-1, \
              corpus_ko1y.shape[1]-1, corpus_ko2y.shape[1]-1, corpus_ko3y.shape[1]-1]

# As DataFrame
corpus_all_old = pd.concat([corpus_c1, corpus_c2.iloc[:,1:], corpus_c3.iloc[:,1:], corpus_c4.iloc[:,1:], \
                       corpus_ko1.iloc[:,1:], corpus_ko2.iloc[:,1:], corpus_ko3.iloc[:,1:], corpus_ko4.iloc[:,1:]], axis=1)
corpus_all_young = pd.concat([corpus_c1y, corpus_c2y.iloc[:,1:], corpus_c3y.iloc[:,1:], \
                       corpus_ko1y.iloc[:,1:], corpus_ko2y.iloc[:,1:], corpus_ko3y.iloc[:,1:]], axis=1)

corpus_all_old = corpus_all_old.replace(0, np.nan)
corpus_imputed_all_old = corpus_all_old.fillna(corpus_all_old.median())
corpus_log_all_old = np.log2(corpus_imputed_all_old)
corpus_all_scale_old = scale_data(corpus_imputed_all_old)

data_proc_old = corpus_log_all_old

output_old = pd.DataFrame(np.copy(data_proc_old))
old_cutoff = np.sum(corpus_index_ref_old[0:5])

ko_wt_old = []
for i in range(np.sum(corpus_index_ref_old)):
    if i < old_cutoff:
        ko_wt_old.append(str(0)) #str('control'))
    if i >= old_cutoff:
        ko_wt_old.append(str(1)) #str('knockout'))

mouse_id_old = []
track = 0
for i in range(size_set_old):
    for j in range(corpus_index_ref_old[i+1]):
        mouse_id_old.append(str(i+1))
        track+=1

data_annot = name_annot.to_list()
data_annot.insert(0, 'condition')
data_annot.insert(0, 'mouse')
output_old.loc[-1] = ko_wt_old
output_old.loc[-2] = mouse_id_old
output_old.index = output_old.index + 2
output_old = output_old.sort_index()
output_old.index = data_annot

data_out_old = output_old.transpose()
data_out_old.to_csv('corpus_data_old.csv')


corpus_all_young = corpus_all_young.replace(0, np.nan)
corpus_imputed_all_young = corpus_all_young.fillna(corpus_all_young.median())
corpus_log_all_young = np.log2(corpus_imputed_all_young)
corpus_all_scale_young = scale_data(corpus_imputed_all_young)

data_proc_young = corpus_log_all_young

output_young = pd.DataFrame(np.copy(data_proc_young))
young_cutoff = np.sum(corpus_index_ref_young[0:4])

ko_wt_young = []
for i in range(np.sum(corpus_index_ref_young)+1):
    if i < young_cutoff:
        ko_wt_young.append(str(0)) #str('control'))
    if i > young_cutoff:
        ko_wt_young.append(str(1)) #str('knockout'))

print(len(ko_wt_young))
mouse_id_young = []
track = 0
for i in range(size_set_young):
    for j in range(corpus_index_ref_young[i+1]):
        mouse_id_young.append(str(i+1))
        track+=1
print(len(mouse_id_young))
print(output_young.shape)
data_annot = name_annot.to_list()
data_annot.insert(0, 'condition')
data_annot.insert(0, 'mouse')
output_young.loc[-1] = ko_wt_young
output_young.loc[-2] = mouse_id_young
output_young.index = output_young.index + 2
output_young = output_young.sort_index()
output_young.index = data_annot

corpus_out_young = output_young.transpose()
corpus_out_young.to_csv('corpus_data_young.csv')
