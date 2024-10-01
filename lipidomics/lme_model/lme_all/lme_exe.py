# # Libraries
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import numpy as np
import pandas as pd

import statsmodels.api as sm
import statsmodels.formula.api as smf

import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
warnings.simplefilter('ignore', ConvergenceWarning)

# # Formatting
# Seaborn (matplotlib)
sns.set_theme(style='ticks', font_scale=2)

# Reference
ref = pd.read_csv('lipids-correlation-fullset-youngmice.csv', delimiter=',')
grp_annot = ref.iloc[:,8]
name_annot = ref.iloc[:,7]

input_old = pd.read_csv('corpus_data_old.csv')
print(input_old)
input_young = pd.read_csv('corpus_data_young.csv')
print(input_young)

def LME_model(data_in):
    out_list = []
    llr_test = 0
    output1, output2 = np.zeros(data_in.shape[1]-3), np.zeros(data_in.shape[1]-3)
    for i in range(data_in.shape[1]-3):
        df = pd.concat([data_in.iloc[:,:3], data_in.iloc[:,i+3]], axis=1)
        df.columns = ['element', 'mouse', 'condition', 'lipid']
        model_with_random = smf.mixedlm("lipid ~ condition", data=df, groups="mouse", re_formula="~(element+mouse)").fit()
        #model_without = smf.ols("lipid ~ condition", df).fit()

        print(model_with_random.summary())
        print(name_annot[i])

        output1[i] = model_with_random.pvalues.iloc[1]
        output2[i] = np.mean(model_with_random.pvalues.iloc[:])


        from scipy.stats import chi2

        # ll_rand = model_with_random.llf
        # ll_norm = model_without.llf
        # lr_stat = 2 * (ll_rand - ll_norm)
        # df_diff = (model_with_random.df_modelwc - model_without.df_model)
        # p_value = chi2.sf(lr_stat, df_diff)
        # llr_test += p_value

        if (output1[i] < 0.05):
            out_list.append(name_annot[i])
    print(llr_test / (i+1))

    return output1, output2, out_list

def compute_differential(data_in):
    output = np.zeros(data_in.shape[1]-3)
    indexer = 0
    for j in range(data_in.shape[0]):
        if data_in.iloc[j]['condition'] == 0:
            indexer+=1
    for i in range(data_in.shape[1]-3):
        mean_in = np.mean(data_in.iloc[indexer:, i+3])
        mean_ref = np.mean(data_in.iloc[0:indexer, i+3])
        output[i] = (mean_in - mean_ref) #/ mean_ref
    return output



# -------------------------- Example code - Ventricle ----------------------#
pval1_old, pval2_old, sig_old = LME_model(input_old)
diff_old = compute_differential(input_old)
print(sig_old)
pval1_young, pval2_young, sig_young = LME_model(input_young)
diff_young = compute_differential(input_young)
print(sig_young)


fig, ax = plt.subplots(layout='constrained')
ax.scatter(x=diff_old, y=-np.log10(pval1_old), s=16, label=None, color='darkred', alpha=0.8)
# ax.scatter(x=diff_old, y=-np.log10(pval2_old), s=16, label=None, color='red')
ax.scatter(x=diff_old, y=-np.log10(pval1_young), s=16, label=None, color='darkgray', alpha=0.6)

plt.axvline(-0.1,color="grey",linestyle="--")
plt.axvline(0.1,color="grey",linestyle="--")
plt.axhline(1.3,color="grey",linestyle="--")

# for i in range(len(diff_old)):
#     if -np.log10(pval1_old[i]) > 1.3:
#         plt.text(x=diff_old[i], y=-np.log10(pval1_old[i]), s=name_annot[i], fontsize=9)

plt.show()