import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues
import matplotlib.pyplot as plt
import seaborn as sns

#--- Formatting
sns.set_context('talk', font_scale=1.5)
palette = {'Control': 'darkgrey', 'GLUT1cKO': 'teal'}

#--- Load and Clean
def load_and_clean_csv(path, delimiter=';', comment='#'):
    df = pd.read_csv(path, delimiter=delimiter, comment=comment)
    return df.apply(pd.to_numeric, errors='coerce')  # Convert all to numeric

#--- Permutation Test
def permutation_test(x, y, n_resamples=10000, statistic='mean', subset_size=100, seed=None):
    rng = np.random.default_rng(seed)

    x = np.asarray(x)
    y = np.asarray(y)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]

    # Compute observed difference
    observed = np.mean(x) - np.mean(y) if statistic == 'mean' else np.median(x) - np.median(y)

    # Combine all for permutation
    combined = np.concatenate([x, y])
    n_total = len(combined)

    # Subsampling size must be ≤ group size
    subset_size = min(subset_size, len(x), len(y))

    resampled_stats = []
    for _ in range(n_resamples):
        permuted = rng.permutation(combined)
        x_perm = permuted[:subset_size]
        y_perm = permuted[subset_size:2*subset_size]
        stat = np.mean(x_perm) - np.mean(y_perm) if statistic == 'mean' else np.median(x_perm) - np.median(y_perm)
        resampled_stats.append(stat)

    p_value = np.mean(np.abs(resampled_stats) >= np.abs(observed))
    return observed, p_value

def combine_pvalues_harmonic(pvals):
    pvals = np.array(pvals)
    return len(pvals) / np.sum(1.0 / pvals)


#--- Normalize Entire Batch by Max
def normalize_batch(dataframes):
    all_vals = pd.concat([df.select_dtypes(include=[np.number]) for df in dataframes])
    max_val = np.nanmax(all_vals.values)
    return [df / max_val for df in dataframes], max_val

import os

#--- KDE Plotting Function
def plot_kde_comparison(wt_vals, ko_vals, wt_id, ko_id, batch_name, lipid_name):
    fig, ax = plt.subplots(figsize=(6, 4))

    sns.kdeplot(wt_vals, ax=ax, label=f'WT {wt_id}', color='darkgrey', linewidth=2)
    sns.kdeplot(ko_vals, ax=ax, label=f'KO {ko_id}', color='teal', linewidth=2)

    ax.set_xlabel('Normalized Ion Intensity')
    ax.set_ylabel('Density')
    ax.set_title(f'{lipid_name}\n{batch_name}: WT {wt_id} vs KO {ko_id}')
    ax.legend()

    sns.despine()
    plt.tight_layout()

    # Ensure output directory exists
    os.makedirs('fig_output', exist_ok=True)
    safe_lipid = lipid_name.replace(':', '-').replace('/', '-')
    fname = f'fig_output/{batch_name.replace(" ", "")}_WT{wt_id}_KO{ko_id}_{safe_lipid}.png'
    fig.savefig(fname)
    plt.close(fig)

#--- Lipid Name
lipid_name = 'SHexCer 36:1;O2'
lipid_idx = 81

# lipid_name = 'PI 36:4, PG 38:5;O'
# lipid_idx = 86

#--- Define batches and subjects
batches = [
    {
        'name': 'Batch 1',
        'WT': {
            '565': 'ion_counts/565-cortex-Total Ion Count.csv',
        },
        'KO': {
            '271': 'ion_counts/cortex-271-Total Ion Count.csv',
            '373': 'ion_counts/373-cortex-Total Ion Count.csv',
        }
    },
    {
        'name': 'Batch 2',
        'WT': {
            '163': 'ion_counts/163-cortex-Total Ion Count.csv',
            '164': 'ion_counts/164-cortex-Total Ion Count.csv',
            '162': 'ion_counts/162-cortex-Total Ion Count.csv',
        },
        'KO': {
            '272': 'ion_counts/272-cortex-Total Ion Count.csv',
            '292': 'ion_counts/292-cortex-Total Ion Count.csv',
        }
    }
]

#--- Process batches
batch_results = []
ctrl_all, ko_all = [], []

#--- Process batches
batch_results = []
ctrl_all, ko_all = [], []

for batch in batches:
    batch_name = batch['name']

    # Load all WT and KO data first
    wt_ids, wt_paths = zip(*batch['WT'].items())
    ko_ids, ko_paths = zip(*batch['KO'].items())

    wt_dfs = [load_and_clean_csv(path) for path in wt_paths]
    ko_dfs = [load_and_clean_csv(path) for path in ko_paths]

    # Normalize all samples in batch by the global max of this batch
    all_dfs = wt_dfs + ko_dfs
    norm_dfs, batch_max = normalize_batch(all_dfs)
    wt_dfs = norm_dfs[:len(wt_dfs)]
    ko_dfs = norm_dfs[len(wt_dfs):]

    # Pairwise permutation tests
    for i, (wt_id, wt_df) in enumerate(zip(wt_ids, wt_dfs)):
        for j, (ko_id, ko_df) in enumerate(zip(ko_ids, ko_dfs)):
            wt_vals = wt_df.iloc[lipid_idx, :].to_numpy()
            ko_vals = ko_df.iloc[lipid_idx, :].to_numpy()

            # --- Print means prior to permutation test
            wt_mean = np.nanmean(wt_vals)
            ko_mean = np.nanmean(ko_vals)
            print(f"{batch_name} | WT {wt_id} mean: {wt_mean:.4f}, KO {ko_id} mean: {ko_mean:.4f}")

            obs_diff, pval = permutation_test(wt_vals, ko_vals)
            plot_kde_comparison(wt_vals, ko_vals, wt_id, ko_id, batch_name, lipid_name)

            batch_results.append({
                'batch': batch_name,
                'wt_id': wt_id,
                'ko_id': ko_id,
                'obs_diff': obs_diff,
                'p_value': pval
            })

            ctrl_all.append(wt_vals)
            ko_all.append(ko_vals)

#--- Combine p-values
pvals = [res['p_value'] for res in batch_results]
combined_p = combine_pvalues_harmonic(pvals)

#--- Print Results
print("\nPairwise Results:")
for res in batch_results:
    print(f"{res['batch']} | WT {res['wt_id']} vs KO {res['ko_id']}: diff = {res['obs_diff']:.4f}, p = {res['p_value']:.4g}")
print(f"\nCombined p-value (Fisher’s method): {combined_p:.4g}")

#--- Plot Summary
ctrl_means = np.array([np.nanmean(arr) for arr in ctrl_all])
ko_means = np.array([np.nanmean(arr) for arr in ko_all])

means = [np.mean(ctrl_means), np.mean(ko_means)]
stds = [np.std(ctrl_means), np.std(ko_means)]

group_labels = ['Control'] * len(ctrl_means) + ['GLUT1cKO'] * len(ko_means)
subject_vals = np.concatenate([ctrl_means, ko_means])

fig, ax = plt.subplots(figsize=(6, 5))

sns.barplot(x=['Control', 'GLUT1cKO'], y=means, ax=ax, palette=palette,
            errorbar=None, edgecolor='black', linewidth=3, zorder=1)
ax.errorbar(x=[0, 1], y=means, yerr=stds, fmt='none', ecolor='black',
            elinewidth=1.5, capsize=5, zorder=2)

sns.stripplot(x=group_labels, y=subject_vals, ax=ax, palette=palette,
              size=8, jitter=0.15, edgecolor='black', linewidth=0.5, zorder=3)

ax.set_ylabel(f'Normalized Ion Intensity\n({lipid_name})')
ax.set_ylim(0, max(means) + max(stds) + 0.1)
ax.text(0.3, max(means) + 0.05, f'p = {combined_p:.4f}', fontsize=14, weight='bold')

sns.despine()
plt.tight_layout()

import re

safe_lipid_name = re.sub(r'[\\/*?:"<>|]', "_", lipid_name)
plt.savefig(f'fig_output/{safe_lipid_name}_barplot.png')

#--- Save
result_df = pd.DataFrame(batch_results)
result_df['lipid'] = lipid_name
result_df.to_csv('permutation_test_results.csv', index=False)