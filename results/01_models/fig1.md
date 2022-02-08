---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Figure 1

## Method

1. Train 13 (species) x 5 (sample size) x 2 (within and across) models

```{code-cell} ipython3
from collections import OrderedDict
from glob import glob
import hashlib
from io import StringIO
import itertools
import json
from pathlib import Path

from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyfaidx import Fasta
import scipy.optimize
import scipy.stats
import sklearn
import tensorflow as tf

from a2z.data import IntervalReferenceDataset
from a2z.models import allow_growth

allow_growth()

clades = OrderedDict([
    ('Vitis vinifera', {'Angiosperm', 'Dicot'}),
    ('Arabidopsis thaliana', {'Angiosperm', 'Dicot'}),
    ('Eutrema salsugineum', {'Angiosperm', 'Dicot'}),
    ('Populus trichocarpa', {'Angiosperm', 'Dicot'}),
    ('Phaseolus vulgaris', {'Angiosperm', 'Dicot'}),
    ('Glycine max', {'Angiosperm', 'Dicot'}),
    ('Asparagus officinalis', {'Angiosperm', 'Monocot'}),
    ('Brachypodium distachyon', {'Angiosperm', 'Monocot', 'Grass'}),
    ('Hordeum vulgare', {'Angiosperm', 'Monocot', 'Grass'}),
    ('Oryza sativa', {'Angiosperm', 'Monocot', 'Grass'}),
    ('Setaria viridis', {'Angiosperm', 'Monocot', 'Grass'}),
    ('Sorghum bicolor', {'Angiosperm', 'Monocot', 'Grass'}),
    ('Zea mays', {'Angiosperm', 'Monocot', 'Grass'})
])
```

```{code-cell} ipython3
results_acc = {
    'species': [],
    'model': [],
    'auPR': [],
    'auROC': [],
    'f1_score': []
}

configuration_across = "accessibility_base"
configuration_within = "accessibility_base_within"
for metrics_path in map(lambda p: Path(p), glob(f"tmp/results/{configuration_across}/*/*/metrics.json")):
    with open(metrics_path) as metrics_file:
        try:
            metrics_json = json.load(metrics_file)
        except json.JSONDecodeError:
            continue
    results_acc['species'].append(metrics_path.parts[3].replace("_", " "))
    results_acc['model'].append('across')
    results_acc['auPR'].append(metrics_json['auPR'])
    results_acc['auROC'].append(metrics_json['auROC'])
    results_acc['f1_score'].append(metrics_json['f1_score'])
    
for metrics_path in map(lambda p: Path(p), glob(f"tmp/results/{configuration_within}/*/*/metrics.json")):
    with open(metrics_path) as metrics_file:
        try:
            metrics_json = json.load(metrics_file)
        except json.JSONDecodeError:
            continue
    results_acc['species'].append(metrics_path.parts[3].replace("_", " "))
    results_acc['model'].append('within')
    results_acc['auPR'].append(metrics_json['auPR'])
    results_acc['auROC'].append(metrics_json['auROC'])
    results_acc['f1_score'].append(metrics_json['f1_score'])

results_acc = pd.DataFrame(results_acc)
```

```{code-cell} ipython3
def shorten_species(scientific_name: str) -> str:
    genus, species = scientific_name.split(" ")
    return f"{genus[0]}. {species}"
```

```{code-cell} ipython3
def plot_results_bar(results: pd.DataFrame, ax: matplotlib.axes.Axes, label_y = True, label_all = False):
    results_agg = results.groupby(['species', 'model'])['auPR'].aggregate([np.mean, scipy.stats.sem]).reindex(clades, level = 'species')

    species = results_agg.index.get_level_values('species').unique()

    clade_list = list(clades.keys())
    y_pos = np.array([clade_list.index(s) for s in species])
    width = 0.35
    ax.barh(y_pos + (width / 2) + 1, results_agg.xs('across', level = 'model')['mean'], xerr = results_agg.xs('across', level = 'model')['sem'], height = width, label = 'Across')
    ax.barh(y_pos - (width / 2) + 1, results_agg.xs('within', level = 'model')['mean'], xerr = results_agg.xs('within', level = 'model')['sem'], height = width, label = 'Within')
    
    ax.set_yticks(y_pos + 1)
    ax.set_yticklabels([shorten_species(s) for s in clade_list if s in species])
    ax.set_xlabel('auPR')
    ax.set_xlim(0, 1)
    ax.set_ylim(1 - width, y_pos.max() + width + 1)
    ax.legend()
    # flip the axes
    ax.invert_yaxis()
    
fig, ax = plt.subplots()
plot_results_bar(results_acc, ax)
```

```{code-cell} ipython3
results_hypo = {
    'species': [],
    'model': [],
    'auPR': [],
    'auROC': [],
    'f1_score': []
}

configuration_across = "hypomethylation_base"
configuration_within = "hypomethylation_base_within"
for metrics_path in map(lambda p: Path(p), glob(f"tmp/results/{configuration_across}/*/*/metrics.json")):
    with open(metrics_path) as metrics_file:
        try:
            metrics_json = json.load(metrics_file)
        except json.JSONDecodeError:
            continue
    results_hypo['species'].append(metrics_path.parts[3].replace("_", " "))
    results_hypo['model'].append('across')
    results_hypo['auPR'].append(metrics_json['auPR'])
    results_hypo['auROC'].append(metrics_json['auROC'])
    results_hypo['f1_score'].append(metrics_json['f1_score'])
    
for metrics_path in map(lambda p: Path(p), glob(f"tmp/results/{configuration_within}/*/*/metrics.json")):
    with open(metrics_path) as metrics_file:
        try:
            metrics_json = json.load(metrics_file)
        except json.JSONDecodeError:
            continue
    results_hypo['species'].append(metrics_path.parts[3].replace("_", " "))
    results_hypo['model'].append('within')
    results_hypo['auPR'].append(metrics_json['auPR'])
    results_hypo['auROC'].append(metrics_json['auROC'])
    results_hypo['f1_score'].append(metrics_json['f1_score'])

results_hypo = pd.DataFrame(results_hypo)
```

```{code-cell} ipython3
fig, ax = plt.subplots()
plot_results_bar(results_hypo, ax)
```

```{code-cell} ipython3
tree_newick = StringIO("((Vitis_vinifera:100,((Arabidopsis_thaliana:30,Eutrema_salsugineum:30):60,Populus_trichocarpa:90,(Phaseolus_vulgaris:25,Glycine_max:25):65):10):60,(Asparagus_officinalis:110,((Brachypodium_distachyon:35,Hordeum_vulgare:35):15,Oryza_sativa:50,(Setaria_viridis:25,(Sorghum_bicolor:10,Zea_mays:10):15):25):60):50):3;")
tree = Phylo.read(tree_newick, "newick")

def plot_tree(tree: Phylo.Newick.Tree ,ax: matplotlib.axes.Axes):
    Phylo.draw(
        tree = tree,
        label_func = lambda clade: None,
        do_show = False,
        axes = ax
    )
    ax.set_xlabel("MYA")
    ax.set_xticks([3, 43, 83, 123, 163])
    ax.set_xticklabels([160, 120, 80, 40, 0])
    ax.set_xlim(0, 163)
    ax.set_yticks(range(1, 13))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_visible(False)

fig, ax = plt.subplots()
plot_tree(tree, ax)
```

```{code-cell} ipython3
def plot_genome_size_performance(ax: matplotlib.axes.Axes, configuration_across: str, results: pd.DataFrame):
    results_agg = results.groupby(['species', 'model'])['auPR'].aggregate([np.mean, scipy.stats.sem]).reindex(clades, level = 'species')

    with open(f"configs/{configuration_across}.json") as config_file:
        config = json.load(config_file)

    genome_dir = Path(config['paths']['genomes'])
    references = {}
    for reference_path in filter(lambda p: p.is_file() and (p.suffix == '.fa'),genome_dir.iterdir()):
        species = reference_path.stem.replace('_', ' ')
        references[species] = Fasta(str(reference_path), as_raw = True, rebuild = False)

    genome_sizes = {species: sum([len(seq) for seq in references[species]]) for species in references}

    x = np.array([genome_sizes[species] for species in references])
    y_across = np.array([results_agg.loc[(species, 'across'), 'mean'] for species in references])
    y_within = np.array([results_agg.loc[(species, 'within'), 'mean'] for species in references])

    x1, b = np.polyfit(np.log(x), y_across, 1)

    sct_across = ax.scatter(x = x, y = y_across, label = "Across")
    sct_within = ax.scatter(x = x, y = y_within, label = "Within")
    xx = np.arange(10_000_000, max(genome_sizes.values()) + 500_000_000, step = 10_000_000)
    ax.plot(xx, ((x1 * np.log(xx)) + b), linestyle = "--", color = sct_across.get_facecolor())
    ax.set_ylim(0, 1)
    ax.set_ylabel('auPR')
    ax.set_xlim(0, max(genome_sizes.values()) + 100_000_000)
    ax.set_xlabel('Genome Size')
    ax.legend()
```

```{code-cell} ipython3
fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 5), constrained_layout = True)
plot_genome_size_performance(axs[0], 'accessibility_base', results_acc)
axs[0].set_title("Accessibility")
plot_genome_size_performance(axs[1], 'hypomethylation_base', results_hypo)
axs[1].set_title("Methylation")
fig.savefig("figs/genome-size-auPR.png", dpi = 300, facecolor = "white")
```

```{code-cell} ipython3
def predict_test(configuration_across: str) -> pd.DataFrame:
    with open(f"configs/{configuration_across}.json") as config_file:
        config = json.load(config_file)

    split_hash = hashlib.sha256((json.dumps(config['paths'], indent = 2) + "\n" + json.dumps(config['preprocessing'], indent = 2) + "\n" + json.dumps(config['splitting'], indent = 2) + "\n").encode('utf-8')).hexdigest()
    config_hash = hashlib.sha256((json.dumps(config, indent = 2) + "\n").encode('utf-8')).hexdigest()
    split_data_path = f"tmp/data/split/{split_hash}.tsv"

    intervals = pd.read_table(
        split_data_path,
        header = 0,
        dtype = {
            'seqid': str,
            'start': int,
            'end': int,
            'species': str,
            'target': np.byte,
            'is_train': np.byte,
            'is_test': np.byte,
            'is_val': np.byte
        }
    )
    
    sp = intervals['species'].unique()

    references = {}
    for species in sp:
        reference_path = Path(config['paths']['genomes']) / f"{species.replace(' ', '_')}.fa"
        references[species] = Fasta(str(reference_path), as_raw = True, rebuild = False)

    test_intervals_sp = []
    for test_species in sp:
        test_intervals = intervals.loc[(intervals['species'] == test_species) & (intervals['is_test'] == 1), :].copy()
        test_data = IntervalReferenceDataset(test_intervals, references)
        model_path = f"tmp/results/{configuration_across}/{test_species.replace(' ', '_')}/0/model"
        model = tf.keras.models.load_model(model_path)
        predictions = model.predict(test_data)
        test_intervals['prediction'] = predictions
        test_intervals_sp.append(test_intervals)
    test_intervals = pd.concat(test_intervals_sp).reset_index(drop = True)
    del test_intervals_sp
    
    return test_intervals
```

```{code-cell} ipython3
%%time

cache_path = Path("tmp/accessibility_base_preds.tsv")
if not cache_path.exists():
    # ~5 mins
    test_intervals_acc = predict_test('accessibility_base')
    test_intervals_acc.to_csv(
        cache_path,
        sep = "\t",
        index = False
    )
else:
    test_intervals_acc = pd.read_table(
        cache_path,
        header = 0,
        dtype = {
            'seqid': 'category',
            'target': np.byte,
            'distance_class': 'category',
            'species': 'category',
            'is_test': np.byte,
            'is_val': np.byte,
            'is_train': np.byte
        }
    )
```

```{code-cell} ipython3
%%time

cache_path = Path("tmp/hypomethylation_base_preds.tsv")
if not cache_path.exists():
    # ~4 mins
    test_intervals_hypo = predict_test('hypomethylation_base')
    test_intervals_hypo.to_csv(
        cache_path,
        sep = "\t",
        index = False
    )
else:
    test_intervals_hypo = pd.read_table(
        cache_path,
        header = 0,
        dtype = {
            'seqid': 'category',
            'target': np.byte,
            'distance_class': 'category',
            'species': 'category',
            'is_test': np.byte,
            'is_val': np.byte,
            'is_train': np.byte
        }
    )
```

```{code-cell} ipython3
def plot_distance_class_pr(test_intervals: pd.DataFrame, ax: matplotlib.axes.Axes, linestyle: str, colors: tuple[str, str, str]):
    precision = {}
    recall = {}
    baseline = {}
    thresholds = {}

    for gene_distance_class in ['distal', 'proximal', 'genic']:
        # false positives are double-counted
        intervals = test_intervals.loc[(test_intervals['distance_class'] == gene_distance_class), :]
        y = intervals['target']
        y_hat = intervals['prediction']

        baseline[gene_distance_class] = y.mean()

        precision[gene_distance_class], recall[gene_distance_class], thresholds[gene_distance_class] = sklearn.metrics.precision_recall_curve(y, y_hat)

    ax_line_pr_distal = ax.plot(recall['distal'], precision['distal'], linestyle = linestyle, color = colors[0], label = 'Distal')
    ax_line_pr_proximal = ax.plot(recall['proximal'], precision['proximal'], linestyle = linestyle, color = colors[1], label = 'Proximal')
    ax_line_pr_genic = ax.plot(recall['genic'], precision['genic'], linestyle = linestyle, color = colors[2], label = 'Genic')
    ax.set_xlim(0, 1)
    ax.set_xlabel('Recall')
    ax.set_ylim(0, 1)
    ax.set_ylabel('Precision')
```

```{code-cell} ipython3
fig, ax = plt.subplots()
distance_class_colors = ['C0', 'C1', 'C2']
plot_distance_class_pr(test_intervals_acc, ax, linestyle = '-', colors = distance_class_colors)
plot_distance_class_pr(test_intervals_hypo, ax, linestyle = '--', colors = distance_class_colors)

ax.legend(handles = [
    matplotlib.patches.Patch(color = distance_class_colors[0], label = 'Distal'),
    matplotlib.patches.Patch(color = distance_class_colors[1], label = 'Proximal'),
    matplotlib.patches.Patch(color = distance_class_colors[2], label = 'Genic'),
    matplotlib.lines.Line2D(xdata = [], ydata = [], color = 'black', linestyle = '-', label = 'Accessibility'),
    matplotlib.lines.Line2D(xdata = [], ydata = [], color = 'black', linestyle = '--', label = 'Hypomethylation')
])
```

```{code-cell} ipython3
def get_distance_class_results(test_intervals: pd.DataFrame, threshold: float = 0.9):
    class_results = {
        'species': [],
        'distance_class': [],
        'TP': [],
        'FP': [],
        'FN': [],
        'TN': []
    }
    for (species, gene_distance_class), data in test_intervals.groupby(['species', 'distance_class']):
        class_results['species'].append(species)
        class_results['distance_class'].append(gene_distance_class)
        class_results['TP'].append(len(data.loc[(data['target'] == 1) & (data['prediction'] >= threshold)]))
        class_results['FP'].append(len(data.loc[(data['target'] == 0) & (data['prediction'] >= threshold)]))
        class_results['FN'].append(len(data.loc[(data['target'] == 1) & (data['prediction'] < threshold)]))
        class_results['TN'].append(len(data.loc[(data['target'] == 0) & (data['prediction'] < threshold)]))
    class_results = pd.DataFrame(class_results)
    class_results['recall'] = class_results['TP'] / (class_results['TP'] + class_results['FN'])
    
    return class_results
```

```{code-cell} ipython3
class_results_acc = get_distance_class_results(test_intervals_acc)
```

```{code-cell} ipython3
class_results_hypo = get_distance_class_results(test_intervals_hypo)
```

```{code-cell} ipython3
def plot_distance_class_species_count(ax: matplotlib.axes.Axes, class_results: pd.DataFrame):
    species = [s for s in clades if s in class_results['species'].unique()]
    x_pos = np.arange(len(species))
    width = 0.25
    ax.bar(x_pos - width, height = class_results.loc[class_results['distance_class'] == 'distal'].set_index('species').reindex(clades).dropna().apply(lambda row: row['TP'] + row['FN'], axis = 1), width = width, label = 'Distal')
    ax.bar(x_pos, height = class_results.loc[class_results['distance_class'] == 'proximal'].set_index('species').reindex(clades).dropna().apply(lambda row: row['TP'] + row['FN'], axis = 1), width = width, label = 'Proximal')
    ax.bar(x_pos + width, height = class_results.loc[class_results['distance_class'] == 'genic'].set_index('species').reindex(clades).dropna().apply(lambda row: row['TP'] + row['FN'], axis = 1), width = width, label = 'Genic')
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(x_pos))
    ax.set_xticklabels(list(map(lambda s: ''.join([e[0] for e in s.split(" ")]), species)))
    ax.set_xlabel('Species')
    ax.set_ylabel('Count')
    ax.legend()
```

```{code-cell} ipython3
fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 5))
plot_distance_class_species_count(axs[0], class_results_acc)
axs[0].set_title("Accessibility")
plot_distance_class_species_count(axs[1], class_results_hypo)
axs[1].set_title("Methylation")
fig.savefig("figs/distance_class_species_count.png", dpi = 300, facecolor = "white")
```

```{code-cell} ipython3
def plot_fdr_vs_for(ax: matplotlib.axes.Axes, class_results: pd.DataFrame):
    species_results = class_results.groupby('species').aggregate(sum).drop(columns = ['recall'])

    species = [s for s in clades if s in species_results.index.unique()]
    x_pos = np.arange(len(species))
    width = 0.25
    ax.bar(x_pos - (width / 2), height = species_results['FP'] / (species_results['TP'] + species_results['FP']), width = width, label = 'FDR')
    ax.bar(x_pos + (width / 2), height = species_results['FN'] / (species_results['TN'] + species_results['FN']), width = width, label = 'FOR')
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(x_pos))
    ax.set_xticklabels(list(map(lambda s: ''.join([e[0] for e in s.split(" ")]), species)))
    ax.set_xlabel('Species')
    ax.set_ylim(0, 1)
    ax.legend()
```

```{code-cell} ipython3
fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 5))
plot_fdr_vs_for(axs[0], class_results_acc)
axs[0].set_title("Accessibility")
plot_fdr_vs_for(axs[1], class_results_hypo)
axs[1].set_title("Methylation")
fig.savefig("figs/fdr_vs_for.png", dpi = 300, facecolor = "white")
```

```{code-cell} ipython3
def plot_pr(ax: matplotlib.axes.Axes, test_intervals: pd.DataFrame):
    ax.set_prop_cycle('color',[plt.cm.tab20(i) for i in range(len(plt.cm.tab20.colors))])
    for species, data in test_intervals.groupby('species'):
        y = data['target']
        y_hat = data['prediction']

        baseline = y.mean()
        precision, recall, _ = sklearn.metrics.precision_recall_curve(y, y_hat)
        aupr = sklearn.metrics.average_precision_score(y, y_hat)
        ax.plot(recall, precision, label = f"{''.join((w[0] for w in species.split(' ')))} ({aupr:.2f})")
    ax.set_xlim(0, 1)
    ax.set_xlabel('Recall')
    ax.set_ylim(0, 1)
    ax.set_ylabel('Precision')
    ax.legend(ncol = 2)
```

```{code-cell} ipython3
def plot_roc(ax: matplotlib.axes.Axes, test_intervals: pd.DataFrame):
    ax.set_prop_cycle('color',[plt.cm.tab20(i) for i in range(len(plt.cm.tab20.colors))])
    for species, data in test_intervals.groupby('species'):
        y = data['target']
        y_hat = data['prediction']

        fpr, tpr, _ = sklearn.metrics.roc_curve(y, y_hat)
        auroc = sklearn.metrics.roc_auc_score(y, y_hat)
        ax.plot(fpr, tpr, label = f"{shorten_species(species)} ({auroc:.2f})")
    ax.plot([0, 1], [0, 1], color="gray", linestyle="--")
    ax.set_xlim(0, 1)
    ax.set_xlabel('FPR')
    ax.set_ylim(0, 1)
    ax.set_ylabel('TPR')
    ax.legend()
```

```{code-cell} ipython3
fig, axs = plt.subplots(constrained_layout = True, nrows = 1, ncols = 2, figsize = (12, 5))
plot_roc(axs[0], test_intervals_acc)
axs[0].set_title("Accessibility")
plot_roc(axs[1], test_intervals_hypo)
axs[1].set_title("Methylation")
fig.savefig("figs/roc.png", dpi = 300, facecolor = "white")
```

```{code-cell} ipython3
fig1 = plt.figure(figsize = (18, 10), constrained_layout = True)

gs_outer = matplotlib.gridspec.GridSpec(2, 1, figure = fig1)
gs_inner_top = matplotlib.gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec = gs_outer[0, 0])
gs_inner_bot = matplotlib.gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec = gs_outer[1, 0])
ax1a = fig1.add_subplot(gs_inner_top[0, 0])
ax1b = fig1.add_subplot(gs_inner_top[0, 1], sharey = ax1a)
ax1c = fig1.add_subplot(gs_inner_top[0, 2])
ax1d = fig1.add_subplot(gs_inner_bot[0, 0])
ax1e = fig1.add_subplot(gs_inner_bot[0, 1])
ax1f = fig1.add_subplot(gs_inner_bot[0, 2])

# Fig 1A: Tree
plot_tree(tree, ax1a)

# Fig 1B: auPR for accessibility
plot_results_bar(results_acc, ax1b)
ax1b.set_title("Accessibility")

# Fig 1C: auPR for hypomethylation
plot_results_bar(results_hypo, ax1c)
ax1c.set_title("Methylation")

# Fig 1F: Distance class recall
plot_distance_class_pr(test_intervals_acc, ax1d, linestyle = '-', colors = distance_class_colors)
plot_distance_class_pr(test_intervals_hypo, ax1d, linestyle = '--', colors = distance_class_colors)

ax1d.legend(handles = [
    matplotlib.patches.Patch(color = distance_class_colors[0], label = 'Distal'),
    matplotlib.patches.Patch(color = distance_class_colors[1], label = 'Proximal'),
    matplotlib.patches.Patch(color = distance_class_colors[2], label = 'Genic'),
    matplotlib.lines.Line2D(xdata = [], ydata = [], color = 'black', linestyle = '-', label = 'Accessibility'),
    matplotlib.lines.Line2D(xdata = [], ydata = [], color = 'black', linestyle = '--', label = 'Methylation')
])

# Fig 1E/F: PR curves for accessibility and hypomethylation
plot_pr(ax1e, test_intervals_acc)
plot_pr(ax1f, test_intervals_hypo)

fig1.savefig("figs/fig1.png", dpi = 300, facecolor = 'white')
```

```{code-cell} ipython3

```
