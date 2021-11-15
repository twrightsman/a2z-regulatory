---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.0
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

```{code-cell} ipython3
from collections import OrderedDict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
```

```{code-cell} ipython3
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
JASPAR_pGIA = pd.read_table(
    "tmp/JASPAR_pGIA.all.melt.tsv",
    index_col = [0, 1, 2, 3, 4]
)
```

```{code-cell} ipython3
top_jaspar_max = JASPAR_pGIA.groupby(
    ["model", "species", "motif"]
).aggregate(
    max
).reset_index(
    level = ['model', 'species']
).sort_values(
    by = ["model", "species", "global_importance"],
    ascending = False
).set_index(
    ["model", "species"],
    append = True
).reorder_levels(
    ["model", "species", "motif"]
).groupby(
    ["model", "species"]
).head(3)
```

```{code-cell} ipython3
medoids_embedded = pd.read_table(
    "tmp/medoids_embedded.tsv"
)
```

```{code-cell} ipython3
medoids_embedded.head()
```

# Plotting

```{code-cell} ipython3
def plot_pGIA(model, species, ax, legend_motif_pos = None, legend_strand_pos = None):
    data = top_jaspar_max.xs((model, species), level = ['model', 'species'])
    colors = []
    labels = []
    for motif in data.index.get_level_values('motif'):
        pGIA_pos = JASPAR_pGIA.xs(key = (model, species, motif, '+'), level = ["model", "species", "motif", "strand"])
        pGIA_pos_color = ax.plot(pGIA_pos, alpha = 0.4)[0].get_color()
        colors.append(pGIA_pos_color)
        labels.append(motif)
        pGIA_neg = JASPAR_pGIA.xs(key = (model, species, motif, '-'), level = ["model", "species", "motif", "strand"])
        ax.plot(pGIA_neg, linestyle = 'dotted', color = pGIA_pos_color, alpha = 0.4)

    ax.set_xlabel("Position")
    ax.set_ylabel("Global Importance")
    legend_motifs = ax.legend(
        loc = "best" if legend_motif_pos is None else legend_motif_pos,
        handles = [
            matplotlib.lines.Line2D(xdata = [], ydata = [], color = colors[0], label = labels[0]),
            matplotlib.lines.Line2D(xdata = [], ydata = [], color = colors[1], label = labels[1]),
            matplotlib.lines.Line2D(xdata = [], ydata = [], color = colors[2], label = labels[2])
        ]
    )
    ax.add_artist(legend_motifs)
    if legend_strand_pos is not None:
        legend_strand = ax.legend(
            loc = legend_strand_pos,
            handles = [
                matplotlib.lines.Line2D(xdata = [], ydata = [], color = 'black', linestyle = '-', label = 'Forward'),
                matplotlib.lines.Line2D(xdata = [], ydata = [], color = 'black', linestyle = 'dotted', label = 'Reverse')
            ]
        )
        ax.add_artist(legend_strand)
```

```{code-cell} ipython3
fig, ax = plt.subplots()
plot_pGIA("accessibility_base", "Arabidopsis thaliana", ax)
```

```{code-cell} ipython3
def plot_embedding_by_config(ax):
    labels = {
        'accessibility_base': 'Accessibility',
        'hypomethylation_base': 'Hypomethylation'
    }

    for config_name, data in medoids_embedded.groupby('config'):
        ax.scatter(data['embedded1'], data['embedded2'], alpha = 0.5, label = labels[config_name])
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlabel("Dimension 1")
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_ylabel("Dimension 2")
    ax.legend()
```

```{code-cell} ipython3
def plot_embedding_by_species(ax, config: str):
    for species, data in medoids_embedded.loc[medoids_embedded['config'] == config].groupby('species'):
        ax.scatter(data['embedded1'], data['embedded2'], alpha = 0.5, label = shorten_species(species))
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlabel("Dimension 1")
    ax.set_xlim(xlim[0], 2. * xlim[1])
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_ylabel("Dimension 2")
    ax.set_ylim(ylim[0], 2. * ylim[1])
    ax.legend()
```

```{code-cell} ipython3
def plot_embedding_by_clade(ax, config: str):
    for clade, color in (('Monocot', 'goldenrod'), ('Dicot', 'darkgreen')):
        clade_set = set(filter(lambda s: clade in clades[s], clades.keys()))
        data = medoids_embedded.loc[medoids_embedded['species'].isin(clade_set) & (medoids_embedded['config'] == config)]
        ax.scatter(data['embedded1'], data['embedded2'], alpha = 0.5, label = clade, color = color)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlabel("Dimension 1")
    #ax.set_xlim(xlim[0], 2. * xlim[1])
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_ylabel("Dimension 2")
    #ax.set_ylim(ylim[0], 2. * ylim[1])
    ax.legend()
```

```{code-cell} ipython3
def shorten_species(scientific_name: str) -> str:
    genus, species = scientific_name.split(" ")
    return f"{genus[0]}. {species}"

def plot_jaspar_prop_match(ax):
    clade_list = list(clades.keys())
    for config in ['accessibility_base', 'hypomethylation_base']:
        data = all_high_effect_jaspar_prop.xs(config, level = 'config')
        species = data.index.get_level_values('species').unique()
        x_pos = np.array([clade_list.index(s) for s in species])
        width = 0.35

        if config == "accessibility_base":
            ax.bar(x_pos - (width / 2), data['JASPAR_proportion'], width = width, label = "Accessibility")
        elif config == "hypomethylation_base":
            ax.bar(x_pos + (width / 2), data['JASPAR_proportion'], width = width, label = "Methylation")
        else:
            raise ValueError(f"Unknown config name '{config}'")
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(len(clades.keys()))))
    ax.set_xlabel("Species")
    ax.set_xticklabels([shorten_species(s) for s in clades.keys()], rotation = 45, ha = 'right')
    ax.set_ylabel("JASPAR-matched proportion")
    ax.legend()
```

```{code-cell} ipython3
all_high_effect_jaspar_prop = pd.read_table(
    "tmp/all_high_effect_jaspar_prop.tsv",
    index_col = [0, 1]
)
```

```{code-cell} ipython3
fig3, axs = plt.subplots(2, 2, figsize = (12, 8), constrained_layout = True, facecolor = "white")

plot_pGIA("accessibility_base", "Arabidopsis thaliana", axs[0][0], legend_motif_pos = "upper left", legend_strand_pos = "upper right")
axs[0][0].set_title("Arabidopsis thaliana")
axs[0][0].set_ylabel("Accessibility")
plot_pGIA("accessibility_base", "Zea mays", axs[0][1])
axs[0][1].set_title("Zea mays")
axs[0][1].set_ylabel(None)
plot_pGIA("hypomethylation_base", "Arabidopsis thaliana", axs[1][0])
axs[1][0].set_ylabel("Methylation")
plot_pGIA("hypomethylation_base", "Zea mays", axs[1][1])
axs[1][1].set_ylabel(None)
fig3.supylabel("Global Importance")

fig3.savefig("figs/fig3.png", dpi = 300)
```

```{code-cell} ipython3
fig4, axs = plt.subplots(nrows = 1, ncols = 3, figsize = (14, 4), constrained_layout = True, facecolor = "white")

plot_embedding_by_config(axs[0])
plot_embedding_by_clade(axs[1], "accessibility_base")
axs[1].set_title("Accessibility")
plot_embedding_by_clade(axs[2], "hypomethylation_base")
axs[2].set_title("Methylation")

fig4.savefig("figs/fig4.png", dpi = 300)
```

```{code-cell} ipython3
fig, axs = plt.subplots(constrained_layout = True, nrows = 1, ncols = 2, figsize = (9, 4), facecolor = "white")

plot_embedding_by_species(axs[0], "accessibility_base")
axs[0].set_title("Accessibility")
plot_embedding_by_species(axs[1], "hypomethylation_base")
axs[1].set_title("Methylation")

fig.savefig("figs/mds-species.png", dpi = 300)
```

```{code-cell} ipython3
fig, ax = plt.subplots(constrained_layout = True, facecolor = "white", figsize = (6, 5))
plot_jaspar_prop_match(ax)

fig.savefig("figs/jaspar-match.png", dpi = 300)
```

```{code-cell} ipython3

```
