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
import functools
import hashlib
import itertools
from pathlib import Path
import random

import Bio.Align
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.manifold

from a2z.utils import kMedoids
```

```{code-cell} ipython3
def shorten_species(scientific_name: str) -> str:
    genus, species = scientific_name.split(" ")
    return f"{genus[0]}. {species}"

def get_config_species_medoids(config_name: str, species: str, kmer_cutoff: int = 1000, medoids: int = 100, kmer_size: int = 10) -> list[str]:
    high_effect_kmers = pd.read_table(
        Path(f"tmp/kmer-occlusion/{config_name}/{species.replace(' ', '_')}/{kmer_size}/high_effect_kmers.tsv"),
        index_col = 0,
        dtype = {
            'masked_seq': str,
            'delta_prediction': float
        }
    ).rename(columns = {'delta_pred': 'effect_size'}).iloc[:kmer_cutoff]
    kmers = high_effect_kmers.index.to_list()
    
    aligner = Bio.Align.PairwiseAligner()
    
    # perform all-by-all alignment and get similarity scores
    proximity_matrix = np.diag(np.full(shape = (len(kmers)), fill_value = len(kmers[0]), dtype = np.byte))

    for i, idx in enumerate(itertools.combinations(range(len(kmers)), 2)):
        proximity_matrix[idx[0], idx[1]] = aligner.score(kmers[idx[0]], kmers[idx[1]])
        
    # mirror across the diagonal
    proximity_matrix = proximity_matrix + proximity_matrix.T - np.diag(np.diag(proximity_matrix))
    # convert to distance matrix
    distance_matrix = proximity_matrix.max() - proximity_matrix
    
    medoid_idx, clusters = kMedoids(distance_matrix, k = medoids)
    
    return [kmers[i] for i in medoid_idx]
```

```{code-cell} ipython3
config_list = [
    "accessibility_base",
    "hypomethylation_base"
]
```

```{code-cell} ipython3
%%time
medoids = None
for config_name in config_list:
    species_list = [child.stem.replace('_', ' ') for child in Path(f"tmp/kmer-occlusion/{config_name}").iterdir()]
    # subset to top 10 medoids ranked by effect size
    df = pd.DataFrame({species: get_config_species_medoids(config_name, species)[:10] for species in species_list}).melt(var_name = "species", value_name = "kmer")
    df['config'] = config_name
    if medoids is None:
        medoids = df
    else:
        medoids = medoids.append(df)
medoids = medoids.reset_index(drop = True)
```

```{code-cell} ipython3
medoids
```

```{code-cell} ipython3
aligner = Bio.Align.PairwiseAligner()
# perform all-by-all alignment and get similarity scores
proximity_matrix = np.diag(np.full(shape = (len(medoids)), fill_value = len(medoids.iloc[0]), dtype = np.byte))

kmers = medoids['kmer']
for i, idx in enumerate(itertools.combinations(range(len(kmers)), 2)):
    proximity_matrix[idx[0], idx[1]] = aligner.score(kmers[idx[0]], kmers[idx[1]])
    
# mirror across the diagonal
proximity_matrix = proximity_matrix + proximity_matrix.T - np.diag(np.diag(proximity_matrix))
# convert to distance matrix
distance_matrix = proximity_matrix.max() - proximity_matrix
```

```{code-cell} ipython3
distances_flattened = np.triu(distance_matrix).flatten()
plt.hist(distances_flattened[distances_flattened > 0])
```

```{code-cell} ipython3
%%time
embedding = sklearn.manifold.MDS(dissimilarity = 'precomputed', random_state = 42)
embedded_kmers = embedding.fit_transform(distance_matrix)
```

```{code-cell} ipython3
medoids_embedded = medoids.join(pd.DataFrame(embedded_kmers, columns = ['embedded1', 'embedded2']))
```

```{code-cell} ipython3
medoids_embedded.to_csv(
    "tmp/medoids_embedded.tsv",
    sep = "\t",
    index = False
)
```

```{code-cell} ipython3
fig, ax = plt.subplots(tight_layout = True, facecolor = "white")

labels = {
    'accessibility_base': 'Accessibility',
    'hypomethylation_base': 'Methylation'
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
fig.savefig("tmp/mds_accessibility_methylation.png", dpi = 300)
```

```{code-cell} ipython3
fig, ax = plt.subplots(constrained_layout = True, facecolor = "white")

for species in medoids.loc[medoids['config'] == "accessibility_base"]['species'].unique():
    species_rows = medoids.loc[(medoids['config'] == "accessibility_base") & (medoids['species'] == species)].index.to_numpy()[:, None]
    ax.scatter(embedded_kmers[species_rows, 0], embedded_kmers[species_rows, 1], alpha = 0.5, label = shorten_species(species))
x_lim = ax.get_xlim()
#ax.set_xlim(x_lim[0], (x_lim[1] - x_lim[0]) * 1.)
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_xlabel("Dimension 1")
y_lim = ax.get_ylim()
#ax.set_ylim(y_lim[0], (y_lim[1] - y_lim[0]) * 1.)
ax.set_yticks([])
ax.set_yticklabels([])
ax.set_ylabel("Dimension 2")
ax.legend()
fig.savefig("tmp/mds_accessibility_species.png", dpi = 300)
```

```{code-cell} ipython3
fig, ax = plt.subplots(constrained_layout = True, facecolor = "white")

for species in medoids.loc[medoids['config'] == "hypomethylation_base"]['species'].unique():
    species_rows = medoids.loc[(medoids['config'] == "hypomethylation_base") & (medoids['species'] == species)].index.to_numpy()[:, None]
    ax.scatter(embedded_kmers[species_rows, 0], embedded_kmers[species_rows, 1], alpha = 0.5, label = shorten_species(species))
x_lim = ax.get_xlim()
#ax.set_xlim(x_lim[0], (x_lim[1] - x_lim[0]) * 1.)
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_xlabel("Dimension 1")
y_lim = ax.get_ylim()
#ax.set_ylim(y_lim[0], (y_lim[1] - y_lim[0]) * 1.)
ax.set_yticks([])
ax.set_yticklabels([])
ax.set_ylabel("Dimension 2")
ax.legend()
fig.savefig("tmp/mds_methylation_species.png", dpi = 300)
```

```{code-cell} ipython3

```
