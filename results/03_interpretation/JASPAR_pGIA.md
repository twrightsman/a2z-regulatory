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
from glob import glob
from itertools import repeat
import json
from pathlib import Path
from statistics import mean
import urllib.request

import IPython
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from a2z.data.motifs import MinimalMEME
```

```{code-cell} ipython3
results = None
for path in map(lambda s: Path(s), glob("tmp/kmer-occlusion/*/*/JASPAR_gia.tsv")):
    df = {
        'motif': [],
        'position': [],
        'strand': [],
        'global_importance': []
    }
    for line in open(path):
        motif, *pGIA = line.rstrip().split("\t")
        assert (len(pGIA) % 2) == 0
        n_pos = len(pGIA) // 2
        for i, global_importance in enumerate(pGIA):
            df['motif'].append(motif)
            df['position'].append(i % n_pos)
            if i < n_pos:
                df['strand'].append('+')
            else:
                df['strand'].append('-')
            df['global_importance'].append(global_importance)
    model, species = path.parts[2:4]
    df = pd.DataFrame(data = df)
    df['model'] = model
    df['species'] = species.replace("_", " ")
    df = df.set_index(['model', 'species', 'motif', 'strand', 'position'])
    df['global_importance'] = df['global_importance'].astype(np.float)
    if results is None:
        results = df
    else:
        results = results.append(df)
```

```{code-cell} ipython3
results
```

```{code-cell} ipython3
# make sure each strand has same number of positions
for (e, s, m), subdata in results.groupby(['model', 'species', 'motif']):
    assert len(subdata.xs('+', level = 'strand')) == len(subdata.xs('-', level = 'strand'))
```

```{code-cell} ipython3
results.to_csv(
    "tmp/JASPAR_pGIA.all.melt.tsv",
    sep = "\t"
)
```

```{code-cell} ipython3
top_jaspar_max = results.groupby(
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
top_jaspar_max
```

```{code-cell} ipython3
JASPAR_name2id = {motif.name: motif.identifier for key, motif in MinimalMEME('../../data/JASPAR/JASPAR2020_CORE_plants_non-redundant_pfms_meme.txt').motifs.items()}

def shorten_species(scientific_name: str) -> str:
    genus, species = scientific_name.split(" ")
    return f"{genus[0]}. {species}"

def get_JASPAR_matrix_info(matrix_id: str) -> dict:
    matrix_cache_path = Path(f"tmp/jaspar-cache/matrix/{matrix_id}.json")
    if not matrix_cache_path.exists():
        matrix_cache_path.parent.mkdir(parents = True, exist_ok = True)
        with urllib.request.urlopen(f"http://jaspar.genereg.net/api/v1/matrix/{matrix_id}") as response:
            if (response.status == 200) and (response.headers['Content-Type'] == 'application/json'):
                matrix_data = json.loads(response.read())
                with open(matrix_cache_path, 'w') as matrix_cache_file:
                    json.dump(matrix_data, matrix_cache_file, indent = 2)
            else:
                raise RuntimeError(f"Failed to get data from JASPAR matrix API for '{matrix_id}'")
    else:
        with open(matrix_cache_path) as matrix_cache_file:
            matrix_data = json.load(matrix_cache_file)

    return matrix_data

def JASPAR_name_to_family_or_class(JASPAR_TF_name: str) -> str:
    matrix_info = get_JASPAR_matrix_info(JASPAR_name2id[JASPAR_TF_name])
    tf_families = matrix_info['family']
    tf_classes = matrix_info['class']
    if len(tf_families) > 0:
        return tf_families[0]
    if len(tf_classes) > 0:
        return tf_classes[0]
    return 'Unknown'
```

```{code-cell} ipython3
def format_cell(tf: str) -> str:
    family_or_class = JASPAR_name_to_family_or_class(tf)
    if family_or_class.endswith(' domain'):
        family_or_class = family_or_class[:-7]
    return f"{tf}<br>({family_or_class})"

JASPAR_ranked_table = results.groupby(
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
).head(10).reset_index().drop(
    columns = 'global_importance'
)
JASPAR_ranked_table['Rank'] = (JASPAR_ranked_table.index % 10) + 1
JASPAR_ranked_table.loc[JASPAR_ranked_table['model'] == 'hypomethylation_base', 'model'] = 'Methylation'
JASPAR_ranked_table.loc[JASPAR_ranked_table['model'] == 'accessibility_base', 'model'] = 'Accessibility'
JASPAR_ranked_table['species'] = list(map(shorten_species, JASPAR_ranked_table['species']))
JASPAR_ranked_table = JASPAR_ranked_table.sort_values(['model', 'species'])
JASPAR_ranked_table = JASPAR_ranked_table.pivot(
    index = 'Rank',
    columns = ['model', 'species'],
    values = 'motif'
)
JASPAR_ranked_table.columns = JASPAR_ranked_table.columns.rename({'model': 'Feature', 'species': 'Species'})
JASPAR_ranked_table = JASPAR_ranked_table.applymap(format_cell)
JASPAR_ranked_table.columns = [f"{model}<br>{species}" for model, species in JASPAR_ranked_table.columns.to_flat_index()]
JASPAR_ranked_table.to_markdown(
    "tables/JASPAR_ranked_table.md",
    stralign = "center"
)
```

```{code-cell} ipython3

```
