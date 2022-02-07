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

```{code-cell} ipython3
from glob import glob
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
```

```{code-cell} ipython3
model_metrics = {
    'model': [],
    'test_species': [],
    'replicate': [],
    'auPR': [],
    'auROC': []
}
for metrics_path in map(lambda x: Path(x), glob("tmp/results/accessibility_*/*/?/metrics.json")):
    _, _, model_name, species, replicate, _ = metrics_path.parts
    if model_name.split('_', maxsplit=1)[1] not in {'base', 'Basset', 'CharPlant', 'DeeperDeepSEA'}:
        continue
    with open(metrics_path) as metrics_file:
        metrics = json.load(metrics_file)

    model_name = model_name.split("_")[1]
    model_name = "DanQ" if model_name == "base" else model_name

    model_metrics['model'].append(model_name)
    model_metrics['test_species'].append(species)
    model_metrics['replicate'].append(replicate)
    model_metrics['auPR'].append(metrics['auPR'])
    model_metrics['auROC'].append(metrics['auROC'])

model_metrics = pd.DataFrame(model_metrics)
```

```{code-cell} ipython3
model_metrics_agg = model_metrics.groupby(['test_species', 'model']).mean()
```

```{code-cell} ipython3
# order the model architectures along the x-axis by mean auPR across test species
model_plot_order = model_metrics_agg.reset_index().groupby('model').mean().sort_values('auPR').index.to_list()
```

```{code-cell} ipython3
model_metrics_agg_sorted = model_metrics_agg.sort_index(level = 'model', key = lambda i: pd.Index([model_plot_order.index(m) for m in i], name = i.name))
```

```{code-cell} ipython3
fig, ax = plt.subplots(constrained_layout = True, facecolor = 'white')

for i, test_species in enumerate(sorted(model_metrics['test_species'].unique())):
    auPRs = model_metrics_agg_sorted.xs(key = test_species, level = 'test_species')['auPR']
    x = [model_plot_order.index(i) for i in auPRs.index]
    label = ''.join((w[0] for w in test_species.split("_")))
    ax.plot(x, auPRs, label = label, color = plt.cm.tab20(i), zorder = 1)
    ax.scatter(x, auPRs, color = plt.cm.tab20(i), zorder = 2)

ax.set_ylim((0, 0.7))
ax.set_xticks(range(len(model_plot_order)))
ax.set_xticklabels(model_plot_order)
ax.set_xlabel('Model')
ax.set_ylabel('auPR')
ax.legend(ncol = 3, loc = 'upper left')

fig.savefig('figs/model_cmp.png', dpi = 300)
```

```{code-cell} ipython3

```
