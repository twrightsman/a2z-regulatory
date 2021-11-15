---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

```{code-cell} ipython3
import hashlib
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pybedtools import BedTool, Interval
from pyfaidx import Fasta
from sklearn.metrics import average_precision_score, precision_recall_curve
import tensorflow as tf
import tensorflow.keras as k

from a2z.data import IntervalReferenceDataset
from a2z.models import allow_growth

allow_growth()
```

```{code-cell} ipython3
def plot_masked_model_performance(config_name: str, ax: matplotlib.axes.Axes):
    with open(f"configs/{config_name}.json") as config_file:
        config = json.load(config_file)
        
    no_maize_model_path = (f"tmp/results/{config_name}/Zea_mays/0/model")
    no_maize_model = k.models.load_model(no_maize_model_path)
    
    split_hash = hashlib.sha256((json.dumps(config['paths'], indent = 2) + "\n" + json.dumps(config['preprocessing'], indent = 2) + "\n" + json.dumps(config['splitting'], indent = 2) + "\n").encode('utf-8')).hexdigest()
    
    intervals = pd.read_table(
        f"tmp/data/split/{split_hash}.tsv",
        header = 0,
        dtype = {
            'seqid': str,
            'start': int,
            'end': int,
            'distance_class': 'category',
            'species': str,
            'target': np.byte,
            'is_train': np.byte,
            'is_test': np.byte,
            'is_val': np.byte
        }
    )
    test_intervals = intervals.loc[(intervals['species'] == 'Zea mays') & (intervals['is_test'] == 1)].copy()
    references = {'Zea mays': Fasta(str(Path(config['paths']['genomes']) / f"Zea_mays.fa"), as_raw = True, rebuild = False)}
    test_data = IntervalReferenceDataset(test_intervals, references)
    test_intervals['prediction'] = no_maize_model.predict(test_data).flatten()
    
    precision, recall, thresholds = precision_recall_curve(test_intervals['target'], test_intervals['prediction'])
    auPR = average_precision_score(test_intervals['target'], test_intervals['prediction'])
    ax.plot(recall, precision, label = f"Standard (auPR = {auPR:.2f})")
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 1))
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    
    repeats = BedTool("../../data/maize_RM/v4.repeats.renamed.bed")
    test_regions = BedTool((Interval(row.seqid, row.start, row.end) for row in test_intervals.itertuples()))
    repeat_regions = test_regions.intersect(repeats, wa = True, u = True, f = 0.5)
    
    repeat_df = repeat_regions.to_dataframe().rename(columns = {'chrom': 'seqid'}).astype({'seqid': str}).set_index(['seqid', 'start', 'end'])
    test_df = test_intervals.reset_index().set_index(['seqid', 'start', 'end'])
    joined = test_df.join(repeat_df, how = 'inner')
    repeat_overlapping_idx = joined['index'].values
    
    masked_test_intervals = test_intervals.copy()
    masked_test_intervals.loc[repeat_overlapping_idx, 'prediction'] = 0
    
    masked_precision, masked_recall, masked_thresholds = precision_recall_curve(masked_test_intervals['target'], masked_test_intervals['prediction'])
    masked_auPR = average_precision_score(masked_test_intervals['target'], masked_test_intervals['prediction'])
    ax.plot(masked_recall, masked_precision, label = f"Masked (auPR = {masked_auPR:.2f})")
    ax.legend()
```

```{code-cell} ipython3
%%time

fig, axs = plt.subplots(constrained_layout = True, nrows = 1, ncols = 2, figsize = (12, 5))
plot_masked_model_performance('accessibility_base', axs[0])
axs[0].set_title("Accessibility")
plot_masked_model_performance('hypomethylation_base', axs[1])
axs[1].set_title("Methylation")
```

```{code-cell} ipython3
fig.savefig('figs/maize-masked.png', dpi = 300, facecolor = "white")
```

```{code-cell} ipython3

```
