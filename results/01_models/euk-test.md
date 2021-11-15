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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pybedtools import BedTool, Interval
from pyfaidx import Fasta
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve
import tensorflow as tf
import tensorflow.keras as k

from a2z.data import IntervalReferenceDataset
from a2z.models import allow_growth

allow_growth()
```

```{code-cell} ipython3
config_name = "accessibility_base"
with open(f"../01_models/configs/{config_name}.json") as config_file:
    config = json.load(config_file)
```

```{code-cell} ipython3
references = {
    'Saccharomyces cerevisiae': Fasta("../../data/schep_2015_yeast_ATAC/yeast.renamed.fa", as_raw = True, rebuild = False),
    'Homo sapiens': Fasta("../../data/schep_2015_yeast_ATAC/hg19.fa", as_raw = True, rebuild = False)
}
```

```{code-cell} ipython3
# load in pan-Angiosperm model trained on everything but maize
# can take a random model (here it's just the first) since the auPR doesn't change much between different samples of the same model
no_maize_model_path = (f"../01_models/tmp/results/{config_name}/Zea_mays/0/model")
no_maize_model = k.models.load_model(no_maize_model_path)
```

```{code-cell} ipython3
yeast_intervals = pd.read_table(
    "tmp/yeast.unlabeled.bed",
    names = ['seqid', 'start', 'end', 'target'],
    dtype = {
        'seqid': str,
        'start': int,
        'end': int,
        'target': np.byte
    }
)
yeast_intervals['species'] = "Saccharomyces cerevisiae"
```

```{code-cell} ipython3
yeast_dataset = IntervalReferenceDataset(yeast_intervals, references)
```

```{code-cell} ipython3
%%time
yeast_intervals['prediction'] = no_maize_model.predict(yeast_dataset).flatten()
```

```{code-cell} ipython3
human_intervals = pd.read_table(
    "tmp/GM.unlabeled.bed",
    names = ['seqid', 'start', 'end', 'target'],
    dtype = {
        'seqid': str,
        'start': int,
        'end': int,
        'target': np.byte
    }
)
human_intervals['species'] = "Homo sapiens"
```

```{code-cell} ipython3
human_dataset = IntervalReferenceDataset(human_intervals, references)
```

```{code-cell} ipython3
%%time

human_cache = Path("tmp/GM.preds.tsv")
if human_cache.exists():
    preds = pd.read_table(
        human_cache,
        names = ['prediction'],
        dtype = {'prediction': float}
    )
    human_intervals['prediction'] = preds['prediction']
else:
    # ~30 minutes on RTX 2080
    human_intervals['prediction'] = no_maize_model.predict(human_dataset).flatten()
    human_intervals['prediction'].to_csv(
        human_cache,
        sep = "\t",
        header = False,
        index = False
    )
```

```{code-cell} ipython3
fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (11, 5), constrained_layout = True, facecolor = "white")

yeast_precision, yeast_recall, yeast_thresholds = precision_recall_curve(yeast_intervals['target'], yeast_intervals['prediction'])
yeast_auPR = average_precision_score(yeast_intervals['target'], yeast_intervals['prediction'])
human_precision, human_recall, human_thresholds = precision_recall_curve(human_intervals['target'], human_intervals['prediction'])
human_auPR = average_precision_score(human_intervals['target'], human_intervals['prediction'])

axs[0].set_xlim((0, 1))
axs[0].set_ylim((0, 1))
axs[0].set_xlabel("Recall")
axs[0].set_ylabel("Precision")

axs[0].plot(yeast_recall, yeast_precision, label = f"Yeast (auPR = {yeast_auPR:.2f})", color = "orange")
axs[0].plot(human_recall, human_precision, label = f"Human Cell Line (auPR = {human_auPR:.2f})", color = "brown")

axs[0].legend()


yeast_fpr, yeast_tpr, yeast_thresholds = roc_curve(yeast_intervals['target'], yeast_intervals['prediction'])
yeast_auROC = roc_auc_score(yeast_intervals['target'], yeast_intervals['prediction'])
human_fpr, human_tpr, human_thresholds = roc_curve(human_intervals['target'], human_intervals['prediction'])
human_auROC = roc_auc_score(human_intervals['target'], human_intervals['prediction'])

axs[1].plot([0, 1], [0, 1], color="gray", linestyle="--", label="Baseline")
axs[1].set_xlim((0, 1))
axs[1].set_ylim((0, 1))
axs[1].set_xlabel("False Positive Rate")
axs[1].set_ylabel("True Positive Rate")

axs[1].plot(yeast_fpr, yeast_tpr, label = f"Yeast (auROC = {yeast_auROC:.2f})", color = "orange")
axs[1].plot(human_fpr, human_tpr, label = f"Human Cell Line (auROC = {human_auROC:.2f})", color = "brown")

axs[1].legend()

fig.savefig("figs/euk-cmp.png", dpi = 300)
```

```{code-cell} ipython3

```
