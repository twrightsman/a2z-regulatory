## Code

The code used for all analyses is in `a2z-regulatory.tgz`.

```
$ mkdir a2z-regulatory
$ tar -C a2z-regulatory -xvf a2z-regulatory.tgz
```

## UMR calls

The UMR calls are the only data needed to reproduce the entire manuscript from
scratch that aren't publicly-available from other sources.

```
$ mkdir -p a2z-regulatory/data/crisp_2021_umr_angiosperms/raw
$ tar -C a2z-regulatory/data/crisp_2021_umr_angiosperms/raw -xvf umr-calls.tgz
```

## Pre-trained models

`models.tgz` contains all of the models trained for the manuscript.

```
$ mkdir -p a2z-regulatory/results/01_models/tmp/results
$ tar -C a2z-regulatory/results/01_models/tmp/results -xvf models.tgz
```

`model-accessibility-full.h5` and `model-methylation-full.h5` contain the
models trained on all species used in the manuscript.

## JASPAR pGIA results

`JASPAR-pGIA.tgz` contains the pGIA results for all 530 JASPAR 2020 TFs.

```
$ mkdir -p a2z-regulatory/results/03_interpretation/tmp/kmer-occlusion
$ tar -C a2z-regulatory/results/03_interpretation/tmp/kmer-occlusion -xvf JASPAR-pGIA.tgz
```

