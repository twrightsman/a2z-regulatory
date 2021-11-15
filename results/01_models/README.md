# Model Training and Evaluation

Do the following for each config in `configs/`.
`accessibility_unbalanced.json` can be skipped as it takes a very long time and isn't included in the manuscript results.
If you want to use multiple GPUs in parallel you can modify the `simple_gpu_scheduler` command line options.

```sh
./01_preprocess.sh configs/accessibility_base.json
./02_split.py configs/accessibility_base.json
./03_generate_jobs.py configs/accessibility_base.json > tmp/jobs
simple_gpu_scheduler --gpus 0 < tmp/jobs > tmp/jobs.out 2> tmp/jobs.err
```

- `fig1.ipynb` requires the `accessibility_base`, `accessibility_base_within`, `hypomethylation_base`, and `hypomethylation_base_within` configs.

For running a few custom models for the comparison supplemental figure:

```sh
simple_gpu_scheduler --gpus 0 < jobs_custom > tmp/jobs.out 2> tmp/jobs.err
```

## Preprocess yeast and human data

```
./preprocess_single.sh 600 ../../data/schep_2015_yeast_ATAC/yeast.open.bed ../../data/schep_2015_yeast_ATAC/yeast.renamed.fa.fai
./preprocess_single.sh 600 ../../data/schep_2015_yeast_ATAC/GM.peaks.final.bed ../../data/schep_2015_yeast_ATAC/hg19.fa.fai
```

## Get smaller dataset deterministically

Credit: https://stackoverflow.com/questions/60266215/shuffle-output-of-find-with-fixed-seed

```
get_fixed_random() { openssl enc -aes-256-ctr -pass pass:"$1" -nosalt </dev/zero 2>/dev/null; }
tail -n+2 tmp/data.split.tsv | shuf --random-source=<(get_fixed_random 42) -n 100000 | cat <(head -n1 tmp/data.split.tsv) - > tmp/data.small.split.tsv
```

