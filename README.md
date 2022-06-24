# Cross-species chromatin state prediction

## Dependencies

- [conda](https://conda.io)
- [git-lfs](https://git-lfs.github.com)

The authors recommend [mambaforge](https://mamba.readthedocs.io/en/latest/installation.html) as a faster alternative to `conda`.
If using mambaforge, use the `mamba create` command instead of `conda create`.

## Environment Setup

```
git clone https://bitbucket.org/bucklerlab/a2z-regulatory
cd a2z-regulatory
conda env create -f envs/a2z.yml
conda activate a2z
pip install src/python/a2z
```

After setting up the environment, you can then use the scripts in `data/` to download all the necessary data for the experiments in the `results/` folder.

### Developer Setup (optional)

This is only useful for better `git-diff` output on Jupyter Notebooks.

```
git config --local --add diff.jupyternotebook.command "git-nbdiffdriver diff"
git config --local --add merge.jupyternotebook.command "git-nbmergedriver merge %O %A %B %L %P"
echo "*.ipynb diff=jupyternotebook merge=jupyternotebook" >> .git/info/attributes
```

## Predicting on your own sequences

### Running predictions using a2z

Use `predict_genome.py -h` to see all available options.

```
(base)$ curl 'https://zenodo.org/record/5724562/files/model-accessibility-full.h5?download=1' > model-accessibility-full.h5
(base)$ conda activate a2z
(a2z)$ src/python/scripts/predict_genome.py model-accessibility-full.h5 genome.fa > preds.bed
```

Now you will have a BED3+1 file with predictions in 600bp windows sliding 50bp at a time.
You can do simple "peak-calling" with `bedtools`.

```
(a2z)$ samtools faidx genome.fa
(a2z)$ bedtools map -a <(bedtools makewindows -g genome.fa.fai -w 50) -b preds.bed -c 4 -o max | awk '$4 >= 0.9' | bedtools merge > preds.gte90.merged.bed
```

The above will take the max prediction of all overlapping windows every 50bp and then merge touching 50bp chunks with predictions greater than or equal to 0.90 into whole regions.
You can also mask out repeats with bedtools and an annotation.

```
(a2z)$ grep -v '^#' repeats.gff3 | awk -v 'OFS=\t' '{print $1, $4 - 1, $5}' > repeats.bed
(a2z)$ bedtools subtract -A -f 0.5 -a preds.gte90.merged.bed -b repeats.bed > preds.gte90.merged.masked.bed
```

### Running predictions using Kipoi

The a2z leaf chromatin and DNA methylation models are also available in the [Kipoi](https://kipoi.org) model zoo.

The following will make leaf chromatin accessibility predictions on 600bp sliding windows with a stride of 50bp.
You can adjust the window stride by modifying the `-s` parameter to `bedtools`.

```
(base)$ ls
a2z-regulatory
(base)$ mkdir work
(base)$ cd work
(base)$ curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz | gzip -cd > tair.fa
(base)$ conda activate a2z
(a2z)$ samtools faidx tair.fa
(a2z)$ bedtools makewindows -g tair.fa.fai -w 600 -s 50 | awk '$3-$2 == 600' > windows.bed
(a2z)$ conda env create -f ../a2z-regulatory/envs/kipoi.yml
(a2z)$ conda activate a2z-kipoi
(a2z-kipoi)$ kipoi env export a2z-chromatin/a2z-accessibility -o env.yml
(a2z-kipoi)$ conda env create -f env.yml
(a2z-kipoi)$ conda activate kipoi-a2z-chromatin__a2z-accessibility
(kipoi-a2z-chromatin__a2z-accessibility)$ kipoi predict a2z-chromatin/a2z-accessibility -n 8 --batch_size=128 --dataloader_args='{"intervals_file": "windows.bed", "fasta_file": "tair.fa", "num_chr_fasta": true}' -o a2z.acc.preds.tsv
```

