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
pip install --editable src/python/a2z
```

After setting up the environment, you can then use the scripts in `data/` to download all the necessary data for the experiments in the `results/` folder.

### Optional Developer Setup

```
git config --local --add diff.jupyternotebook.command "git-nbdiffdriver diff"
git config --local --add merge.jupyternotebook.command "git-nbmergedriver merge %O %A %B %L %P"
echo "*.ipynb diff=jupyternotebook merge=jupyternotebook" >> .git/info/attributes
```

