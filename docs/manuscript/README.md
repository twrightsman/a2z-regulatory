# Manuscript

Based on the [Manubot Rootstock](https://github.com/manubot/rootstock).


## Setup

```
$ conda env create -f build/environment.yml
$ conda activate manubot
(manubot) $ wget --directory-prefix=build/pandoc/filters https://github.com/pandoc/lua-filters/raw/13c3fa7e97206413609a48a82575cb43137e037f/spellcheck/spellcheck.lua
```


## Editing

Simply edit files under `content/` and render (see below) when you want to see changes.


## Rendering

```
(manubot) $ SPELLCHECK=true ./build/build.sh
```

The rendered manuscript will be in the `output/` directory in HTML format.
To view it, open up `manuscript.html` in a browser.

If you prefer PDF or DOCX files, you can build them by defining the appropriate environment variables while running the build script, as shown below.

```
(manubot) $ SPELLCHECK=true BUILD_PDF=true BUILD_DOCX=true ./build/build.sh
```


## Serve generated webpage

```
python -m http.server 8080 --bind 127.0.0.1 --directory webpage
```


## Generate differences between manuscript versions

### Install pandiff using conda

```
conda install -c conda-forge nodejs
npm install -g pandiff
```

### Generate manuscripts to compare

```
git checkout $OLD_REVISION
BUILD_PDF=true BUILD_DOCX=true ./build/build.sh
cp output/manuscript.md output/manuscript.old.md
git checkout $NEW_REVISION
BUILD_PDF=true BUILD_DOCX=true ./build/build.sh
pandiff --output=output/pandiff.html --reference-links --standalone --resource-path=content --wrap=preserve output/manuscript.old.md output/manuscript.md
```
