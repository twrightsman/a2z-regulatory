# Manuscript

Based on the [Manubot Rootstock](https://github.com/manubot/rootstock).


## Setup

```
$ conda env create -f build/environment.yml
$ conda activate manubot
(manubot) $
```


## Editing

Simply edit files under `content/` and render (see below) when you want to see changes.


## Rendering

```
(manubot) $ ./build/build.sh
```

The rendered manuscript will be in the `output/` directory in HTML format.
To view it, open up `manuscript.html` in a browser.

If you prefer PDF or DOCX files, you can build them by defining the appropriate environment variables while running the build script, as shown below.

```
(manubot) $ BUILD_PDF=true BUILD_DOCX=true ./build/build.sh
```


## Serve generated webpage

```
python -m http.server 8080 --bind 127.0.0.1 --directory webpage
```

