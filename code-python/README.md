## Tools used for the analysis of neural recording data

## Getting started

* Download Miniconda from [Continuum Analytics](http://conda.pydata.org/miniconda.html). vdmlab\code-python uses Python 3.5 (verified working on Windows and Linux x64).
* In the terminal, create and activate a new conda environment.

```
conda create -n yourenv python=3.5
activate yourenv [Windows] or source activate yourenv [Linux]
```

* Install package dependencies (it's possible to install multiple packages at once or individually). If conda doesn't have a package of interest (eg. shapely), on in the terminal try: `pip install shapely`. For Windows, download the most recent *.whl file [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely) and install using `pip install Shapely-1.5.16-cp35-cp35m-win_amd64.whl`

```
conda install numpy scipy shapely
```

* Clone the analysis code from Github and developer installation.

```
git clone https://github.com/mvdm/vandermeerlab.git
cd vandermeerlab/code-python
python setup.py develop
```

* **All set!**


## Documentation

Build documentation (eg. from docstrings) using `python setup.py build_sphinx` in the vdmlab directory


