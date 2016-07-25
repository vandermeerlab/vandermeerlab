## Tools used for the analysis of neural recording data

## Getting started

* Download Miniconda from
  [Continuum Analytics](http://conda.pydata.org/miniconda.html).
  We recommend the Python 3 version.
* Open a *new* terminal, create and activate a new conda environment.

  ```
  conda create -n yourenv python=3.5
  activate yourenv [Windows] or source activate yourenv [Linux]
  ```

* Install package dependencies (it's possible to
  install multiple packages at once or individually).
  If conda doesn't have a package of interest (eg. shapely),
  in the terminal try: `pip install shapely`.
  In Windows, download the most recent `*.whl` file
  [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely)
  and install using `pip install yourshapelyinstall.whl`
  (must be in the directory where this .whl is located).

  ```
  conda install numpy scipy shapely pytest matplotlib
  ```

* Clone the analysis code from Github and developer installation.

  ```
  git clone https://github.com/mvdm/vandermeerlab.git
  cd vandermeerlab/code-python
  python setup.py develop
  ```

* **All set!**

## Documentation

Build documentation (eg. from docstrings) using 
`python setup.py build_sphinx` in the vdmlab directory.

## Testing

Run tests with [pytest](http://docs.pytest.org/en/latest/usage.html).
