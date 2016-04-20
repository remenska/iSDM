# iSDM
Towards a species-by-species approach to global biodiversity modelling.


##### Step 1:  Download [Anaconda](https://www.continuum.io/downloads) for Python 3.

NOTE for Windows:
Do not install as Administrator unless admin privileges are required. If you encounter any issues during installation, please temporarily disable your anti-virus software during install, then immediately re-enable it.

##### Step 2: Create a python environment and install the following packages:
```
conda create --name=biodiversity six pandas ipython-notebook scikit-learn git basemap matplotlib xlrd numba python=3
```
Linux: &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```source activate biodiversity```

Windows:&nbsp;&nbsp; ```activate biodiversity```

##### Step 3: Install these additional packages:
```
pip install pygbif rpy2 geopy geopandas
```

##### Step 4: Install iSDM: 
```
git clone https://github.com/remenska/iSDM.git
cd iSDM/
python setup.py install
```

##### Step 5: Test if it works:
```
$ ipython
Python 3.5.1 |Continuum Analytics, Inc.| (default, Dec  7 2015, 11:16:01) 
[...]
In [1]: import iSDM

In [2]: iSDM.__version__
Out[2]: '0.0.1'
```


