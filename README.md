# iSDM [![Build Status](https://travis-ci.org/remenska/iSDM.svg?branch=master)](https://travis-ci.org/remenska/iSDM) [![Join the chat at https://gitter.im/remenska/iSDM](https://badges.gitter.im/remenska/iSDM.svg)](https://gitter.im/remenska/iSDM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
Towards a species-by-species approach to global biodiversity modelling.


##### Step 1:  Download [Anaconda](https://www.continuum.io/downloads) for Python 3.

NOTE for Windows:
Do not install as Administrator unless admin privileges are required. If you encounter any issues during installation, please temporarily disable your anti-virus software during install, then immediately re-enable it.
In the dialog **Advanced Options**, check both options: *Add Anaconda to my PATH environment variable*, and *Register Anaconda as my default Python*. 

NOTE for Linux (and Mac OSX): to add Anaconda to your PATH (environment variable) permanently, it's best to say 'yes' to appending to .bashrc, during installation. Otherwise, to temporarily make Anaconda visible in your environment, do:

```export PATH=~/anaconda3/bin:$PATH  # assuming you installed Anaconda at ~/anaconda3/```
                                   
##### Step 2: Create a python environment and install the following packages:
```
conda create --name=biodiversity six pandas ipython-notebook scikit-learn git basemap matplotlib xlrd numba gdal rasterio python=3
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


