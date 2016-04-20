# iSDM
Towards a species-by-species approach to global biodiversity modelling.


Step 1: Download [Anaconda](https://www.continuum.io/downloads) for Python 3.

NOTE for Windows:
Do not install as Administrator unless admin privileges are required. If you encounter any issues during installation, please temporarily disable your anti-virus software during install, then immediately re-enable it.

Step 2: Create a python environment and install the following packages:
```
conda create --name=biodiversity six pandas ipython-notebook scikit-learn numpy git basemap rpy2 matplotlib xlrd numba python=3
source activate biodiversity
```

Step 2.1 install pygbif
```
pip install pygbif
```

Step 3: install iSDM 
```
git clone https://github.com/remenska/iSDM.git
cd iSDM/
python setup.py install
```

