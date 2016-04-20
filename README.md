# iSDM
Towards a species-by-species approach to global biodiversity modelling.


Step 1: Get Anaconda
```
wget http://repo.continuum.io/archive/Anaconda3-2.4.0-Linux-x86_64.sh  # Linux
bash Anaconda3-2.4.0-Linux-x86_64.sh  # Say 'yes' to appending to .bashrc and specify the installation directory
```

Step 2: install packages
```
conda install six pandas ipython-notebook scikit-learn numpy pygbif basemap rpy2 matplotlib xlrd numba python=3
```

Step 2.1 install pygbif
```
pip install pygbif
```

Step 3: install iSDM 
```
python setup.py install
```

