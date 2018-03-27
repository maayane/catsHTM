# catsHTM
catsHTM is a tool for fast accessing and cross-matching of large astronomical catalogs, originally written in Matlab by Eran O. Ofek. Here we present the Python version. 

## Full Documenation
Preliminary documentation as well as the Matlab version are available here: https://webhome.weizmann.ac.il/home/eofek/matlab/doc/catsHTM.html

## Installation
`pip install catsHTM`

### Python version
* `python 2`: higher than 2.7.10
* a `python 3` version will soon be available.

### Required python packages
* math
* numpy
* scipy
* h5py

## Before running catsHTM ...
### Specification of the default parameters 

Before running, specify the path where the HDF5 catalogs are stored, in `params.py`. 
Other default parameters can be edited (but do not have to, for the code to run) in `params.py`

### Example 

An example is given in the code example.py
To run it:
```python
python example.py
```
