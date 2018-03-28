# catsHTM
catsHTM is a tool for fast accessing and cross-matching of large astronomical catalogs, originally written in Matlab by Eran O. Ofek. Here we present the Python version. 

```python
>>> import catsHTM
>>> catsHTM.cone_search('FIRST',0,0,500)
```
## Full Documenation
Preliminary documentation, as well as the Matlab version are available here: https://webhome.weizmann.ac.il/home/eofek/matlab/doc/catsHTM.html

## Installation
`pip install catsHTM`

### Python version
* `python 2`: higher than `2.7.10`
* a `python 3` version will soon be available.

### Required python packages
* `math`
* `numpy`
* `scipy`
* `h5py`


## How to make a cone search with ``catsHTM``?

First, you need to specify the path to the directory where the HDF5 formatted catalogs where downloaded (default is `./data`):
```python
>>> import catsHTM
>>> path='path/to/directory'
```
You can then call the `cone_search` function. For example, to look for the sources in a cone of 100 arcsec centered on RA=0 rad and DEC=0 rad, in the FIRST catalog, type: 
```python
>>> cat,colcell, colunits=catsHTM.cone_search('FIRST',0,0,100,catalogs_dir=path)
```
The catalog lines corresponding to the sources within the cone are stored in the `numpy` array `cat`:
```python
>>> print cat
*************
Catalog: FIRST; cone radius: 100 arcsec; cone center: (RA,DEC)=(0,0)
*************
[[  6.28300408e+00  -7.24311612e-05   1.68128476e-01   7.59999990e-01
    6.70588255e-01   9.23669338e-02   4.30000019e+00   0.00000000e+00
    6.20000000e+01   7.13999987e+00   4.32999992e+00   5.32000008e+01
    2.45002071e+06   2.45247496e+06]]
```
The names of the catalog columns are stored in the `numpy` array `colcell`

```python
>>> print colcell
['RA' 'Dec' 'SideProb' 'Fpeak' 'Fint' 'rms' 'Major' 'Minor' 'PosAng'
 'FitMajor' 'FitMinor' 'FitPosAng' 'StartMJD' 'StopMJD']
```
The units of the catalog columns are stored in the `numpy` array `colunits`

```python
>>> print colcell
[u'rad' u'rad' ' ' u'mJy' u'mJy' u'mJy' u'arcsec' u'arcsec' u'deg'
 u'arcsec' u'arcsec' u'deg' u'MJD' u'MJD']
```

### Specification of other default parameters 

Other default parameters, such as the files and datasets naming format, can be edited in the python file `params.py` (type `pip show catsHTM` in the comand line to see where `params.py` is stored).

