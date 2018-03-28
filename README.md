# catsHTM
The `catsHTM` package is a tool for fast accessing and cross-matching large astronomical catalogs, originally written in `Matlab` by Eran O. Ofek. Here we present the Python version. 

```python
>>> import catsHTM
>>> catsHTM.cone_search('FIRST',0,0,500)
```
## Documenation

The HDF5/HTM format, designed to store and provide fast access for large astronomical catalogs (with >10^6 rows) is described in the [preliminary documentation](https://webhome.weizmann.ac.il/home/eofek/matlab/doc/catsHTM.html), together with the `Matlab` version.

The `catsHTM` package is also described in a paper by Soumagnac & Ofek (2018, in preparation).

## Credit
If you are using one of the large catalogs, or this tool, please give the specific reference and acknowledgments to the catalogs you used (see [list](https://webhome.weizmann.ac.il/home/eofek/matlab/doc/catsHTMcredit.html)) and add the following acknowledgment:

*The XXX catalog we use was formatted into the HDF5/HTM large catalog format as described in Soumagnac & Ofek (2018) and was developed as part of the Matlab Astronomy & Astrophysics Toolbox (Ofek 2014; ascl.soft 07005).*

## How to install `catsHTM`?

### pip

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


## How to download the `catsHTM` catalogs?

The catsHTM directory is very large and therefore available on request. The HDF5/HTM catalogs requires about 1.6TB of disk space.

Data is available from:

* Web download (link will be provided soon)
* Shared (and updated automatically) via Dropbox on request from eran dot ofek at weizmann dot ac dot il.


The catalog format is based on the HDF5 file format and HDF5 file access utilities, which are available on many platforms. The catalog format is designed to allow fast access for cone searches in the range of 1 arcsec to about 1 deg. For fast access, the sources are sorted into hierarchical triangular mesh (HTM). Currently, the following catalogs are available in this format (alphabetical order):

* 2MASS
* 2MASSxsc - 2MASS extended source catalog
* AKARI - 
* APASS - AAVSO All Sky Photometric Sky Survey (~5.5x10^7 sources)
* Cosmos - Sources in the Cosmos field
* DECaLS - DECaLS/DR5
* FIRST - (~9.5x10^5 sources)
* GAIA/DR1 -  (~1.1x10^9 sources).
* GALEX -  GALAEX/GR6Plus7 (~1.7x10^8 sources).
* HSC (not yet available)
* IPHAS/DR2 
* NVSS - (~1.8x10^6 sources)
* PS1 - Pan-STARRS (~2.6x10^9 sources; A cleaned version of the PS1 stack catalog; some missing tiles below declination of zero [being corrected])
* ROSATfsc - ROSAT faint source catalog
* SDSSDR10 - SDSS/DR10 - Primary sources from SDSS/DR10 (last photometric release)
* SDSS/DR14offset - RA/Dec offsets for SDSS primary sources with z<20 mag.
* SpecSDSS/DR14 - SDSS spectroscopic catalog
* UCAC4  (~1.1x10^8 sources)
* UKIDSS/DR10
* USNOB1 (not yet available)
* VISTA/Viking/DR3
* VST/ATLAS/DR3
* VST/KiDS/DR3
* WISE (~5.6x10^8 sources)
* XMM - 7.3x10^5 sources) 3XMM-DR7 (Rosen et al. 2016; A&A 26, 590)
