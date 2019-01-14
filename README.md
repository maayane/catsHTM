# catsHTM
The `catsHTM` package is a tool for fast accessing and cross-matching large astronomical catalogs, originally written in `Matlab` by Eran O. Ofek. Here we present the Python version. 

[![PyPI](https://img.shields.io/pypi/v/catsHTM.svg?style=flat-square)](https://pypi.python.org/pypi/catsHTM)

```python
>>> import catsHTM
>>> catsHTM.cone_search('FIRST',0,0,500)
```
## Documenation

The HDF5/HTM format, designed to store and provide fast access for large astronomical catalogs (with >10^6 rows) is described in the [preliminary documentation](https://webhome.weizmann.ac.il/home/eofek/matlab/doc/catsHTM.html), together with the `Matlab` version.

The `catsHTM` package is also described in a paper by [Soumagnac & Ofek 2018](https://arxiv.org/abs/1805.02666).

## Credit
If you are using one of the large catalogs, or this tool, please give the specific reference and acknowledgments to the catalogs you used (see [list](https://webhome.weizmann.ac.il/home/eofek/matlab/doc/catsHTMcredit.html)) and add the following acknowledgment:

*The XXX catalog we use was formatted into the HDF5/HTM large catalog format as described in Soumagnac & Ofek (2018) and was developed as part of the Matlab Astronomy & Astrophysics Toolbox (Ofek 2014; ascl.soft 07005).*

[Bibtex entry for Soumagnac & Ofek 2018](http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2018arXiv180502666S&data_type=BIBTEX&db_key=PRE&nocookieset=1)

## How to obtain the formatted catalogs?

The catalogs are available from:
* Web download is now available. See instructions [here](http://euler1.weizmann.ac.il/catsHTM/)
* If you have any question or encounter any prblem while trying to download, email **eran dot ofek at weizmann dot ac dot il** or **maayane dot soumagnac at weizmann dot ac dot il** .

## How to install the `catsHTM` code?

These instruction are for installing the catsHTM **code**, i.e. do not include the installation of the catalogs in HDF5 format. In order to download the catalogs, see section on 'How to obtain the formatted catalogs'.

### pip

`pip install catsHTM`

### Python version
* `python 2`: higher than `2.7.10`
* `python 3`

### Required python packages
* `math`
* `numpy`
* `scipy`
* `h5py`

## How to make a cone search with ``catsHTM``?

First, you need to specify the path to the directory where the HDF5 formatted catalogs where downloaded (default is `./data`). This will only work if you have previously downloaded the catalogs in HDF5 format (in order to download them, see section on ['How to obtain the formatted catalogs'](https://github.com/maayane/catsHTM#how-to-make-a-cone-search-with-catshtm)):

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
[[  6.28300408e+00  -7.24311612e-05   1.68128476e-01   7.59999990e-01
    6.70588255e-01   9.23669338e-02   4.30000019e+00   0.00000000e+00
    6.20000000e+01   7.13999987e+00   4.32999992e+00   5.32000008e+01
    2.45002071e+06   2.45247496e+06]]
```
If you set the optionnal argument verbose to be true, you will get a line with a summary of your search:
```python
>>> cat,colcell, colunits=catsHTM.cone_search('FIRST',0,0,100,catalogs_dir=path, verbose=True)
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
>>> print colunits
[u'rad' u'rad' ' ' u'mJy' u'mJy' u'mJy' u'arcsec' u'arcsec' u'deg'
 u'arcsec' u'arcsec' u'deg' u'MJD' u'MJD']
```
## How to cross-match two catalogs with ``catsHTM``?

First, you need to specify the path to the directory where the HDF5 formatted catalogs where downloaded (default is `./data`). This will only work if you have previously downloaded the catalogs in HDF5 format (in order to download them, see section on ['How to obtain the formatted catalogs'](https://github.com/maayane/catsHTM#how-to-make-a-cone-search-with-catshtm)):

```python
>>> import catsHTM
>>> path='path/to/directory'

```
You then need to call the `xmatch_2cats` function. For example, to look for overlaps between the `FIRST` and `NVSS` catalogs:

```python
>>> catsHTM.xmatch_2cats('FIRST','NVSS',catalogs_dir=path)

Catalog_1 is FIRST (43688 trixels)
Catalog_2 is NVSS (43688 trixels)
************** I am building all the trixels relevant to our search **************
The number of trixels in the highest level, for FIRST is 32768
The number of trixels in the highest level, for NVSS is 32768
************** I am looking for overlapping trixels **************
I am looking for Catalog_2 (NVSS) trixels overlapping with the non-empty trixel #10921 of Catalog_1 (FIRST)
...
```

By default, this will create a directory `./cross-matching_results`, where it will save three files:
1. `cross-matching_result_[name of catalog 1].txt`: the catalog entries of catalog 1 (e.g. FIRST) for which one or mors counterparts were found in catalog 2 (e.g. NVSS), within the search radius.
2. `cross-matching_result_[name of catalog 2].txt`: the catalog entries corresponding to the closest counterpart found in catalog 2 (e.g. NVSS)
3. `cross-matching_result_full.txt`: a file where the two above files were merged.

The header of all these files specify the catalog columns. Examples of such files, obtained when running the code for `FIRST` and `NVSS` as in the commands above can be found in the directory `cross-matching_results_test/`.

You can modify the location of the output files with the `output` keyword:
```python
>>> catsHTM.xmatch_2cats('FIRST','NVSS',catalogs_dir=path)
```
You can also choose to only save the two separate files (1. and 2. in the list above), by setting the `save_in_one_file` keyword to `False`, or save only the large file (3. in the list above) by setting the `save_in_separate_files` keyword to `False`. E.g.:

```python
>>> catsHTM.xmatch_2cats('FIRST','NVSS',catalogs_dir=path,save_in_one_file=False)
```
You can modify the search radius (default is 2 arcsec) with the `Search_radius` keyword. E.g.

```python
>>> catsHTM.xmatch_2cats('FIRST','NVSS',Search_radius=5)
```

### Specification of other default parameters 

Other default parameters, such as the files and datasets naming format, can be edited in the python file `params.py` (type `pip show catsHTM` in the comand line to see where `params.py` is stored).

## What are the available catalogs?
Currently, the following catalogs are available in this format (alphabetical order):

* 2MASS (input name: `TMASS`)
* 2MASSxsc  (input name: `TMASSxsc`) - 2MASS extended source catalog
* AKARI  (input name: `AKARI`)
* APASS  (input name: `APASS`) - AAVSO All Sky Photometric Sky Survey (~5.5x10^7 sources)
* Cosmos (input name: `Cosmos`) - Sources in the Cosmos field
* DECaLS (input name: `DECaLS`) - DECaLS DR5 release
* FIRST (input name: `FIRST`) - (~9.5x10^5 sources)
* GAIA/DR1 (input name: `GAIADR1`) -  (~1.1x10^9 sources).
* GAIA/DR2 (input name: `GAIADR2`) - NEW! (~1.6x10^9 sources)
* GALEX (input name: `GALEX`) -  GALAEX/GR6Plus7 (~1.7x10^8 sources).
* HSC/v2 (input name: `HSCv2`)- Hubble source catalog
* IPHAS/DR2 (input name: `IPHAS`)
* NED redshifts (input name: `NEDz`)
* NVSS (input name: `NVSS`) - (~1.8x10^6 sources)
* PS1 (input name: `PS1`) - Pan-STARRS (~2.6x10^9 sources; A cleaned version of the PS1 stack catalog; some missing tiles below declination of zero [being corrected])
* PTFpc (input name: `PTFpc`) - PTF photometric catalog 
* ROSATfsc (input name: `ROSATfsc`) - ROSAT faint source catalog
* SDSS/DR10 (input name: `SDSSDR10`)- Primary sources from SDSS/DR10 (last photometric release)
* Skymapper - will be added soon.
* SpecSDSS/DR14 (input name: `SpecSDSS`) - SDSS spectroscopic catalog
* Spitzer/SAGE (input name `SAGE`)
* Spitzer/IRAC (input name `IRACgc`) - Spitzer IRAC galactic center survey
* UCAC4 (input name: `UCAC4`) - (~1.1x10^8 sources)
* UKIDSS/DR10 (input name: `UKIDSS`)
* USNOB1 (not yet available)
* VISTA/Viking/DR3 (not yet available)
* VST/ATLAS/DR3 (input name: `VSTatlas`)
* VST/KiDS/DR3 (input name: `VSTkids`)
* WISE (input name: `WISE`) - ~5.6x10^8 sources
* XMM (input name: `XMM`)- 7.3x10^5 sources 3XMM-DR7 (Rosen et al. 2016; A&A 26, 590)

## How to download the `catsHTM` catalogs?

The catsHTM directory is very large and therefore available on request. The HDF5/HTM catalogs requires about 1.6TB of disk space.

Data is available from:

* Web download (link will be provided soon)
* Shared (and updated automatically) via Dropbox on request from eran dot ofek at weizmann dot ac dot il.


The catalog format is based on the HDF5 file format and HDF5 file access utilities, which are available on many platforms. The catalog format is designed to allow fast access for cone searches in the range of 1 arcsec to about 1 deg. For fast access, the sources are sorted into hierarchical triangular mesh (HTM). 


