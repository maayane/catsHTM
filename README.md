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

The catalog format is based on the HDF5 file format and HDF5 file access utilities, which are available on many platforms. The catalog format is designed to allow fast access for cone searches in the range of 1 arcsec to about 1 deg. For fast access, the sources are sorted into hierarchical triangular mesh (HTM).
The HDF5/HTM catalogs requires about 1.6TB of disk space.

The catalogs are available from:
* Web download is now available. See instructions [here](http://euler1.weizmann.ac.il/catsHTM/)
* If you have any question or encounter any prblem while trying to download, email **eran dot ofek at weizmann dot ac dot il** or **maayane dot soumagnac at weizmann dot ac dot il** .

## How to install the `catsHTM` code?

These instruction are for installing the catsHTM **code**, i.e. do not include the installation of the catalogs in HDF5 format. In order to download the catalogs, see section on 'How to obtain the formatted catalogs'.

### pip

`pip install catsHTM`

### Python version
* `python 2`: higher than `2.7.10`
* `python 3` (required for the cross-matcher)

### Required python packages
* `math`
* `numpy`
* `scipy`
* `h5py`
* `tqdm`

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
['rad' 'rad' ' ' 'mJy' 'mJy' 'mJy' 'arcsec' 'arcsec' 'deg'
 'arcsec' 'arcsec' 'deg' 'MJD' 'MJD']
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

**To save the results, you need to specify `save_results==True`!**

By default, this will create a directory `./cross-matching_results`, where it will save three files:
1. `cross-matching_result_[name of catalog 1].txt`: the catalog entries of catalog 1 (e.g. FIRST) for which one or mors counterparts were found in catalog 2 (e.g. NVSS), within the search radius.
2. `cross-matching_result_[name of catalog 2].txt`: the catalog entries corresponding to the closest counterpart found in catalog 2 (e.g. NVSS)
3. `cross-matching_result_full.txt`: a file where the two above files were merged.

The header of all these files specify the catalog columns. Examples of such files, obtained when running the code for `FIRST` and `NVSS` as in the commands above can be found in the directory `cross-matching_results_test/`.

You can modify the location of the output files with the `output` keyword:
```python
>>> catsHTM.xmatch_2cats('FIRST','NVSS',catalogs_dir=path)
```

You can speed-up the run by leaving `save_results` to its default value (`False`), e.g., if you do not need to save the results and would rather [use the output of the cross-matching algorythm in your own function](https://github.com/maayane/catsHTM/blob/master/README.md#getting-information-on-the-cross-matched-sources-and-counterparts-and-running-your-own-function-on-the-outputs-of-the-cross-matcher).

You can also choose to only save the two separate files (1. and 2. in the list above), by setting the `save_in_one_file` keyword to `False`, or save only the large file (3. in the list above) by setting the `save_in_separate_files` keyword to `False`. E.g.:

```python
>>> catsHTM.xmatch_2cats('FIRST','NVSS',catalogs_dir=path,save_in_one_file=False)
```
You can modify the search radius (default is 2 arcsec) with the `Search_radius` keyword. E.g.

```python
>>> catsHTM.xmatch_2cats('FIRST','NVSS',Search_radius=5)
```

### Getting information on the cross-matched sources and counterparts and running your own function on the outputs of the cross-matcher ###

The `QueryAllFun` and `QueryAllFunPar` keywords allow you to define a function to be ran on the outputs of the cross-matching algorythm. `QueryAllFun` takes the following input arguments:
1. `Cat1`: the content of a trixel of catalog 1.
2. `Cat2`: the content of a trixel of catalog 2 overlapping with `Cat1`.
3. `Ind`: a list of dictionnaries, with one dictionnary per `Cat1`'s object having one or more counterparts in `Cat2`:
* `Ind[i]["IndRef"]`: the index of the i-th `Cat1`'s source having one or more counterpart in `Cat2`
* `Ind[i]["IndCat"]`: after sorting Cat2 by declination (``` cat2=Cat2[Cat2[:, 1].argsort(),] ``` if DEC is the second column of the catalog), `Ind[i]["IndCat"]` is the list of indices of the counterparts of the matched `Cat1`'s source. 
* `Ind[i]["Dist"]`: a vector of angular distances (radians) between the i-th `Cat1`'s source and its counterparts in `Cat2`.
4. `IndCatMinDist`: a vector, with as many elements as lines in `Cat1`, with 'nan' at lines where there is no counterpart in `Cat2`, and at line where there is, the index - in `Cat2` - of the closest counterpart.

You can write a `QueryAllFun` function e.g. to save or use one or all of the above informations.

An example of such a function, `Example_QueryAllFun`, is built in the code. Simply structure yours in the same way:
```python
def Your_QueryAllFun(Cat1,Ind,Cat2,IndCatMinDist,i,additionnal_args=[1,2,'hi']):
    return [your output]
```
The `QueryAllFunPar` keyword, if not `None`, must be a tuple which will be passed to the `additional_args` keyword of `Your_QueryAllFun`.

For example, runing the code with `QueryAllFun=Example_QueryAllFun` and `QueryAllFunPar=['test']` allows you to print `Cat1`, `Cat2`, `Ind` and `IndCatMinDist` and save the content of `Cat1` in a directory called `test`:

```python
>>> import catsHTM
>>> from catsHTM import Example_QueryAllFun
>>> catsHTM.xmatch_2cats('FIRST','APASS',catalogs_dir=path,QueryAllFun=Example_QueryAllFun,QueryAllFunPar=['test'])
```

### Example: cross-matching two catalogs and selecting objects based on a criterion on both catalogs

In the following example, we cross-match the GAIA and PS1 catalog, and select only objects based on a criterion on the GAIA parallax and the PS1 g-r color.

```python
import catsHTM
import math
import numpy as np
import pandas as pd
import os

def query_function_example(Cat1,Ind,Cat2,IndCatMinDist,i,additionnal_args):

    columns_in_crossmatchfile=additionnal_args[0]+additionnal_args[1]
    Verbose=additionnal_args[2]

    Nind=len(Ind)
    if Verbose==True:
        print('Nind',Nind)

    data_list=[]
    for si in range(Nind):
        if Verbose==True:
            print('Ind',Ind)
        I1=Ind[si]['IndRef'] # (Index of a Cat1's source having one or more counterpart in Cat2)
        I2 = int(Ind[si]['IndCat'][np.argmin(Ind[si]['Dist'])])
        if Verbose==True:
            print('IndCatMinDist',IndCatMinDist)
            print("Ind[si]['IndRef']",Ind[si]['IndRef'])
            print("Ind[si]['IndCat']",Ind[si]['IndCat'])
            print('I2, the index of the closest is',I2) # List of indexes of the Cat2's counterparts.

        gaia3_parallax     = cat1_cols.index('Plx')
        gaia3_parallax_err = cat1_cols.index('ErrPlx')
        ps1_g              = cat2_cols.index('gPSFMag')
        ps1_r              = cat2_cols.index('rPSFMag')

        parallax_constraint = 15   # constraint on the parallax: Parallax/Err > 15
        g_r_constraint      = -0.5 # constraint on the g-r color: (g-r) < -0.5
        
        
        # the criterion based on which we want to select the matches:
        cond=(Cat1[I1,gaia3_parallax]/Cat1[I1,gaia3_parallax_err] > parallax_constraint) & \
             (Cat2[I2,ps1_g]-Cat2[I2,ps1_r] < g_r_constraint) & \
             (Cat2[I2,ps1_g]-Cat2[I2,ps1_r] > -5.0 )

        if cond:
            if Verbose==True:
                print('Success.')
            data_line = {}
            for bli,blu in enumerate(np.hstack((Cat1[I1,:],Cat2[I2,:]))):
                data_line[columns_in_crossmatchfile[bli]]=blu
            data_list.append(data_line)
    data = pd.DataFrame(data_list)

    if len(data)>0:
        data['RA_gaia3'] *= 180./np.pi # Convert RA,Dec back to decimal degrees before saving.
        data['Dec_gaia3'] *= 180./np.pi
        data['RA_ps1'] *= 180./np.pi
        data['Dec_ps1'] *= 180./np.pi
        if not os.path.isfile("data_xmatch.csv"):
            # If the output file does not exist, then create it and save the header and data.
            data.to_csv('data_xmatch.csv'.format(i),index=False, mode='a', header=True)
        else:
            # Else if the output file does exist, then do not print the header line; only print the data.
            data.to_csv('data_xmatch.csv'.format(i),index=False, mode='a', header=False)
    return data
```
Then we run the xmatch_2cats function with this:

```

    cat_path = "the/directory/where/you/have/your/HDF5/GAIA/and/PS1/catalogs"
    verbose=False

    galex_columns = ['RA_galex', 'Dec_galex', 'nuv_mag', 'nuv_magerr', 'fuv_mag', 'fuv_magerr', 'NUV_FWHM_WORLD', 'nuv_weight', 'fuv_weight']
    gaia_columns = ['RA_gaia3', 'Dec_gaia3', 'Epoch', 'ErrRA_gaia3', 'ErrDec_gaia3', 'Plx', 'ErrPlx', 'PMRA', 'ErrPMRA', 'PMDec', 'ErrPMDec', 'RA_Dec_Corr',
                    'NobsAst', 'ExcessNoise', 'ExcessNoiseSig', 'Chi2Ast', 'DofAst', 'NGphot', 'Mag_G', 'ErrMag_G', 'Mag_BP', 'ErrMag_BP',
                    'Mag_RP', 'ErrMag_RP', 'BPRP_Excess', 'RV', 'ErrRV', 'Teff', 'LogG', 'FeH']
    ps1_columns = ['RA_ps1', 'Dec_ps1', 'ErrRA_ps1', 'ErrDec_ps1', 'MeanEpoch', 'posMeanChi2', 'gPSFMag', 'gPSFMagErr', 'gpsfLikelihood', 'gMeanPSFMagStd',
                   'gMeanPSFMagNpt', 'gMeanPSFMagMin', 'gMeanPSFMagMax', 'rPSFMag', 'rPSFMagErr', 'rpsfLikelihood', 'rMeanPSFMagStd', 'rMeanPSFMagNpt',
                   'rMeanPSFMagMin', 'rMeanPSFMagMax', 'iPSFMag', 'iPSFMagErr', 'ipsfLikelihood', 'iMeanPSFMagStd', 'iMeanPSFMagNpt', 'iMeanPSFMagMin',
                   'iMeanPSFMagMax', 'zPSFMag', 'zPSFMagErr', 'zpsfLikelihood', 'zMeanPSFMagStd', 'zMeanPSFMagNpt', 'zMeanPSFMagMin', 'zMeanPSFMagMax',
                   'yPSFMag', 'yPSFMagErr', 'ypsfLikelihood', 'yMeanPSFMagStd', 'yMeanPSFMagNpt', 'yMeanPSFMagMin', 'yMeanPSFMagMax']

    cat1_cols = gaia_columns
    cat2_cols = ps1_columns

    catsHTM.xmatch_2cats('GAIAEDR3','PS1',Search_radius=4,QueryAllFun=query_function_example,QueryAllFunPar=[cat1_cols, cat2_cols, verbose],
                     catalogs_dir=cat_path,Verbose=verbose,save_results=False,save_in_one_file=False,
                     save_in_separate_files=False,output='./cross-matching_results',time_it=False,Debug=False)
```

### Specification of other default parameters 

Other default parameters, such as the files and datasets naming format, can be edited in the python file `params.py` (type `pip show catsHTM` in the comand line to see where `params.py` is stored).

## What are the available catalogs?
Currently, the following catalogs are available in this format (alphabetical order; for a list with the number of sources and references see [here](https://euler1.weizmann.ac.il/catsHTM/catsHTM_catalogs.html)):

* 2MASS (input name: `TMASS`)
* 2MASSxsc  (input name: `TMASSxsc`) - 2MASS extended source catalog
* AKARI  (input name: `AKARI`)
* APASS  (input name: `APASS`) - AAVSO All Sky Photometric Sky Survey (~5.5x10^7 sources)
* Cosmos (input name: `Cosmos`) - Sources in the Cosmos field
* DECaLS (input name: `DECaLS`) - DECaLS DR5 release
* FIRST (input name: `FIRST`) - (~9.5x10^5 sources)
* GAIA/DR1 (input name: `GAIADR1`) -  (~1.1x10^9 sources).
* GAIA/DR2 (input name: `GAIADR2`) - NEW! (~1.6x10^9 sources)
* GAIA/EDR3 (input name: `GAIAEDR3`) - NEW! (~1.8x10^9 sources)
* GALEX (input name: `GALEX`) -  GALAEX/GR6Plus7 (~1.7x10^8 sources).
* HSC/v2 (input name: `HSCv2`)- Hubble source catalog
* IPHAS/DR2 (input name: `IPHAS`)
* NED redshifts (input name: `NEDz`)
* NVSS (input name: `NVSS`) - (~1.8x10^6 sources)
* HYPERLEDA (input name: `PGC`) 
* PS1 (input name: `PS1`) - Pan-STARRS (~2.6x10^9 sources; A cleaned version of the PS1 stack catalog; some missing tiles below declination of zero [being corrected])
* The PTF photometric catalog (input name: `PTFpc`)
* ROSATfsc (input name: `ROSATfsc`) - ROSAT faint source catalog
* SDSS/DR10 (input name: `SDSSDR10`)- Primary sources from SDSS/DR10 (last photometric release)
* Skymapper DR1 (input name: `Skymapper`)
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
* ZTF-DR1 stellar variability catalog (input name: `ztfSrcLCDR1`)
* ZTF-DR1 variable star candidates (input name: `ztfSrcLCDR1`)

## How to access the catalog of ZTF variable sources from Ofek et al. 2020 with ``catsHTM``?

The routine `read_ztf_HDF_matched` allows access the catalogs presented in Ofek et al 2020.

The path to the directory where the HDF5 formatted light curves where downloaded needs to be specified first (default is `.`).

The other arguments of the routine are (1) the ZTF field number, (2) an array with the [start end] lines to read. The lines for a given source are available in I1 and I2 in the `ztfSrcLCDR1` catsHTM catalog.

```python
>>> import catsHTM
>>> path='path/to/directory'
>>> Cat,ColCel=catsHTM.read_ztf_HDF_matched(815,[10,25],path=path)
>>> print(Cat)
[[  5.84743309e+04   1.84250000e+01   4.90000000e-02   9.80000000e-02
    0.00000000e+00]
 [  5.84405160e+04   1.83110000e+01   4.50000000e-02   8.90000000e-02
    0.00000000e+00]
 [  5.83755087e+04   1.83670000e+01   4.70000000e-02   9.40000000e-02
    0.00000000e+00]
 [  5.84643027e+04   1.83350000e+01   4.60000000e-02   1.09000000e-01
    0.00000000e+00]
 [  5.84254711e+04   1.83920000e+01   4.80000000e-02   1.02000000e-01
    0.00000000e+00]
 [  5.84434111e+04   1.83710000e+01   4.70000000e-02   8.90000000e-02
    0.00000000e+00]
 [  5.84615229e+04   1.83350000e+01   4.60000000e-02   9.30000000e-02
    0.00000000e+00]
 [  5.83695057e+04   1.84010000e+01   4.80000000e-02   1.01000000e-01
    0.00000000e+00]
 [  5.84562872e+04   1.83250000e+01   4.50000000e-02   8.90000000e-02
    0.00000000e+00]
 [  5.84503844e+04   1.83100000e+01   4.50000000e-02   1.00000000e-01
    0.00000000e+00]
 [  5.82311746e+04   1.83100000e+01   4.50000000e-02   8.10000000e-02
    0.00000000e+00]
 [  5.82211967e+04   1.83870000e+01   4.80000000e-02   1.07000000e-01
    3.27680000e+04]
 [  5.82431606e+04   1.83440000e+01   4.60000000e-02   7.70000000e-02
    0.00000000e+00]
 [  5.84373895e+04   1.83480000e+01   4.60000000e-02   9.70000000e-02
    0.00000000e+00]
 [  5.82241671e+04   1.83120000e+01   4.50000000e-02   7.80000000e-02
    3.27680000e+04]
 [  5.84344523e+04   1.84180000e+01   4.90000000e-02   9.80000000e-02
    0.00000000e+00]]
>>> print(ColCel)
['HMJD' 'Mag' 'MagErr' 'ColorCoef' 'Flags']
```
