# Pisco Pipeline

This is a current version of PISCO data reduction pipeline, developed by [Taweewat Somboonpanyakul](http://leogulus.github.io/) at MIT in order to reduce PISCO raw images to calibrated images with calibrated photometry.

## Installation:

Required softwares before using the pipeline:
- [Astrometry.net](http://astrometry.net/use.html) (including necessary index files for your field of interest)
- [Sextractor](http://www.astromatic.net/software/sextractor)
- [SCAMP](https://www.astromatic.net/software/scamp)
- [SWARP](https://www.astromatic.net/software/swarp)

Included packages:
- Cosmics ray removal python packages: *comics.py* python package (developed from [LA-Cosmic](http://www.astro.yale.edu/dokkum/lacosmic/). The package is also included here, but required to be installed before.
- Stellar Locus Regression python package: [Big-macs-calibrate](https://github.com/patkel/big-macs-calibrate) by Kelly et al. 2014 (http://arxiv.org/abs/1208.0602)

Other Required Python Packages:
- [astropy](www.astropy.org/)
- [photutils](https://photutils.readthedocs.io/) (an affiliated package of Astropy) (required only for the Photometry pipeline)
- numpy
- matplotlib
- subprocess, shlex, os, sys

## Usage:

### Running the Astrometry pipeline
```python
python pisco_pipeline/pisco_combine.py data/ Field026
```
where *data/* is the directory wherer the PISCO raw data Field024_A_82.fits and Field024_A_83.fits are located, and Field024 is the prefix for all the files that you want to combine together. Run python script outside *pisco_pipeline/* directory next to *data/* directory is located.

The main outputs (4 fits images based on 4 different bands) will be located in *final/* directory.

### Running the Photometry pipeline
```python
python pisco_pipeline/pisco_photometry.py Field026
```
Need to run Astrometry pipeline first as the input for the photometry pipeline. The main output for the pipeline (csv file with corrected magnitude in different bands with their uncertainty) is located in *slr_output/* directory.

## License

MIT
