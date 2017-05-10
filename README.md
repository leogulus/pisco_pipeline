# Pisco Pipeline

This is a current version of PISCO data reduction pipeline, developed by [Taweewat Somboonpanyakul](http://leogulus.github.io/) at MIT in order to reduce PISCO raw images to calibrated images with calibrated photometry.

## Installation:

Required softwares before using the pipeline
- [Astrometry.net](http://astrometry.net/use.html) (including necessary index files for your field of interest) 
- [Sextractor](http://www.astromatic.net/software/sextractor)
- [SCAMP](https://www.astromatic.net/software/scamp)
- [SWARP](https://www.astromatic.net/software/swarp)
- *comics.py* python package (developed from [LA-Cosmic](http://www.astro.yale.edu/dokkum/lacosmic/). The package is also included here. 

## Usage:

```python
python pisco_pipeline/pisco_combine.py data/ Field024
python pisco_pipeline/pisco_photometry.py Field026
```
where 'data/' is the directory wherer the data Field024_A_82.fits and Field024_A_83.fits are located, and Field024 is the prefix for all the files that you want to combine together. 

Run python script outside pisco_pipeline/ directory where data/ directory is also located.

## License

MIT
