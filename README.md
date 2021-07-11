# Color from CRISM
Python scripts to replicate human color vision from CRISM images.

## Requirements
- requirements.txt

## Description

This Python package converts a Compact Reconnaissance Imaging Spectrometer for Mars (CRISM) data cube to human perceptual color. It integrates across a given wavelength range using the CIE color matching functions before conversion into sRGB color space. For the visible wavelength range, this is a close approximation to standard human vision. For other wavelength ranges, this emulates the standard human visual response for that wavelength range. Currently, this package is only capable of processing Map-projected Targed Reduced Data Records (MTRDR), which represent the highest level of processing by the CRISM team. MTRDR images are map-projected, have the instruments 3 detectors joined into a single image, and are processed to reduce signal from atmospheric features (water ice/dust) and instrumental artifacts.

This code was developed to aid visualization of hyperspectral imaging data. It is free for personal use and academic presentations and publications. Please provide an acknowledgement in your visualization/presentation/publication when using this work.

## Setup and Use

Install the dependent packages with pip. Copy all files to a new directory. This program uses the Fire package to implement a command-line interface. To use it, open a terminal in the directory where files are stored and type `python crism.py --file=[cube name].lbl --name=[output name].png`

- 'VIS' - Similar to the official CRISM VIS parameter, covering the wavelength range from 380 - 780 nm
- 'FAL' - Similar to the official CRISM FAL parameter, covering the wavelength range from 1.1 - 2.2 microns.
- 'FEM' - "Fe-mineral" capturing the 750 nm - 1.2 micron range due to variations in iron oxidation and mineralogy.
- 'MAF' - "Mafic minerals" capturing the 800 nm - 2.0 micron range, primarily due to variations in unweathered basaltic minerals.
- 'PHY' - "Phyllosilicates" capturing the 1.8 - 2.2 micron range, with color strongly dependent on H20 and OH vibrational bands.
- 'FAR' - Integrates across the 2.8 - 3.9 micron spectral range of the CRISM longwave sensor. Sensor was prone to light leaks, not all images are usable.
- 'CAR' - "Carbonates" capturing the 2.9 - 3.4 micron spectral range, color somewhat dependent on the presence of carbonates. Same caveat as above.

If you would not like the above standard outputs, add `--standard_params=False` to the mtrdr_to_color inputs. 

If you would like to use your own wavelength ranges, add `--new_params=[[wave1, wave2], [wave1, wave2]...]` to the mtrdr_to_color inputs. To find wavelength ranges which may highlight interesting mineralogy, see [Vivano-Beck (2014)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JE004627), which discusses the spectral signals of various common Mars minerals.

## Future improvements

- Develop color-matching functions for various spacecraft filters.

## Acknowledgements
This code makes heavy use of the 'ColourSystem' class [described on the SciPython blog](https://scipython.com/blog/converting-a-spectrum-to-a-colour/). It also uses the [SpectRes package](https://spectres.readthedocs.io/en/latest/) to convert CIE matching functions to different wavelength ranges. 

