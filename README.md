# Color from CRISM
Python scripts to replicate human and spacecraft color vision from CRISM images.

## Requirements
- requirements.txt

## Description

This Python package converts a Compact Reconnaissance Imaging Spectrometer for Mars (CRISM) data cube to various color spaces, including human perceptual color and emulated spacecraft filters. It integrates across a given wavelength range using the CIE color matching functions (for human perceptual color) or spacecraft filter functions (for spacecraft color) before conversion into sRGB color space. Currently, this package is only capable of processing Map-projected Targed Reduced Data Records (MTRDR), which represent the highest level of processing by the CRISM team. MTRDR images are map-projected, have the instruments 3 detectors joined into a single image, and are processed to reduce signal from atmospheric features (water ice/dust) and instrumental artifacts.

This code was developed to aid visualization of hyperspectral imaging data. It is free for personal use and academic presentations and publications. Please provide an acknowledgement in your visualization/presentation/publication when using this work.

## Setup and Use

Install the dependent packages with pip. Copy all files to a new directory. This program uses the Fire package to implement a command-line interface. To use it, open a terminal in the directory where you copied the files. Type `python crism.py [function] --arg1=[$1] --arg2=[$2]` to run a command. 

### For Human Perceptual Color

The `mtrdr_to_color()` function uses the CIE color matching functions as a close approximation to standard human vision. Additional keyword arguments can be used to apply the matching function to a user-specified wavelength range, emulating the standard human visual response for that wavelength range. 

User input: `python crism.py mtrdr_to_color --file="[cube_name].lbl" --name="[output_name]"`. The standard outputs from this function include the following:

- 'VIS' - CIE human visual response, covering the wavelength range from 380 - 780 nm
- 'FAL' - Human visual response shifted to 1.1 - 2.2 microns, covering range of CRISM FAL parameter.
- 'FEM' - Human visual response shifted to 750 - 1.2 microns, capturing color variability in iron oxidation and Fe-rich mineralogy.
- 'MAF' - Human visual response shifted to 800 - 2.0 microns, capturing color variability in unweathered basaltic minerals.
- 'PHY' - Human visual response shifted to 1.8 - 2.2 microns, with color strongly dependent on presence of H20 and OH vibrational bands (hydrated minerals). 
- 'FAR' - Human visual response shifted to 2.8 - 3.9 microns, spanning range of CRISM longwave sensor. Sensor prone to light leaks, not all images are usable. 
- 'CAR' - Human visual response shifted to 2.9 - 3.4 microns, color somewhat dependent on presence of carbonate minerals. Same caveat as "FAR" applies. 

If you would not like the above standard outputs, add `--standard_params=False` to the mtrdr_to_color inputs. 

If you would like to specify a custom wavelength range, add `--new_params=[[wave1, wave2], [wave1, wave2]...]` to the mtrdr_to_color inputs. To find wavelength ranges which may highlight interesting mineralogy, see [Vivano-Beck (2014)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JE004627), which discusses the spectral signals of various common Mars minerals.

### For Spacecraft Filter Color
***Work in Progress***
Several functions are provided to calculate color data returned from various spacecraft. As of the current version, color information can be calculated for the following spacecraft instruments.

- Mars Reconnaissance Orbiter HiRISE
- Mars Express HRSC
- Curiosity MastCam

This functionality is still experimental, and may not be radiometrically-accurate. In addition, function control flow is still being worked out.

#### Mars Express HRSC
To calculate HRSC images, use `mtrdr_to_hrsc --file="[cube_name].lbl" --fname="[output_name]"`. Standard output returns an IGB (NIR-GRN-BLU) color image and grayscale images through the spacecraft's Nadir, NIR, RED, GRN, BLU, P1, and S1 filters. 

To change color image combinations, add `--color=[kwarg]` to the command. Kwargs are "IGB" (NIR-GRN-BLU), "IRB" (NIR-RED-BLU), and "RGB" (RED-GRN-BLU)

If grayscale images are not desired, add `--lumin=False` to the command.

#### Mars Reconnaissance Orbiter HiRISE
To calculate HiRISE images, use `mtrdr_to_hirise --file="[cube_name].lbl" --fname="[output_name]"`. Standard output returns an IRB (NIR-RED-BGR) color image.

To output an RGB color product, add `--color="RGB"` to the command.

#### Curiosity MastCam
To calculate MastCam images, use `mtrdr_to_mastcam --file="[cube_name].lbl" --fname="output_name"`. Standard output returns an RGB image equivalent to standard bayer-filtered MastCam color images and grayscale images equivalent to images taken through narrowband filters L1-L6 and R1-R6. 

If narrowband filter images are not desired, add `--narrowband=False` to the `mtrdr_to_mastcam` inputs.


## Future improvements

- Continue implementing color-matching functions for various spacecraft filters, double-check to make sure method is radiometrically accurate.

## Acknowledgements
This code makes heavy use of the 'ColourSystem' class [described on the SciPython blog](https://scipython.com/blog/converting-a-spectrum-to-a-colour/). It also uses the [SpectRes package](https://spectres.readthedocs.io/en/latest/) to convert CIE matching functions to different wavelength ranges. 
