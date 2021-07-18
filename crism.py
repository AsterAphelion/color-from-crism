import rasterio
import numpy as np
import spectres as spec
import fire

#Initialize color_system.py (this segment of code by 'christian' on the SciPython blog)
#See: https://scipython.com/blog/converting-a-spectrum-to-a-colour/

#This is the core of the code, tinker at your own risk.
def xyz_from_xy(x, y):
        """Return the vector (x, y, 1-x-y)."""
        return np.array((x, y, 1-x-y))

class ColourSystem:
    """A class representing a colour system.

    A colour system defined by the CIE x, y and z=1-x-y coordinates of
    its three primary illuminants and its "white point"."""

    # The CIE colour matching function for 380 - 780 nm in 5 nm intervals
    cmf = np.loadtxt('matching_functions/cie-cmf.txt', usecols=(1,2,3))

    def __init__(self, red, green, blue, white):
        """Initialise the ColourSystem object.

        Pass vectors (ie NumPy arrays of shape (3,)) for each of the
        red, green, blue  chromaticities and the white illuminant
        defining the colour system."""

        # Chromaticities
        self.red, self.green, self.blue = red, green, blue
        self.white = white
        # The chromaticity matrix (rgb -> xyz) and its inverse
        self.M = np.vstack((self.red, self.green, self.blue)).T 
        self.MI = np.linalg.inv(self.M)
        # White scaling array
        self.wscale = self.MI.dot(self.white)
        # xyz -> rgb transformation matrix
        self.T = self.MI / self.wscale[:, np.newaxis]

    def xyz_to_rgb(self, xyz, out_fmt=None):
        """Transform from xyz to rgb representation of colour.

        The output rgb components are normalized on their maximum
        value. If xyz is out the rgb gamut, it is desaturated until it
        comes into gamut."""

        rgb = self.T.dot(xyz)
        if np.any(rgb < 0):
            # We're not in the RGB gamut: approximate by desaturating
            w = - np.min(rgb)
            rgb += w

        return rgb

    def spec_to_xyz(self, spec):
        """Convert a spectrum to an xyz point.

        The spectrum must be on the same grid of points as the colour-matching
        function, self.cmf: 380-780 nm in 5 nm steps."""

        XYZ = np.sum(spec[:, np.newaxis] * self.cmf, axis=0)
        den = np.sum(XYZ)
        if den == 0.:
            return XYZ
        return XYZ / den

    def spec_to_rgb(self, spec, out_fmt=None):
        """Convert a spectrum to an rgb value."""

        xyz = self.spec_to_xyz(spec)
        return self.xyz_to_rgb(xyz, out_fmt)

illuminant_D65 = xyz_from_xy(0.3127, 0.3291)
cs_hdtv = ColourSystem(red=xyz_from_xy(0.67, 0.33),
                       green=xyz_from_xy(0.21, 0.71),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)

cs_smpte = ColourSystem(red=xyz_from_xy(0.63, 0.34),
                        green=xyz_from_xy(0.31, 0.595),
                        blue=xyz_from_xy(0.155, 0.070),
                        white=illuminant_D65)

cs_srgb = ColourSystem(red=xyz_from_xy(0.64, 0.33),
                       green=xyz_from_xy(0.30, 0.60),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)


##Defining a few internal functions to help us on our journey.

#Some frequently-used few-liner functions
def find_band(array, value):
    """One-liner to find the index value of the nearest band to a given 
    wavelength value."""
    idx = (np.abs(array - value)).argmin()
    return idx

def quicknorm(data):
    data = (data-np.amin(data))/np.amax(data)
    return data

def calculate_luminance(weights, cube):
    """Function to calculate an image through a filter given the filter transmission properties
    (weights) from a cube."""
    weights = weights/np.sum(weights)
    lumin = np.average(cube, axis=2, weights=weights)    
    return(lumin)

def convert_uint16(cube):
    """Converts cube data (float format) to 16-bit unsigned integer."""
    cube = cube * 65535
    cube = cube.astype(np.uint16)
    return(cube)

#MTRDR pre-processing functions
def modify_mtrdr_axis():
    """Crops the image cube to the given wavelength range."""
    mtrdr_axis = np.genfromtxt("matching_functions/mtrdr_axis.tab", delimiter=",")
    mtrdr_axis = mtrdr_axis[:,2]
    
    #Fill in the gaps where bad bands are present
    #Blue gap - 380-436 nm
    add_waves = np.linspace(377.58, 429.62, num=9)
    add_waves = np.around(add_waves, decimals=2)
    mtrdr_axis = np.concatenate((add_waves, mtrdr_axis))
    
    #NIR bad bands 637-710 nm
    bad_band_fill = np.linspace(637.96, 703.1, num=10)
    bad_band_fill = np.around(bad_band_fill, decimals=2)
    mtrdr_axis = np.insert(mtrdr_axis, 40, bad_band_fill, axis=0)
    
    return(mtrdr_axis)

def mtrdr_crop_bands(image_cube, wave_list):
    """Crops the image cube to the given wavelength range."""
    mtrdr_axis = modify_mtrdr_axis()
    
    short = find_band(mtrdr_axis, wave_list[0])
    long = find_band(mtrdr_axis, wave_list[1])
    crop_cube = image_cube[short:long, :, :]
    
    return(crop_cube)
    
def mtrdr_color_matching(wave_list):
    """Adjusts the CIE color matching function to span the given wavelength range."""
    #Import CIE color matching function
    #Index 0 - wavelengths, Index 1 - red matching function
    #Index 2 - green matching function, Index 3 - blue matching function
    cie_matrix = np.genfromtxt("matching_functions/cie-cmf.txt")

    #Import tab-delimited file of wavelength axis
    mtrdr_axis = modify_mtrdr_axis()
    
    #Find mtrdr axis indices with closest values to user-specified values.
    short = find_band(mtrdr_axis, wave_list[0])
    long = find_band(mtrdr_axis, wave_list[1])
    
    ##Now use normalization to rescale wavelength axis of CIE color matching functions
    #to user-specified wavelength range...
    cie_matrix[:,0] = (mtrdr_axis[long] - mtrdr_axis[short]) / (cie_matrix[-1,0] - cie_matrix[0,0]) * (cie_matrix[:,0]-cie_matrix[-1,0]) + mtrdr_axis[long]
    
    #..then resample CIE function values using MTRDR axis values
    red = spec.spectres(mtrdr_axis[short:long], cie_matrix[:,0], cie_matrix[:,1], fill=0, verbose=False)
    green = spec.spectres(mtrdr_axis[short:long], cie_matrix[:,0], cie_matrix[:,2], fill=0, verbose=False)
    blue = spec.spectres(mtrdr_axis[short:long], cie_matrix[:,0], cie_matrix[:,3], fill=0, verbose=False)
    
    #Concatenate the results
    new_mat = np.stack([red, green, blue], axis=-1)
    return(new_mat)


def format_mtrdr(cube):
    """Prepare MTRDR data cube for color production by filling in missing bands."""
    
    ##Grab the modified MTRDR axis. To document some index values I'm using in this function:
    #Indices 0-9 in this axis represent the 377-436 nm channels, which need to be calculated 
    #through extrapolation and are needed to produce the blue channel in the color output. 
    #Indices 40-50 represent the missing 631-709 nm channels, needed to produce the red channel 
    #in the color output. 
    mtrdr_axis = modify_mtrdr_axis()
    
    ##Extrapolate the missing bands (377-436 nm) necessary for blue. To do this, I am 
    #creating a dummy channel by averaging the values from the first six bands. This creates
    #an array of 9 copies of the average of the first 6 valid bands. Later I will subtract
    #a slope constant subtracted from each band to extrapolate the radiance of each band in
    #this wavelength range.
    interp_channel = np.average(cube[0:6], axis=0)
    interp_channel = np.tile(interp_channel, (9,1,1))
    
    #Band-to-band noise is reduced by calculating the slope from three channel pairs and 
    #averaging the result. This step produces a blue slope for each pixel in the image.
    slope = (cube[2] - cube[6]) / (mtrdr_axis[6] - mtrdr_axis[2])
    slope2 = (cube[1] - cube[5]) / (mtrdr_axis[5] - mtrdr_axis[1])
    slope3 = (cube[0] - cube[4]) / (mtrdr_axis[4] - mtrdr_axis[0])
    slope = (slope + slope2 + slope3) / 3
    
    #Next, we tile the slope so that the array shape matches the number of bands we need to
    #fill in, then multiply the slope by the distance from the first good band. The resulting
    #array is then added to the dummy bands to produce the extrapolated data array.
    slope = np.tile(slope, (9,1,1))
    multiplier = mtrdr_axis[9] - mtrdr_axis[0:9]
    multiplier = multiplier[:, np.newaxis, np.newaxis]
    slope = slope * multiplier
    interp_channel += slope
    
    cube = np.concatenate((interp_channel, cube), axis=0)
    
    #Now repeat the process to fill in the VIS-NIR bad bands.
    interp_channel = np.tile(cube[39], (10,1,1))

    slope = (cube[37] - cube[40])/(mtrdr_axis[50] - mtrdr_axis[39])
    slope2 = (cube[38] - cube[41])/(mtrdr_axis[51]- mtrdr_axis[38])
    slope3 = (cube[39] - cube[42])/(mtrdr_axis[52]- mtrdr_axis[37])
    slope = (slope + slope2 + slope3) / 3

    slope = np.tile(slope, (10, 1, 1))
    multiplier = mtrdr_axis[40:50] - mtrdr_axis[39]
    multiplier = multiplier[:, np.newaxis, np.newaxis]
    slope = slope * multiplier

    interp_channel += slope
    
    #Now to insert the VIS-NIR bad bands into the array. This might be faster and more 
    #memory efficient using np.insert, but my brain hurts trying to figure out that function. 
    #So for now, using the inefficient way.
    
    #Blue side of the bad bands
    short = cube[0:40]
    #Red side of the bad bands
    long = cube[40::]

    intermed = np.concatenate((short, interp_channel), axis=0)
    cube = np.concatenate((intermed, long), axis=0)
    
    return(cube)

def color_from_cube(cube, cs, mode="raw"):
    """Core functionality for calculating human perceptual color from CRISM MTRDR."""
    #Transpose array to put the wavelength axis last - personal preference
    cube = cube.transpose(1,2,0)
    
    #We will lose luminance data once we calculate chromaticity, so before transforming the shape
    #of the data cube, I'm going to calculate luminance images by scaling the brightness of each band
    #by the CIE scaling factor at that band, then integrating across the entire wavelength range.
    
    #I'm doing this step here because when calculating color from spacecraft filters in other functions
    #the data cubes remain three-dimensional. It's easier to run this step here while the cube is still
    #three-dimensional than it is to add a dimensionality argument to calculate_luminance() and specify
    #the dimensionality of the data every time.
    
    weights = cs.cmf.copy()
    
    blu_lumin = calculate_luminance(weights[:,0], cube)
    grn_lumin = calculate_luminance(weights[:,1], cube)
    red_lumin = calculate_luminance(weights[:,2], cube)
    
    #Merge luminance cubes together, then perform a contrast stretch. Adding 2% buffers to the minimum
    #and maximum values to avoid histogram clipping.
    lumin = np.stack((blu_lumin, grn_lumin, red_lumin), axis=0)
    lumin = (lumin - (np.amin(lumin) - (0.02*np.amin(lumin)))) / ((np.amax(lumin) + (0.02*np.amax(lumin))))

    #Now reshape the data array so that it's one dimensional and runs more quickly in the loop that
    #calculates chromaticity values (I'm not sure if the chromaticity calculation can be set up to take
    #advantge of broadcasting ufuncs. 
    rows = cube.shape[0]
    cols = cube.shape[1]
    pixels = rows*cols
    cube = cube.reshape(pixels, cube.shape[2])

    #Create a new cube to handle the interpolated color data
    clone_cube = np.empty((cube.shape[0], 3))

    #Convert the wavelength range to RGB values. If this can be broadcast this would
    #run much more quickly.
    # todo: this operation can be done must faster with numpy and broadcasting ufuncs
    for pixel in range(0, pixels):
        clone_cube[pixel] = cs.spec_to_rgb(cube[pixel])
        
    #When chromaticity values integrate outside of the [0-1] range, they need to be scaled back to 
    #that range to be displayed within the chosen colorspace. The ColourSystem class as written by
    #"Christian" normalized on a per-pixel basis, which destroys relative color information. This was
    #dealt with by removing a normalization statement from the xyz_to_rgb function within the class
    #definition.
    
    #We still need to normalize back to the [0-1 range], which we're doing here with the quicknorm()
    #function. "Raw" normalization preserves the relative color channel brightnesses by simply stretching
    #between the highest chromaticity(typically red) and lowest chromaticity (typically blue) values.
    #"WB" independently normalizes each color channel, similar to the output provided in the official
    #CRISM parameter products. 

    if mode=="raw":
        clone_cube = quicknorm(clone_cube)

    if mode=="wb":
        for channel in range(0, clone_cube.shape[1]):
            clone_cube[:,channel] = quicknorm(clone_cube[:,channel])
        
    #Reshape pixels back to original x,y orientation
    cube = clone_cube.reshape(rows, cols, 3).transpose(2, 0, 1)
    
    #Add luminance data to cube
    cube[0,:,:] = cube[0,:,:] * lumin[0]
    cube[1,:,:] = cube[1,:,:] * lumin[1]
    cube[2,:,:] = cube[2,:,:] * lumin[2]
    
    #Convert to unsigned 16-bit
    cube = convert_uint16(cube)
    
    return(cube)

def mtrdr_to_color(file, name, standard_params=True, new_params=None):
    """Function to produce perceptually-accurate color from CRISM MTRDR data."""

    with rasterio.open(file) as src:
        profile = src.profile
        img = src.read()
    cs = cs_srgb

    #Make null values = 0 so that it doesn't break when doing rgb conversion
    #Also need to convert the null pixels outside of image to 0.
    img[img < 0] = 0
    img[img >= 1] = 0

    img = format_mtrdr(img)

    profile.update(
        dtype = rasterio.uint16,
        count = 3,
        driver = 'PNG'
    )
    
    if standard_params == True:
        process_list = ["VIS", "FAL", "FEM", "MAF", "PHY", "FAR", "CAR"]
        mode_list = ["raw", "raw", "wb", "wb", "wb", "wb", "wb"]
        
        for param, mode in zip(process_list, mode_list):

            if param == "VIS":
                #Imitate VIS browse product summarizing wavelength range from 380 to 780 nm
                wave_range = [380, 780]

            if param == "FAL":
                #Imitate FAL browse product summarizing wavelength range from 1.01 to 2.60 microns
                wave_range = [1010, 2600]

            if param == "FEM":
                #Integrate over 750 nm to 1200nm to capture variability in Fe oxidation state/mineralogy
                wave_range = [750, 1200]

            if param == "MAF":
                #Integrate over 800 nm to 2 micron wavelength range capturing variability in
                #primary basaltic minerals.
                wave_range = [800, 2000]

            if param == "PHY":
                #Integrate over 1.8 to 2.3 micron wavelength range capturing variability in 
                #clay mineralogy.
                wave_range = [1800, 2300]

            if param == "FAR":
                #Integrate over the longwave detector (2.8 microns to 3.6 microns)
                wave_range = [2800, 3900]

            if param == "CAR":
                #Integrate from 2.8 microns to 3.4 microns capturing region of water and carbonate
                wave_range = [2900, 3400]

            
            cube = mtrdr_crop_bands(img, wave_range)
            ColourSystem.cmf = mtrdr_color_matching(wave_range)
            cube = color_from_cube(cube, cs, mode=mode)
            #Export PNG file
            with rasterio.open(name+"_"+param+".png", 'w', **profile) as out:
                out.write(cube)
            
    if new_params != None:
        
        for item in new_params:
            
            if len(item) != 2:
                print("Error: Wavelength list appears to be incorrectly formatted.")
                print("New parameters should be in form [[wave1, wave2], [wave1, wave2], ...]")
                
            else:
                cube = mtrdr_crop_bands(img, item)
                ColourSystem.cmf = mtrdr_color_matching(item)
                cube = color_from_cube(cube, cs, mode=mode)
                with rasterio.open(name+"_"+str(item[0])+"_"+str(item[1])+".png", 'w', **profile) as out:
                    out.write(cube)
    
    pass


def mtrdr_to_mastcam(file, fname, narrowband=True):
    
    ##Data I/O and formatting
    with rasterio.open(file) as src:
        profile = src.profile
        cube = src.read()
    
    cube = format_mtrdr(cube)
    cube = mtrdr_crop_bands(cube, [380, 1200])
    cube = cube.transpose(1,2,0)
    cube = np.ma.masked_values(cube, 65535)
    
    #Developer note: Mastcam filter responses are stored in the following order:
    #[0] - wavelength
    #[1-4] - bayer filters (blue, green, red)
    #[4] - Left IR-bandcut
    #[5-12] - Left narrowband filters (L1-L7)
    #[12] - Right IR-bandcut
    #[13:] - Right narrowband filters (R1-R7)
    
    filter_response = np.genfromtxt("matching_functions/mastcam-response-mtrdr.txt", delimiter="\t")
    
    ##Calculate filter images filters via integration
    
    blue = calculate_luminance((filter_response[:, 1] * filter_response[:, 4]), cube)
    green = calculate_luminance((filter_response[:, 2] * filter_response[:,4]), cube)
    red = calculate_luminance((filter_response[:, 3] * filter_response[:,4]), cube)
    
    export = np.stack((red, green, blue))
    export = (export - (np.amin(export) - (0.02*np.amin(export)))) / ((np.amax(export) + (0.02*np.amax(export))))
    filter_name = "RGB"
    
    export = convert_uint16(export)
    
    #Update profile for color export
    profile.update(
        dtype = rasterio.uint16,
        count = 3,
        driver = 'PNG'
    )
    
    with rasterio.open(fname+"_"+filter_name+".png", 'w', **profile) as out:
                out.write(export)
    
    if narrowband == False:
        return
    
    #This section probably does not have the cleanest setup. Would prefer to execute this by iterating through 
    #filters, but MastCam narrowband filters are obtained by discarding two of the Bayer filters (see Bell 
    #et al. 2016 for documentation). The Bayer filters which get dropped change filter to filter, so I'm not 
    #sure I can cleanly iterate through this in a loop.
    
    #The Bayer filters are effectively transparent in the NIR and are treated as identically transparent.
    #Here I will emulate the interpolation by averaging the three Bayer filter bandpasses before applying 
    #it to the narrowband filter.
    bayer_response = np.average(filter_response[:, 1:4], axis=1)
    
    l1 = calculate_luminance((filter_response[:, 5] * filter_response[:,2]), cube)
    l2 = calculate_luminance((filter_response[:, 6] * filter_response[:,1]), cube)
    l3 = calculate_luminance((filter_response[:, 7] * filter_response[:,3]), cube)
    l4 = calculate_luminance((filter_response[:, 8] * filter_response[:,3]), cube)
    l5 = calculate_luminance((filter_response[:, 9] * bayer_response), cube)
    l6 = calculate_luminance((filter_response[:, 10]* bayer_response), cube)
    
    r1 = calculate_luminance((filter_response[:,13] * filter_response[:,2]), cube)
    r2 = calculate_luminance((filter_response[:,14] * filter_response[:,1]), cube)
    r3 = calculate_luminance((filter_response[:,15] * filter_response[:,3]), cube)
    r4 = calculate_luminance((filter_response[:, 16]* bayer_response), cube)
    r5 = calculate_luminance((filter_response[:, 17]* bayer_response), cube)
    r6 = calculate_luminance((filter_response[:, 18]* bayer_response), cube)
    
    
    filter_list = [l1, l2, l3, l4, l5, l6, r1, r2, r3, r4, r5, r6]
    filter_names = ["L1_527nm", "L2_445nm", "L3_751nm", "L4_676nm", "L5_867nm",
                   "L6_1012nm", "R1_527nm", "R2_447nm", "R3_805nm", "R4_908nm",
                   "R5_937nm", "R6_1013nm"]
    
    profile.update(
        dtype = rasterio.uint16,
        count = 1,
        driver = 'PNG'
    )
    
    for item, name in zip(filter_list, filter_names):
        item = np.expand_dims(item, 0)
        item = convert_uint16(item)
        with rasterio.open(fname+"_"+name+".png", 'w', **profile) as out:
                out.write(item)
                
    return


def mtrdr_to_hrsc(file, fname, color="IGB", lumin=False):
    
    ##Data I/O and formatting
    with rasterio.open(file) as src:
        profile = src.profile
        cube = src.read()
    
    cube = format_mtrdr(cube)
    cube = mtrdr_crop_bands(cube, [380, 1100])
    cube = cube.transpose(1,2,0)
    cube = np.ma.masked_values(cube, 65535)
    
    #Developer note: HRSC filter responses are stored in the following order:
    #[0] - MTRDR wavelength; [1] - Nadir; [2] - NIR; [3] - Red; [4] - Green; [5] - Blue;
    #[6] - Photometry; [7] - Stereo
    
    #Approximately 15% of the light entering the HRSC blue filter in the N-UV is not visible
    #to CRISM, so the filter response is a little different from reality.
    
    #Filter information retrieved from the Spanish Virtual Observatory Filter Profile Repository
    filter_response = np.genfromtxt("matching_functions/hrsc-response-mtrdr.txt", delimiter="\t")
    
    ##Calculate filter images filters via integration
    
    nad = calculate_luminance((filter_response[:, 1]), cube)
    nir = calculate_luminance((filter_response[:, 2]), cube)
    red = calculate_luminance((filter_response[:, 3]), cube)
    grn = calculate_luminance((filter_response[:, 4]), cube)
    blu = calculate_luminance((filter_response[:, 5]), cube)
    pho = calculate_luminance((filter_response[:, 6]), cube)
    ste = calculate_luminance((filter_response[:, 7]), cube)
    
    if color == "IGB":
        export = np.stack((nir, grn, blu))
        
    elif color == "IRB":
        export = np.stack((nir, red, blu))
        
    elif color == "RGB":
        export = np.stack((red, grn, blu))
        
    else:
        print("Invalid color keyword, use 'IGB', 'IRB', or 'RGB'.")
    
    export = (export - (np.amin(export) - (0.02*np.amin(export)))) / ((np.amax(export) + (0.02*np.amax(export))))    
    export = convert_uint16(export)
    
    #Update profile for color export
    profile.update(
        dtype = rasterio.uint16,
        count = 3,
        driver = 'PNG'
    )
    
    with rasterio.open(fname+"_"+color+".png", 'w', **profile) as out:
                out.write(export)
    
    if lumin == False:
        return
    
    filter_list = [nad, nir, red, grn, blu, pho, ste]
    filter_names = ["ND", "IR", "RED", "GRN", "BLU", "P1", "S1"]
    
    profile.update(
        dtype = rasterio.uint16,
        count = 1,
        driver = 'PNG'
    )
    
    for item, name in zip(filter_list, filter_names):
        item = np.expand_dims(item, 0)
        item = convert_uint16(item)
        with rasterio.open(fname+"_"+name+".png", 'w', **profile) as out:
                out.write(item)
                
    return


def mtrdr_to_hirise(file, fname, color="IRB"):
    
    ##Data I/O and formatting
    with rasterio.open(file) as src:
        profile = src.profile
        cube = src.read()
    
    cube = format_mtrdr(cube)
    cube = mtrdr_crop_bands(cube, [380, 1100])
    cube = cube.transpose(1,2,0)
    cube = np.ma.masked_values(cube, 65535)
    
    #Developer note: HRSC filter responses are stored in the following order:
    #[0] - MTRDR wavelength; [1] - NIR; [2] - Red; [3] - Blue-Green
    
    #Filter information retrieved from the Spanish Virtual Observatory Filter Profile Repository
    filter_response = np.genfromtxt("matching_functions/hirise-response-mtrdr.txt", delimiter="\t")
    
    ##Calculate filter images filters via integration
    
    nir = calculate_luminance((filter_response[:, 1]), cube)
    red = calculate_luminance((filter_response[:, 2]), cube)
    bgr = calculate_luminance((filter_response[:, 3]), cube)

    if color == "IRB":
        export = np.stack((nir, red, bgr))
        
    elif color == "RGB":
        #If RGB is requested, calculate synthetic blue filter according to HiRISE team formula
        blu = (bgr * 2) - (0.3 * red)
        export = np.stack((red, bgr, blu))
        
    else:
        print("Invalid color keyword, use 'IRB' or 'RGB'.")
    
    export = (export - (np.amin(export) - (0.02*np.amin(export)))) / ((np.amax(export) + (0.02*np.amax(export))))    
    export = convert_uint16(export)
    
    #Update profile for color export
    profile.update(
        dtype = rasterio.uint16,
        count = 3,
        driver = 'PNG'
    )
    
    with rasterio.open(fname+"_"+color+".png", 'w', **profile) as out:
                out.write(export)
                
    return

if __name__ == '__main__':
  fire.Fire()