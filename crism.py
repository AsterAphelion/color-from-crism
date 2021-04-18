import rasterio
import numpy as np
import spectres as spec
import png

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
    cmf = np.loadtxt('mtrdr-cie-cmf.txt', usecols=(1,2,3))

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


#Defining a few internal functions to help us on our journey.
def find_band(array, value):
    """One-liner to find the index value of the nearest band to a given 
    wavelength value."""
    idx = (np.abs(array - value)).argmin()
    return idx


def mtrdr_crop_bands(image_cube, wave_list):
    """Crops the image cube to the given wavelength range."""
    mtrdr_axis = np.genfromtxt("mtrdr_axis.tab", delimiter=",")
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
    
    new_short = find_band(mtrdr_axis, wave_list[0])
    new_long = find_band(mtrdr_axis, wave_list[1])
    crop_cube = image_cube[new_short:new_long, :, :]
    
    return(crop_cube)
    
def mtrdr_color_matching(wave_list):
    """Adjusts the CIE color matching function to span the given wavelength range."""
    #Import CIE color matching function
    #Index 0 - wavelengths, Index 1 - red matching function
    #Index 2 - green matching function, Index 3 - blue matching function
    cie_matrix = np.genfromtxt("cie-cmf.txt")

    #Import tab-delimited file of wavelength axis
    mtrdr_axis = np.genfromtxt("mtrdr_axis.tab", delimiter=",")
    mtrdr_axis = mtrdr_axis[:,2]
    
    ##Fill in the gaps where bad bands are present
    #Blue gap - 380-436 nm
    add_waves = np.linspace(377.58, 429.62, num=9)
    mtrdr_axis = np.concatenate((add_waves, mtrdr_axis))
    
    #NIR bad bands 637-710 nm
    bad_band_fill = np.linspace(637.96, 703.1, num=10)
    mtrdr_axis = np.insert(mtrdr_axis, 40, bad_band_fill, axis=0)
    mtrdr_axis = np.around(mtrdr_axis, decimals=2)
    
    ##Now to figure out where on the user-specified wavelength
    #range falls on the mtrdr axis by finding indices with closest
    #value.
    new_short = find_band(mtrdr_axis, wave_list[0])
    new_long = find_band(mtrdr_axis, wave_list[1])
    
    ##Now rescale CIE color matching function to user-specified wavelength range
    cie_matrix[:,0] = (mtrdr_axis[new_long] - mtrdr_axis[new_short]) / (cie_matrix[-1,0]-cie_matrix[0,0]) * (cie_matrix[:,0]-cie_matrix[-1,0]) + mtrdr_axis[new_long]
    red = spec.spectres(mtrdr_axis[new_short:new_long], cie_matrix[:,0], cie_matrix[:,1], fill=0, verbose=False)
    green = spec.spectres(mtrdr_axis[new_short:new_long], cie_matrix[:,0], cie_matrix[:,2], fill=0, verbose=False)
    blue = spec.spectres(mtrdr_axis[new_short:new_long], cie_matrix[:,0], cie_matrix[:,3], fill=0, verbose=False)
    new_mat = np.stack([red, green, blue], axis=-1)
    return(new_mat)

def format_mtrdr(img_cube):
    """Prepare MTRDR data cube for color production by filling in bad bands."""
    axis = np.genfromtxt("mtrdr_axis.tab", delimiter=",")
    waves = axis[:,2]
    
    #First step - create some new channels mirroring the actual CRISM channel
    #spacing to fill the gap from 436 nm back to 377 nm. Doing this by averaging
    #the first six channels to reduce noise when we start subtracting a slope from
    #these channels later.
    add_waves = np.linspace(377.58, 429.62, num=9)
    add_waves = np.around(add_waves, decimals=2)

    #Assigning the new wavelength axis to a new object for the moment because we
    #still need the old one for the next steps.
    new_waves = np.concatenate((add_waves, waves))

    interp_channels = np.empty((9, img_cube.shape[1], img_cube.shape[2]))
    interp_channels[0:9,:,:] = np.average(img_cube[0:6], axis=0)

    #Second step - find blue slope. Another noise-reduction step I'm taking here is
    #to find several different slopes using different pairs of channels and averaging
    #the result.
    slope = (img_cube[2] - img_cube[6]) / (waves[6] - waves[2])
    slope2 = (img_cube[1] - img_cube[5]) / (waves[5] - waves[1])
    slope3 = (img_cube[0] - img_cube[4]) / (waves[4] - waves[0])
    slope = (slope + slope2 + slope3) / 3

    #Third step - creating a new array with the amount that needs to be added/subtracted
    #from the extrapolated channels. 
    slope = np.broadcast_to(slope, (9, slope.shape[0], slope.shape[1]))
    multiplier = (waves[0] - add_waves)
    multiplier = multiplier[:, np.newaxis, np.newaxis]
    slope = slope * multiplier

    #Fourth step - do the math and slap it in
    interp_channels = interp_channels + slope
    img_cube = np.concatenate((interp_channels, img_cube), axis=0)
    
    #Removing for memory management
    del interp_channels

    #Now we can reassign new_waves back to the original waves object since we're done with it.
    waves = new_waves

    #Repeating the previous steps, but this time filling in the VIS-NIR bad bands.
    bad_band_fill = np.linspace(637.96, 703.1, num=10)
    bad_band_fill = np.around(bad_band_fill, decimals=2)

    interp_channels = np.empty((10, img_cube.shape[1], img_cube.shape[2]))
    interp_channels[0:10,:,:] = img_cube[39]

    slope = (img_cube[40] - img_cube[39])/(waves[40] - waves[39])

    slope = np.broadcast_to(slope, (10, slope.shape[0], slope.shape[1]))
    multiplier = (bad_band_fill - waves[39]) 
    multiplier = multiplier[:, np.newaxis, np.newaxis]
    slope = slope * multiplier

    interp_channels = interp_channels + slope
    
    #Splitting up the array so we can merge everything
    #Blue side of the bad bands
    short = img_cube[0:40]
    #Red side of the bad bands
    long = img_cube[40::]

    intermed = np.concatenate((short,interp_channels), axis=0)
    
    #Removing for memory management
    del interp_channels
    
    img_cube = np.concatenate((intermed, long), axis=0)
    waves = np.insert(waves, 40, bad_band_fill, axis=0)
    
    return(img_cube)

def color_from_cube(cube, cs, mode="raw"):
    #Transpose array to put the wavelength axis last - personal preference
    cube = cube.transpose(1,2,0)

    #Reorganizing the pixels to 1-D to speed up a loop we'll be doing in a moment.
    rows = cube.shape[0]
    cols = cube.shape[1]
    pixels = rows*cols
    cube = cube.reshape(pixels, cube.shape[2])

    #Create a new cube to handle the interpolated color data
    clone_cube = np.empty((cube.shape[0], 3))

    #Convert the wavelength range to RGB values. If this can be broadcast this would
    #run much more quickly.
    for pixel in range(0, pixels):
        clone_cube[pixel] = cs.spec_to_rgb(cube[pixel])
        
    #An issue with the color conversion is that color values typically integrate well out
    #of the [0-1] normalized range of color values (typical output values for red in Mars
    #images typically seem to be in the 2-3 range). If left to the ColourSystem class, each
    #pixel would be normalized independently and we would end up with no relative color 
    #information for the red band. So I stripped that out of the original class definition
    #and we're going to normalize to the [0-1] range differently. "Raw" is a simple normalization
    #based on the highest overall value , "WB" 
    
    def quicknorm(data):
        data = (data-np.amin(data))/np.amax(data)
     
        return data
    
    #Simple normalization based on highest overall value. Provides something resembling what an
    #observer would see.
    if mode=="raw":
        clone_cube = quicknorm(clone_cube)
    
    #Normalize each channel independently - this is similar to how the VIS/FAL parameter products
    #work.
    if mode=="wb":
        for channel in range(0, clone_cube.shape[1]):
            clone_cube[:,channel] = quicknorm(clone_cube[:,channel])
        
    #The output from the previous line is color only - no luminance data at all.
    #So here, we're going to create a luminance frame using the average brightness across
    #the wavelength range of interest.
    
    #Doing a contrast stretch on the data, but making sure that we're not stripping any of the
    #lights or darks. Subtract 2% from the darkest dark and add 2% from the lightest light.
    #That way there's a little bit of cushion.
    lumin = np.average(cube, axis=1)
    lumin = (lumin - (np.amin(lumin) - (0.02*np.amin(lumin)))) / ((np.amax(lumin) + (0.02*np.amax(lumin))))
     
    #Add luminance data to cube
    clone_cube = clone_cube * lumin[:, np.newaxis]

    #Reshape array to original shape (cols multiplied by 3 to hold the RGB color values).
    cube = clone_cube.reshape(rows, cols * 3)
    
    #Convert sRGB space values to unsigned 16-bit
    cube = cube * 65535
    cube = cube.astype(np.uint16)
    
    return(cube)

def mtrdr_to_color(file, name, standard_params=True, new_params=None):

    a = rasterio.open(file,mode='r+',driver='ISIS3')
    img = a.read()
    a.close()
    cs = cs_srgb

    #Make null values = 0 so that it doesn't break when doing rgb conversion
    #Also need to convert the null pixels outside of image to 0.
    img[img < 0] = 0
    img[img >= 1] = 0

    img = format_mtrdr(img)
    
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
            png.fromarray(cube, mode="RGB;16").save(name+"_"+param+".png")
            
    if new_params != None:
        
        for item in new_params:
            
            if len(item) != 2:
                print("Error: Wavelength list appears to be incorrectly formatted.")
                print("New parameters should be in form [[wave1, wave2], [wave1, wave2], ...]")
                
            else:
                cube = mtrdr_crop_bands(img, item)
                ColourSystem.cmf = mtrdr_color_matching(item)
                cube = color_from_cube(cube, cs, mode=mode)
                
                png.fromarray(cube, mode="RGB;16").save(name+"_"+str(item[0])+"_"+str(item[1])+".png")
    
    return
