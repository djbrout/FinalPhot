import galsim
import galsim.des
import numpy

#INPUTS, PARAMS TO CONSTRAIN
model_arr = [23,26,44],[12,1,44],[15,161,23]
SN_flux = 10
SN_x = 10 # arcsec from center of entire image (not stamp)
SN_y = 10 # arsec from center of entire image (not stamp)
##############

array = numpy.array(model_arr, dtype=float)  # You build this. (model)

im = galsim.Image(array=array, scale=0.5) # arcsec/pixel (convert model to galsim image)

gal_model = galsim.InterpolatedImage(image=im, x_interpolant='linear') #c reate interpolated image

sn = galsim.Gaussian(sigma = 1.e-8, flux = SN_flux) # make supernova point source (very small gaussian... ie delta function)

sn = sn.shift(galsim.PositionD( SN_x, SN_y ))  # arcsec (shift SN relative to galaxy center)

total_gal = gal_model + sn #combine galaxy model and supernova

wcs = galsim.FitsWCS(image_file) #get wcs information from real data file (this doesnt change as model changes)
psf_model = galsim.des.DES_PSFEx(psfex_file, wcs=wcs) #this is the psf over the entire ccd

image_pos = galsim.PositionD( image_x, image_y )  # position of galaxy in original image. (pixels)
psf = psf_model.getPSF(image_pos)#get the psf just in the position of the galaxy (we only care about the position at the center of the galaxy because psf doesnt vary much over the size of a galaxy

final = galsim.Convolve( [total_gal, psf] ) # convolve galaxy+sn model with psf

sim_stamp = galsim.ImageF(64, 64, wcs=wcs.local(image_pos=image_pos) ) # create a blank stamp to be used by drawImage

final.drawImage(image=stamp)

full_real_data_image = galsim.fits.read(image_file) # read in real data (initally all we cared about was the wcs info, now we want to compare the real data with the sim)

real_data_stamp = full_rea_data_image[ galsim.BoundsI( image_x-32, image_x+32, image_y-32, image_y+32 ) ] # chop out real data stamp

# compare real_data_stamp with sim_stamp and iterate.
