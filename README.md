# GALSIM Iterator

This is the code suite test environment to ultimately be used for
*DES Final Photometry: Supernova Fluxes*
by Dillon Brout: dbrout@physics.upenn.edu.

Images are located in /global/scratch2/sd/dbrout/SN_images/

PROGRESS:
=========

Currently I'm looking at the first exposure (247897) of galaxy DES13E1sae found here http://dessne.cosmology.illinois.edu/SNWG/web/display/examineCand.php?Name=DES13E1sae

After looking at the corresponding pixel values and stamping out a 32x32 image, I'm finding that there is no galaxy in this area. See image below.

I have verified that it is correcty zooming in on the correct pixel values given by the query in the link above (1416, 1409).

The .fits of this exposure is located here: 

/global/scratch2/sd/dbrout/SN_images/20131028_SN-E1/g_01/SNp1_247897_SN-E1_tile20_g_01.fits

You can verify by eye that there is no galaxy centered at pixel values (1416, 1409)

![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/update1.png)


Speed:
======

At the moment, the code takes 15 seconds to run, but all of that time is taken up by this line of code within the kernel:

self.gal_model = galsim.InterpolatedImage( image = self.im, x_interpolant = 'linear' )

I need to figure out if its necessary to remake the model each time using this code, or if its possible to alter the pixel values of the model without regenerating the interpolated image... Not likely...


TO DO LIST:
===========
* Check that it is correctly on top of a galaxy

* Allow GalsimKernel to accept arrays of galpos_ra's and decs

* Double check stamp to see if it is correctly centered on galaxy

* Figure out how to ravel the Galsim image

* Implement pixelize image

* Implement adjust_model()

* How do I know which images are pre-post SN?

