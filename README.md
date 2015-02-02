# GALSIM Iterator

This is the code suite test environment to ultimately be used for
*DES Final Photometry: Supernova Fluxes*
by Dillon Brout: dbrout@physics.upenn.edu.


PROBLEM:
=========

Currently I'm looking at the first exposure (247897) of galaxy DES13E1sae found here http://dessne.cosmology.illinois.edu/SNWG/web/display/examineCand.php?Name=DES13E1sae

I have verified that the __REAL__ data stamp is correcty zooming in on the galaxy with pixel values given by the query in the link above (1416, 1409).

Will need to trim outside of image in order to avoid poor correlations between the real data and galsim output.

The .fits of this exposure is located here: 

/global/scratch2/sd/dbrout/SN_images/20131028_SN-E1/g_05/SNp1_247897_SN-E1_tile20_g_05.fits

INPUT (LEFT) OUTPUT (RIGHT):
============================

![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/update_stampsworking.png)


Speed:
======

At the moment, the code takes 15 seconds to run, but all of that time is taken up by this line of code within the kernel:

self.gal_model = galsim.InterpolatedImage( image = self.im, x_interpolant = 'linear' )

I need to figure out if its necessary to remake the model each time using this code, or if its possible to alter the pixel values of the model without regenerating the interpolated image... Not likely...


TO DO LIST:
===========
* add trimming feature to stamps

* figure out how to adjust without drawing and interpolating

* Implement adjust_model()

* Allow GalsimKernel to accept arrays of galpos_ra's and decs

* Figure out how to ravel the Galsim image

* Implement pixelize image


* How do I know which images are pre-post SN?

