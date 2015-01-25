# GALSIM Iterator

This is the code suite test environment to ultimately be used for
*DES Final Photometry: Supernova Fluxes*
by Dillon Brout: dbrout@physics.upenn.edu.


PROGRESS:
=========

Currently I'm looking at the first exposure (247897) of galaxy DES13E1sae found here http://dessne.cosmology.illinois.edu/SNWG/web/display/examineCand.php?Name=DES13E1sae

I have verified that it is correcty zooming in on the correct pixel values given by the query in the link above (1416, 1409).

In the image shown below you can see the input model (left) which was initialized to be the "real" observed image. You can also see the output simulated image given by the psf convolution used in Galsim (right).

It seems odd that the location of the galaxy has moved slightly when comparing the input model (left) to the output of the simulation (right). Maybe its an artifact of accidentally flipping the image when feeding RA pixels and DEC pixels into galsim incorrectly.

The .fits of this exposure is located here: 

/global/scratch2/sd/dbrout/SN_images/20131028_SN-E1/g_05/SNp1_247897_SN-E1_tile20_g_05.fits

![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/update2.png)


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

