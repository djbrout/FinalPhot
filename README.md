# GALSIM Iterator

This is the code suite test environment to ultimately be used for
*DES Final Photometry: Supernova Fluxes*
by Dillon Brout: dbrout@physics.upenn.edu.


Update:
=========

Currently I'm looking at the first exposure (247897) of galaxy DES13E1sae found here http://dessne.cosmology.illinois.edu/SNWG/web/display/examineCand.php?Name=DES13E1sae

I have verified that the __REAL__ data stamp is correcty zooming in on the galaxy with pixel values given by the query in the link above (1416, 1409).

Will need to trim outside of image in order to avoid poor correlations between the real data and galsim output.

The .fits of this exposure is located here: 

/global/scratch2/sd/dbrout/SN_images/20131028_SN-E1/g_05/SNp1_247897_SN-E1_tile20_g_05.fits

INPUT (LEFT) OUTPUT (RIGHT):
============================

![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/update_stampsworking.png)


Pixelating and Trimming Edges:
==============================
Pixelating has been implemented using bilinear interpolation. Currently losing RA and DEC values for the pixels so I might want to keep track of that.
Galsim also blurrs the edges of the stamp when convoliving with the psf so I'm further trimming the stamps to avoid that affecting the correlation.

![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/pixelize2.png)


Speed:
======

At the moment, the code takes 0.09 seconds to run. 0.08 of that time is taken by Drawing...


TO DO LIST:
===========
* Dont lose RA and DEC when pixelating.

* figure out how to adjust without drawing and interpolating

* Implement adjust_model()

* Allow GalsimKernel to accept arrays of galpos_ra's and decs

* How do I know which images are pre-post SN?

