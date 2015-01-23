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

The .fits of this exposure is located here /global/scratch2/sd/dbrout/SN_images/20131028_SN-E1/g_01/SNp1_247897_SN-E1_tile20_g_01.fits

You can verify by eye that there is no galaxy centered at pixel values (1416, 1409)

![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/update1.png)


TO DO LIST:
-----------
* Check that it is correctly on top of a galaxy

* Need to convert galpos_ra and galpos_dec from sky pos to pixels
* Allow GalsimKernel to accept arrays of galpos_ra's and decs
* Double check stamp to see if it is correctly centered on galaxy

* figure out how to ravel the Galsim image
* Implement pixelize image
* pixelize real image
* Understand RA vs DEC in images
* figure out how to iterate
* Implement adjust_model()
* How do I get the galaxy positions???
* How do i know which images are pre-post SN?

