# GALSIM Iterator

This is the code suite test environment to ultimately be used for
*DES Final Photometry: Supernova Fluxes*
by Dillon Brout: dbrout@physics.upenn.edu.


Update:
=========

Shown is the first exposure (247897) of galaxy DES13E1sae found here http://dessne.cosmology.illinois.edu/SNWG/web/display/examineCand.php?Name=DES13E1sae

The .fits of this exposure is located here: 

/global/scratch2/sd/dbrout/SN_images/20131028_SN-E1/g_05/SNp1_247897_SN-E1_tile20_g_05.fits

GALSIM INPUT (LEFT) OUTPUT (RIGHT):
============================

![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/update_stampsworking.png)

pixelated and trimmed

![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/sim_and_pix_working.png)



MCMC Burn In Period:
====================
![alt tag](https://raw.github.com/djbrout/FinalPhot/master/readme_files/burn_in.png)


Speed:
======

At the moment, the code takes 0.09 seconds to run. 0.08 of that time is taken by Drawing...


TO DO LIST:
===========

* Zero points using mcmc on known stars

* Calculate backgrounds using sextractor methodology (image mesh with sigma cuts and mode calc)

* Add Flux(Supernova(ra,dec))_epoch to model


