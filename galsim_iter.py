"""
Dillon Brout
dbrout@physics.upenn.edu

TO DO LIST:
1. Implement pixelize image
1b. pixelize real image
2. Understand RA vs DEC in images
3. figure out how to iterate
4. Implement adjust_model()
5. How do I get the galaxy positions???
"""

import galsim
import galsim.des
import numpy as np
import pyfits as pf
from scipy import stats

class GalsimKernel:
    """Pixelize the input image and use as a first guess towards a simulated image. 
    Iterate to find original image to be convolved with psf to produce input image.
    
    a. Input image = real image
    
    b. Input Image --> pixelize --> GalSim --> simulated image
    c. Adjust Input image.
    d. Repeat from part b.
    """

    def __init__(self 
                 , real_img # whole image, not stamp
                 , model_img # whole image, not stamp, WITHOUT SN
                 , galpos_ra
                 , galpos_dec
                 , SN_RA_guess = 0 # arcsec from center of entire image (not stamp)
                 , SN_DEC_guess = 0 # arsec from center of entire image (not stamp)
                 , SN_flux_guess = 0.0
                 , satisfactory = .8 # process is iterated until correlation reaches this value
                 , stamp_RA = 32
                 , stamp_DEC = 32
                 ):

        self.real_img = real_img
        self.model_img = model_img
        self.SN_flux = SN_flux_guess
        self.galpos_ra = galpos_ra
        self.galpos_dec = galpos_dec
        self.stamp_RA = stamp_RA
        self.stamp_DEC = stamp_DEC

        self.satisfactory = satisfactory

        self.model_img_pix = self.pixelize( self.model_img )

        self.model = numpy.array( self.model_img_pix, dtype=float )

        # Create supernova point source (ie. very very small gaussian)
        self.sn = galsim.Gaussian( sigma = 1.e-8, flux = self.SN_flux )

        # Shift SN relative to galaxy center
        self.sn = sn.shift(galsim.PositionD( self.SN_RA_guess, self.SN_DEC_guess ))#NEED TO DOUBLE CHECK RA vs DEC
        
        #get wcs information from real data file (this doesnt change as model changes)                                                
        self.wcs = galsim.FitsWCS( self.real_img )
        # position of galaxy in original image. (pixels) (doesnt iterate) NEED TO FIGURE OUT RA VS DEC
        self.image_pos = galsim.PositionD( self.galpos_ra, self.galpos_dec )

        # Read in real image for comparison to model
        #NEED TO FIGURE OUT HOW/IF THIS NEEDS TO BE PIXELIZED
        full_real_data_image = galsim.fits.read( self.real_img )
        
        # Chop out real data stamp NEED TO DOUBLE CHECK RA VS DEC. THEN RAVEL into single array for corr calc
        self.real_data_stamp = full_rea_data_image[ galsim.BoundsI( self.galpos_ra-self.stamp_RA, 
                                                                    self.galpos_ra+self.stamp_RA,
                                                                    self.galpos_dec-self.stamp_DEC, 
                                                                    self.galpos_dec+self.stamp_DEC 
                                                                    ) ].ravel()

    """
    This will manage the iterating process
    """
    def run(self):
        correlation = 0.0
        while correlation < self.satisfactory:
            self.adjust_sn()
            self.adjust_model()
            self.kernel()
            correlation = self.compare_model_and_sim()


    """                                                                                                                                    
    the kernel gets iterated over...                                                                                       
    """
    def kernel(self):
        # Convert model to galsim image
        self.im = galsim.Image(array=self.model, scale=0.5) # scale is arcsec/pixel

        # Create interpolated image (can mess around with interp methods...)
        self.gal_model = galsim.InterpolatedImage(image=self.im, x_interpolant='linear')

        # Make adjustments to location and flux of SN
        self.adjust_sn()

        # Combine galaxy model and supernova
        self.total_gal = self.gal_model + self.sn

        # Get psf over the entire ccd
        self.psf_model = galsim.des.DES_PSFEx(psfex_file, wcs=self.wcs)                             
        
        # We just care about psf locally at the image pos
        self.psf = psf_model.getPSF(image_pos)
        
        # Convolve galaxy+sn model with psf
        self.final = galsim.Convolve( [total_gal, psf] )
        
        # create a blank stamp to be used by drawImage
        self.sim_stamp = galsim.ImageF(64, 64, wcs=wcs.local(image_pos=image_pos) )

        final.drawImage(image=stamp)

    """
    Adjusting the guess for the location and flux of the supernova
    """
    def adjust_sn(self,flux_adj=0.0,ra_adj=0.0,dec_adj=0.0):        
        self.SN_flux += flux_adj
        self.SN_RA_guess += ra_adj
        self.SN_DEC_guess += dec_adj
        # Shift SN relative to galaxy center                                                                                         
        #NEED TO FIGURE OUT PROPER WAY TO ADJUST (DOES THIS ALWAYS SHIFT FROM CENTER?)
        #NEED TO DOUBLE CHECK RA VS DEC
        self.sn = sn.shift(galsim.PositionD( SN_RA_guess, SN_DEC_guess )) # arcsec  
        

    """
    Adjusting the galaxy model pixel values. Completely empirical!
    """
    def adjust_model(self):
        return

    """
    Use Pearson Correlation to calculate r value (univariate gaussian distr)
    See Ivezic, Connolly, VanderPlas, Gray book p115
    """
    def compare_model_and_sim(self):
        corr_coeff, p_value = stats.pearsonr(self.real_data_stamp, final.ravel()) 
        #string model and sim out into long 1D arrays and correlate
        return

    def pixelize(self,img):
        pix_img = img#NEED TO IMPLEMENT
        return pix_img

if __name__=='__main__':
    image_dir = '/global/scratch2/sd/dbrout/20130829_SN-E1/g_01/'
    real_img_without_SN = 'SNp1_228717_SN-E1_tile20_g_01.fits'
    real_image_with_SN = 'SNp1_228717_SN-E1_tile20_g_01+fakeSN.fits'
    
    # Initial guess for model is real img without SN
    test = GalsimKernel( real_img_without_SN, real_img_without_SN)
    

#model_arr = [23,26,44],[12,1,44],[15,161,23]
