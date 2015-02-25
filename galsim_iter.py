"""
Dillon Brout
dbrout@physics.upenn.edu

Currently only outputs real and sim images for the first epoch provided...

See README.md for To Do List

"""

print 'importing'
import sys
sys.path.append("/global/u1/d/dbrout/FinalPhot/lib/") 

import galsim
import galsim.des
import numpy as np
import pyfits as pf
from scipy import stats
import os
import time
import rdcol
from scipy.ndimage.interpolation import zoom
from scipy.ndimage.filters import median_filter

# def meanclip(indata, clipsig=3.0, maxiter=5, converge_num=0.02, verbose=0):
import sigma_clip

print 'done importing'
print 'innitting'

class GalsimKernel:
    """Pixelize the input image and use as a first guess towards a simulated image. 
    Iterate to find original image to be convolved with psf to produce input image.
    
    a. Input image = real image
    
    b. Input Image --> pixelize --> GalSim --> simulated image
    c. Adjust Input image.
    d. Repeat from part b.
    """

    def __init__( self 
                 , real_images = None
                 , exposure_nums = None
                 , ccd_nums = None
                 , file_path = ''
                 , which_filters = None
                 , galpos_ras = None
                 , galpos_decs = None
                 , SN_RA_guess = None # arcsec from center of entire image (not stamp)
                 , SN_DEC_guess = None # arsec from center of entire image (not stamp)
                 , satisfactory = 39 # process is iterated until chisq reaches this value
                 , stamp_RA = 18
                 , stamp_DEC = 18
                 , psf_files = None
                 , weights_files = None
                 , outdir = None
                 , trim_edges = 1 # num pixels
                 , coarse_pixel_scale = .26 #arcsec
                 , results_tag = 'test' 
                 , run_time = 180
                 , model_img_index = 0
                 , background_mesh_pix_size = 256
                 , background_mesh_median_filter_size = 3 # A value of one means no filter is applied
                 , write_to_file_img_num = 0
                 , SN_counts_guesses = None
                 ):

        if real_images is None:
            raise AttributeError('Must provide real_images in __init__')
        if exposure_nums is None:
            raise AttributeError('Must provide exposure_nums in __init__')
        if ccd_nums is None:
            raise AttributeError('Must provide ccd_nums in __init__')
        if which_filters is None:
            raise AttributeError('Must provide which_filters in __init__')
        if galpos_ras is None:
            raise AttributeError('Must provide galpos_ras in __init__')
        if galpos_decs is None:
            raise AttributeError('Must provide galpos_decs in __init__')
        if psf_files is None:
            raise AttributeError('Must provide psf_files in __init__')
        if weights_files is None:
            raise AttributeError('Must provide weights_files in __init__')     
        if SN_counts_guesses is None:
            raise AttributeError('Must provide SN_flux_guesses in __init__') 

        print exposure_nums
        print len(exposure_nums)
        print len(real_images)
        #NEED TO REQUIRE THAT ALL IMAGE/PSF/WEIGHTS/etc... ARE THE SAME NUMBER OF FILES/length!
        oktogo = False
        if len(real_images) == len(weights_files):
            if len(real_images) == len(psf_files):
                if len(real_images) == len(galpos_ras):
                    if len(real_images) == len(galpos_decs):
                        if len(real_images) == len(which_filters):
                            if len(real_images) == len(exposure_nums):
                                if len(real_images) == len(ccd_nums):
                                    oktogo = True

        if not oktogo:
            raise AttributeError('Require that the dimensions of the following all match: \
            \n\treal_images\n\tweights_files\n\tpsf_files\n\tgalpso_ras\n\tgalpos_decs\n\twhich_filters\n\texposure_nums\n\tccd_nums')

        self.run_time = run_time
        self.exposure_nums = np.array(exposure_nums)
        self.ccd_nums = np.array(ccd_nums)
        self.which_filters = np.array(which_filters,dtype='string')
        self.background_mesh_median_filter_size = background_mesh_median_filter_size
        self.write_to_file_img_num = write_to_file_img_num

        #self.SN_fluxes = np.zeros(len(real_images))#initialize to zero
        self.SN_fluxes = SN_counts_guesses
        self.SN_RA_guesses = np.zeros(len(real_images))#initialize to zero
        self.SN_DEC_guesses = np.zeros(len(real_images))#initialize to zero

        self.galpos_ras = np.array( galpos_ras, dtype = float)
        self.galpos_decs = np.array( galpos_decs, dtype = float )
        self.stamp_ra_center = np.mean(self.galpos_ras)
        self.stamp_dec_center = np.mean(self.galpos_decs)
        self.galpos_backgrounds = np.zeros( len( galpos_ras ) )
        self.stamp_RA = float( stamp_RA )
        self.stamp_DEC = float( stamp_DEC )

        self.satisfactory = satisfactory
        self.trim_edges = trim_edges

        # setting up a dictionary for all the images and their relevant information
        # Unique identifier (exposure_num,ccd_num)
        self.image_library = {}
        for i in np.arange(len(self.exposure_nums)):
            self.image_library[(self.exposure_nums[i],self.ccd_nums[i])] = {}

        print self.image_library.keys()
        #raw_input()


        real_img_files = []
        [ real_img_files.append(os.path.join( file_path, real_image )) for real_image in real_images ]

        weights_files_long = []
        [ weights_files_long.append(os.path.join( file_path, weights_file )) for weights_file in weights_files ] 
        
        #model_img_file = os.path.join( file_path, real_images[model_img_index] )
        
        self.DES_PSFEx_files = []
        [ self.DES_PSFEx_files.append(os.path.join( file_path, psf_file )) for psf_file in psf_files ]

        self.results_npz = os.path.join( outdir, 'RESULTS_'+results_tag+'.npz' )

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        self.outdir = outdir

        self.real_fits = []
        [ self.real_fits.append(pf.open( real_img_file )) for real_img_file in real_img_files ]
        
        self.real_imgs = []
        [ self.real_imgs.append(real_fit[0].data) for real_fit in self.real_fits ]
        self.real_headers = []
        [ self.real_headers.append(real_fit[0].header) for real_fit in self.real_fits ] 
        self.weights_fits = []
        [ self.weights_fits.append(pf.open( weights_file_long )) for weights_file_long in weights_files_long ]
        self.weights = []
        [ self.weights.append(weights_fit[0].data) for weights_fit in self.weights_fits ]



        #FIND BACKGROUNDS
        self.background_mesh()


        ##SET THE MODEL IMAGE TO ONE OF THE REAL IMAGES (SET BY INPUT INDEX)
        #self.model_img = self.real_imgs[model_img_index]

        self.pixel_scale = self.real_headers[0]['PIXSCAL1']
        self.coarse_factor = self.pixel_scale/coarse_pixel_scale
        #self.pixel_scale = 0.2634
        

        # Create supernova point source (ie. very very small gaussian)
        self.sns = []
        [ self.sns.append(galsim.Gaussian( sigma = 1.e-8, flux = SN_flux )) for SN_flux in self.SN_fluxes]

        # Shift SNs relative to galaxy center
        [ self.sns[i].shift(galsim.PositionD( self.galpos_ras[i][0], self.galpos_decs[i][0] )) for i in np.arange(len(self.galpos_ras))]
        
        #get wcs information from real data file (this doesnt change as model changes)                                                
        self.wcss = []
        [ self.wcss.append(galsim.FitsWCS( real_img_file )) for real_img_file in real_img_files ]

        # Get psf over the entire ccd
        self.psf_models = []
        [ self.psf_models.append(galsim.des.DES_PSFEx( self.DES_PSFEx_files[i], wcs=self.wcss[i])) for i in np.arange(len(self.wcss)) ]

        # position of galaxy in original image. (pixels) (doesnt iterate) NEED TO FIGURE OUT RA VS DEC
        self.image_poss = []
        [ self.image_poss.append(galsim.PositionD( self.galpos_ras[i][0], self.galpos_decs[i][0] )) for i in np.arange(len(self.galpos_ras)) ]
        #self.image_pos = galsim.PositionD( self.stamp_ra_center, self.stamp_dec_center )

        # We just care about psf locally at the image pos
        self.psfs = []
        [ self.psfs.append(self.psf_models[i].getPSF( self.image_poss[i] )) for i in np.arange(len(self.psf_models))]

        # DEMO 7. Gets rid of FFT Runtime Error
        # http://stackoverflow.com/questions/24966419/fft-runtime-error-in-running-galsim?newreg=fb140d2381ff47cda008d3ade724ed59
        self.big_fft_params = galsim.GSParams( maximum_fft_size = 10240 )

        # Read in real image for comparison to model
        #NEED TO FIGURE OUT HOW/IF THIS NEEDS TO BE PIXELIZED
        full_real_data_images = []
        [full_real_data_images.append(galsim.fits.read( real_img_file )) for real_img_file in real_img_files]
        full_weights = []
        [ full_weights.append(galsim.fits.read( weights_file_long )) for weights_file_long in weights_files_long ]

        # Chop out real data stamp for each epoch
        self.real_data_stamps = []
        [ self.real_data_stamps.append(full_real_data_images[i][ galsim.BoundsI( int( self.galpos_ras[i][0]-self.stamp_RA ) 
                                                                    , int( self.galpos_ras[i][0]+self.stamp_RA )
                                                                    , int( self.galpos_decs[i][0]-self.stamp_DEC )
                                                                    , int( self.galpos_decs[i][0]+self.stamp_DEC )
                                                                    ) ]) for i in np.arange(len(full_real_data_images)) ]
        # Chop out an uncertainty stamp for each epoch
        self.weights_stamps = []
        [ self.weights_stamps.append(full_weights[i][ galsim.BoundsI( int( self.galpos_ras[i][0]-self.stamp_RA ) 
                                                                    , int( self.galpos_ras[i][0]+self.stamp_RA )
                                                                    , int( self.galpos_decs[i][0]-self.stamp_DEC )
                                                                    , int( self.galpos_decs[i][0]+self.stamp_DEC )
                                                                    ) ]) for i in np.arange(len(full_weights)) ]


        # SET THE MODEL TO THE REAL DATA SPECIFIED BY THE INDEX GIVEN
        self.model = self.real_data_stamps[model_img_index].array - self.galpos_backgrounds[0]
        self.model = np.ascontiguousarray(np.flipud(np.fliplr(self.model.T)))

        self.real_data_stamps_pixelated = []
        [ self.real_data_stamps_pixelated.append(self.pixelize( real_data_stamp.array )) for real_data_stamp in self.real_data_stamps ]
        self.weights_stamps_pixelated = []
        [ self.weights_stamps_pixelated.append(self.pixelize( weights_stamp.array )) for weights_stamp in self.weights_stamps ]
        #self.model_img_pix = self.model_img



        # TRIMMING Boundaries
        self.real_data_stamps_trimmed = []
        [ self.real_data_stamps_trimmed.append( real_data_stamp_pixelated[ self.trim_edges:-self.trim_edges
                                                                        , self.trim_edges:-self.trim_edges
                                                                     ] ) for real_data_stamp_pixelated in self.real_data_stamps_pixelated ]

        self.weights_stamps_trimmed = []
        [ self.weights_stamps_trimmed.append(weights_stamp_pixelated[ self.trim_edges:-self.trim_edges
                                                                        , self.trim_edges:-self.trim_edges
                                                                     ] ) for weights_stamp_pixelated in self.weights_stamps_pixelated ]


        self.real_data_stamps_ravel = []
        [ self.real_data_stamps_ravel.append( real_data_stamp.ravel() ) for real_data_stamp in self.real_data_stamps_trimmed ]
        self.weights_stamps_ravel = []
        [ self.weights_stamps_ravel.append( weights_stamp.ravel() ) for weights_stamp in self.weights_stamps_trimmed ]


        real_data_filename = 'test_data_out.fits'
        real_data_file_out = os.path.join( self.outdir, real_data_filename )
        os.system('rm '+real_data_file_out)
        pf.writeto(real_data_file_out, self.real_data_stamps_trimmed[self.write_to_file_img_num])
        
        #real_data_filename_beforepix = 'test_data_out_before_pix.fits'
        #real_data_file_out_beforepix = os.path.join( self.outdir, real_data_filename_beforepix )

        #os.system('rm '+real_data_file_out_beforepix )
        #pf.writeto(real_data_file_out_beforepix, self.real_data_stamps[0])
        

        #self.real_stamp = pf.open( real_data_file_out )[0].data
        #self.real_stamp_array = self.real_stamp.ravel()
        #self.weights_stamp_array = self.weights_stamp_trimmed.ravel()
        

        self.sim_filename = 'test_sim_out.fits'
        self.sim_full_filename = 'test_sim_pix_out.fits'
        self.model_file = 'test_model_out.fits'
        
        self.simoutfile = os.path.join(self.outdir,self.sim_filename)
        self.simpixout = os.path.join(self.outdir,self.sim_full_filename)
        self.model_file_out = os.path.join(self.outdir,self.model_file)

        '''
        Set iteration parameters
        '''
        self.chisq = []
        self.chisq.append(9999)
        self.model_pixels = []
        [ self.model_pixels.append([]) for i in np.nditer(self.real_data_stamps_ravel[model_img_index])]

        self.pixel_history = []
        self.sn_flux_history = []
        [ self.sn_flux_history.append([]) for i in np.arange(len(self.galpos_ras))]
        self.accepted_history = 0.5
        self.accepted_int = 0
        print 'Done Innitting'

    """
    This will manage the iterating process
    """
    def run( self ):
        self.thischisq = 999999.9
        self.t1 = time.time()
        self.counter = 0
        self.t2 = time.time()

        while self.t2-self.t1 < self.run_time:
            self.t2 = time.time()
            self.counter += 1
            self.accepted_int += 1
            
            #This is it!
            self.mcmc()
            #print counter

        self.summarize_run()

    def summarize_run( self ):
        self.t2 = time.time()
        print 'Total Time: ' + str( self.t2 - self.t1 )
        print 'Num Iterations: ' + str( self.counter )
        print 'Accepted Percentage: ' + str( self.accepted_history )
        np.savez(self.results_npz, pixel_history = self.pixel_history
                                , simulated_stamps = self.simulated_images
                                , data_stamps = self.real_data_stamps_trimmed
                                , sn_flux_history  = self.sn_flux_history
                                )
        os.system( 'rm ' + self.simpixout )
        pf.writeto( self.simpixout, self.simulated_images[self.write_to_file_img_num] )

        os.system( 'rm ' + self.model_file_out )
        pf.writeto( self.model_file_out, np.ascontiguousarray(np.flipud(np.fliplr(self.model.T))) )

    def mcmc(self):
        self.adjust_model() 
        self.adjust_sn()
            
        self.kernel()
        
        self.thischisq = self.compare_sim_and_real()

        #decide whether to accept new values
        accept_bool = self.accept()

        if accept_bool:
            #print 'accepted'
            self.accepted_history = ( self.accepted_history * self.accepted_int + 1.0 ) / ( self.accepted_int + 1 )
            self.copy_adjusted_image_to_model()
            self.update_history()
            self.chisq.append(self.thischisq)
        else:
            self.accepted_history = ( self.accepted_history * self.accepted_int ) / ( self.accepted_int + 1 )
            self.update_unaccepted_history()

    def accept(self):
        alpha = np.exp( self.chisq[-1] - self.thischisq ) / 2.0
        return_bool = False
        if alpha >= 1:
            return_bool = True
        else:
            if np.random.rand() < alpha:
                return_bool = True
        return return_bool

    def copy_adjusted_image_to_model(self):
        self.model = self.kicked_model
        self.sns = self.kicked_sns
        self.SN_fluxes = self.kicked_SN_fluxes
        return

    def update_history(self):
        #if self.thischisq < self.burn_in_chisq : #dont count burn-in period
        self.pixel_history.append( self.kicked_model )
        
        [ self.sn_flux_history[epoch].append( self.kicked_SN_fluxes[epoch] ) for epoch in np.arange(len(self.kicked_SN_fluxes)) ]
        return


    def update_unaccepted_history(self):
        #if self.thischisq < self.burn_in_chisq :#dont count burn in period
        self.pixel_history.append( self.model )
        [ self.sn_flux_history[epoch].append( self.SN_fluxes[epoch] ) for epoch in np.arange(len(self.SN_fluxes)) ]
        return

    """                                                                                                                                    
    the kernel gets iterated over...                                                                                       
    """
    def kernel( self ):

        # Convert model to galsim image
        self.im = galsim.Image( array = self.kicked_model, scale = self.pixel_scale ) # scale is arcsec/pixel

        # Create interpolated image (can mess around with interp methods...)
        self.big_fft_params = galsim.GSParams(maximum_fft_size=10240)

        self.gal_model = galsim.InterpolatedImage( image = self.im, x_interpolant = 'linear')

        self.simulated_images = []

        # Convolve the model with the psf for each epoch and draw to stamp, then pixelize --> self.simulated_image
        for epoch in np.arange(len(self.DES_PSFEx_files)):

            # Combine galaxy model and supernova
            self.total_gal = self.gal_model + self.sns[epoch]
                    

            # Convolve galaxy+sn model with psf
            self.final_big_fft = galsim.Convolve( [self.total_gal, self.psfs[epoch]], gsparams = self.big_fft_params )

            # create a blank stamp to be used by drawImage
            self.sim_stamp = galsim.ImageF( self.stamp_RA*2 + 1, self.stamp_DEC*2 + 1, 
                                        wcs = self.wcss[epoch].local(image_pos=self.image_poss[epoch]) )



            self.final_out_image = self.final_big_fft.drawImage( image = self.sim_stamp, wcs = self.wcss[epoch].local(image_pos=self.image_poss[epoch]) )
        
            #self.final_out_image.write(file_name = self.simoutfile)

            #self.model_img = pf.open(self.simoutfile)[0].data
            self.sim_img_pix = self.pixelize( self.final_out_image.array )

            self.simulated_images.append(self.sim_img_pix[ self.trim_edges:-self.trim_edges
                                                    , self.trim_edges:-self.trim_edges
                                                ])


    """
    Adjusting the guess for the location and flux of the supernova
    """
    def adjust_sn( self, stdev = 10):        
        self.kicked_sns = self.sns
        self.kicked_SN_fluxes = self.SN_fluxes
        for epoch in np.arange(len(self.sns)):
            self.kicked_SN_fluxes[epoch] += np.random.normal(scale = stdev )
            #self.SN_RA_guess += ra_adj
            #self.SN_DEC_guess += dec_adj
            self.kicked_sns[epoch] = galsim.Gaussian( sigma = 1.e-8, flux = self.SN_fluxes[epoch] )

    """
    Adjusting the galaxy model pixel values. Completely empirical!
    """
    def adjust_model( self , stdev = 1):
        self.deltas = np.random.normal(scale= stdev , size= self.model.shape )
        self.kicked_model = self.model + self.deltas
        return

    """
    Use Pearson Correlation to calculate r value (univariate gaussian distr)
    See Ivezic, Connolly, VanderPlas, Gray book p115
    """
    def pearson_corr( self ):
        sim_array = self.simulated_image.ravel()
        corr_coeff, p_value = stats.pearsonr( self.real_stamp_array, sim_array ) 
        #string model and sim out into long 1D arrays and correlate
        return p_value

    """
    Compute Chi Squared for 2 maps and an 1/variance map.
    """
    def compare_sim_and_real( self ):
        chisq = 0.0
        for epoch in np.arange(len(self.simulated_images)):
            this_sim_image_ravel = self.simulated_images[epoch].ravel()
            i = -1
            chisq += np.sum( (this_sim_image_ravel + self.galpos_backgrounds[epoch] - self.real_data_stamps_ravel[epoch])**2 * self.weights_stamps_ravel[epoch])
            #chisq += np.sum( (this_sim_image_ravel - self.real_data_stamps_ravel[epoch])**2 * self.weights_stamps_ravel[epoch])
            #for sim_pixel in this_sim_image_ravel:
            #    i += 1
            #    chisq += (sim_pixel - self.real_data_stamps_ravel[epoch][i])**2 * self.weights_stamps_ravel[epoch][i]

        print chisq
        return chisq

    def pixelize( self, img, zoomxfactor=None ):
        if zoomxfactor is None:
            zoomxfactor = self.coarse_factor
        #print img.shape
        pix_img = zoom(img,zoomxfactor,order=1)# order 1 is bilinear
        #print pix_img.shape
        return pix_img

    def summarize( self ):
        data = np.load(self.results_npz)
        pixel_history = data['pixel_history']
        sim_stamp = data['simulated_stamp']
        real_stamp = data['data_stamp']

        last_model_pixels = pixel_history[-1]
        print last_model_pixels[0,0]

    def plot_pixel_histograms( self ):
        import matplotlib
        matplotlib.use('Agg')
        import pylab as P
        data = np.load(self.results_npz)
        pixel_history = data['pixel_history']
        sim_stamps = data['simulated_stamps']
        real_data_stamps = data['data_stamps']
        sn_flux_history = data['sn_flux_history']

        print pixel_history
        print pixel_history.shape
        pixel1_vec = []
        for step in pixel_history:
            pixel1_vec.append(step[2,1])
        pixel1_vec_np = np.asarray(pixel1_vec)
        pixel1_weight = self.weights_stamps[0].array[2,1]
        pixel1_val = self.real_data_stamps[0].array[2,1]
        pixel2_vec = []
        for step in pixel_history:
            pixel2_vec.append(step[2,5])
        pixel2_vec_np = np.asarray(pixel2_vec)
        pixel2_weight = self.weights_stamps[0].array[2,5]
        pixel2_val = self.real_data_stamps[0].array[2,5]
        pixel3_vec = []
        for step in pixel_history:
            pixel3_vec.append(step[2,4])
        pixel3_vec_np = np.asarray(pixel3_vec)      
        pixel4_vec = []
        for step in pixel_history:
            pixel4_vec.append(step[5,7])
        pixel4_vec_np = np.asarray(pixel4_vec)  
        #n, bins, patches = P.hist(pixel_vec_np, 100, histtype='stepfilled')
        #P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        P.figure(1)
        P.plot(np.arange(0,len(pixel1_vec_np)),pixel1_vec_np)
        P.plot(np.arange(0,len(pixel2_vec_np)),pixel2_vec_np)
        P.plot(np.arange(0,len(pixel3_vec_np)),pixel3_vec_np)
        P.plot(np.arange(0,len(pixel4_vec_np)),pixel4_vec_np)
        out = os.path.join(self.outdir,'pixel_history.png')
        P.savefig(out)
        P.figure(5)
        print sn_flux_history
        raw_input()
        P.plot(np.arange(0,len(sn_flux_history[0])),sn_flux_history[0])
        P.plot(np.arange(0,len(sn_flux_history[1])),sn_flux_history[1])
        P.plot(np.arange(0,len(sn_flux_history[2])),sn_flux_history[2])
        out = os.path.join(self.outdir,'sn_counts_history.png')
        P.savefig(out)
        P.figure(2)
        n, bins, patches = P.hist(pixel1_vec_np[2000:], 100, histtype='stepfilled',alpha=.3)
        P.text(150, 800, 'Hist Mean: '+str(np.mean(pixel1_vec_np[2000:])) + '\n' +
                        'Pix Real Value: ' + str(pixel1_val) + '\n' +
                        'Hist Stdev: '+str(np.std(pixel1_vec_np[2000:])) + '\n' + 
                        'Pixel sigma: '+str((1/pixel1_weight)**.5),
                        style='italic',
                        bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
        out = os.path.join(self.outdir,'pixel1_histogram.png')
        P.savefig(out)
        P.figure(3)
        n, bins, patches = P.hist(pixel2_vec_np[2000:], 100, histtype='stepfilled',alpha=.3)
        P.text(150, 800, 'Hist Mean: '+str(np.mean(pixel2_vec_np[2000:])) + '\n' +
                        'Pix Real Value: ' + str(pixel2_val) + '\n' +
                        'Hist Stdev: '+str(np.std(pixel2_vec_np[2000:])) + '\n' +
                        'Pixel Sigma: '+str((1/pixel2_weight)**.5),
                        style='italic',
                        bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
        out = os.path.join(self.outdir,'pixel2_histogram.png')
        P.savefig(out)
        #n, bins, patches = P.hist(pixel3_vec_np, 100, histtype='stepfilled',alpha=.3)
        #n, bins, patches = P.hist(pixel4_vec_np, 100, histtype='stepfilled',alpha=.3)
        #P.show()

        return


    # Takes in full image, creates mesh, and following the sextractor.pdf
    # "The background estimator is a combination of sigma clipping and mode estimation,
    # Briefly, the local background histogram is clipped iteratively until convergence 
    # at 3sigma around its median; if sigma is changed by less than 20% during that process, 
    # we consider that the field is not crowded and we simply take the mean of the clipped 
    # histogram as a value for the background; otherwise we estimate the mode with:
    # Mode = 2.5 x Median - 1.5 x Mean"
    # MESH PIXEL SIZE MUST DIVIDE INTO IMAGE SIZE CLEANLY!
    def background_mesh( self, mesh_pixel_size = 256 ): 
        self.image_meshes = []
        self.image_meshes_meanRAandDECs = []

        index = -1
        for im in self.real_imgs:
            index += 1
            this_image_SN_RA = self.galpos_ras[ index ]
            this_image_SN_DEC = self.galpos_decs[ index ]
            image_size = np.shape(im)
            imshape = (image_size[0]/mesh_pixel_size - 1, image_size[1]/mesh_pixel_size - 1)
            image_mesh = np.zeros(imshape)
            mean_RAandDecs = np.zeros(imshape)
            x_step_int = -1
            for x_step in np.arange( mesh_pixel_size, image_size[0], mesh_pixel_size ):
                x_step_int += 1
                y_step_int = -1
                for y_step in np.arange( mesh_pixel_size, image_size[1], mesh_pixel_size ):
                    y_step_int += 1
                    #try: #trying becuase when you get to the end if image width is not an exact multiple of the stepsize
                    local_pixel_array = im[ x_step - mesh_pixel_size : x_step, y_step - mesh_pixel_size : y_step ].ravel()
                    start_stdev = np.std( local_pixel_array )
                    final_mean, final_stdev, clipped_local_pixel_array = sigma_clip.meanclip( local_pixel_array, clipsig = 3, maxiter = 8 )
                    if ( start_stdev - final_stdev ) / start_stdev < .2 :
                        background = np.mean( clipped_local_pixel_array )
                    else:
                        background = 2.5 * np.median( clipped_local_pixel_array ) - 1.5 * np.mean( clipped_local_pixel_array )
                    image_mesh[x_step_int,y_step_int] = background
            filtered_image_mesh = median_filter( image_mesh, size = self.background_mesh_median_filter_size )
            self.image_meshes.append(filtered_image_mesh)
            for x_step in np.arange( mesh_pixel_size, image_size[0], mesh_pixel_size ):
                for y_step in np.arange( mesh_pixel_size, image_size[1], mesh_pixel_size ):
                    if this_image_SN_RA > x_step:
                        if this_image_SN_RA < x_step + mesh_pixel_size:
                            if this_image_SN_DEC > y_step:
                                if this_image_SN_DEC < y_step + mesh_pixel_size:
                                    self.galpos_backgrounds[ index ] = background
            #print self.galpos_backgrounds
            #raw_input()
        #print self.galpos_backgrounds
        #raw_input()
        return 


def read_query( file, image_dir, image_nums):
    query = rdcol.read( file, 1, 2, '\t')
    images, exposure_nums, field_ccds = get_all_image_names( image_dir )

    exposures = np.array( query['Exposure'] )
    ccds = np.array( query['Field-CCD'], dtype= 'string' )
    filts = np.array( query['Filter'], dtype= 'string' )

    image_paths = []
    exposures_ = []
    query_ra_pix = np.array( query['x'] )
    query_dec_pix = np.array( query['y'] )

    ra_pix = []
    dec_pix = []
    filter_out = []
    exposure_nums_out = []
    ccd_nums_out = []

    query_wheres = {}

    for exposure in np.sort(np.unique( exposures )):
        #print 'exposure: ' + str(exposure)
        for ccd in np.sort(np.unique( ccds )):
                this_ccd = int(ccd.split('-')[-1])
                query_wheres[ exposure ] = [ (exposures == exposure) & (ccds == ccd) ]

                image_paths.append( images[ (exposure_nums == str(int(exposure))) & (field_ccds == this_ccd) ]) 
                
                ra_pix.append( query_ra_pix[ (exposures == exposure) & (ccds == ccd) ] )
                dec_pix.append( query_dec_pix[ (exposures == exposure) & (ccds == ccd) ] )
                filter_out.append( filts[ (exposures == exposure) & (ccds == ccd) ] )
                exposure_nums_out.append(exposure)
                ccd_nums_out.append(ccd)


    psf_files = []
    real_images = []
    weights_files = []
    filters = []
    exposure_nums = [] 
    ccd_nums = []


    [ real_images.append(str(image_paths[image_num]).strip('[').strip(']').replace("'",'')) for image_num in image_nums ]
    [ psf_files.append(i.split('.')[0]+'.psf') for i in real_images ]
    [ weights_files.append(i.split('.')[0]+'.weight.fits') for i in real_images ]
    [ exposure_nums.append(exposure_nums_out[i]) for i in image_nums ]
    [ ccd_nums.append(ccd_nums_out[i]) for i in image_nums ]


    galpos_ras = []
    galpos_decs = []
    [ galpos_ras.append(ra_pix[i]) for i in image_nums] #in pixels
    [ galpos_decs.append(dec_pix[i]) for i in image_nums]
    [ filters.append(filter_out[i]) for i in image_nums]

    return real_images, weights_files, psf_files, filters, galpos_ras, galpos_decs, exposure_nums, ccd_nums

def get_all_image_names( image_dir ):
    images = []
    exposure_nums = []
    field_ccds = []

    #print image_dir
    
    for (dir, _, files) in os.walk( image_dir ):
        for f in files:
            path = os.path.join( dir, f )
            if os.path.exists( path ):
                if path.split( '.' )[ -1 ] == 'fits':
                    if len( path.split( '+' ) ) == 1:
                        if len( path.split('.') ) == 2:
                            try:
                                field_ccds.append( path.split( '/' )[ -2 ].split( '_' )[ 1 ] )
                                exposure_nums.append( path.split( '/' )[ -1 ].split( '_' )[ 1 ] )
                                images.append(path)
                            except IndexError:
                                a = 'Image Doesnt Belong'
                            #You should only be left with files of the form
                            #/path_to_file/SNp1_228717_SN-E1_tile20_g_01.fits

    print 'Got Images'
    imgs = np.array(images,dtype='string')
    exp_nums = np.array(exposure_nums)
    fld_ccds = np.array(field_ccds,dtype='int')
    sort = np.argsort(exp_nums)
    return imgs[sort], exp_nums[sort], fld_ccds[sort]

if __name__=='__main__':
    
    image_dir = '/global/scratch2/sd/dbrout/'
    outdir = '/global/u1/d/dbrout/FinalPhot/out'
    query_file = './queries/test.txt'

    #image_nums = [0,1,2,3,4]

    #image_nums = [0,1,5,9,13,17]

    #image_nums = [0,1,9,13]
    #SN_counts_guesses = [0,1000,1000,1000]

    image_nums = [0,1,13]
    SN_counts_guesses = [0,5000,5000]

    real_images, weights_files, psf_files, filters, galpos_ras, galpos_decs, exposure_nums, ccd_nums = read_query( query_file, image_dir, image_nums )

    #print real_images
    #print filters
    #raw_input()

    # Initial guess for model is real img without SN
    test = GalsimKernel( real_images = real_images
                        , which_filters = filters
                        , exposure_nums = exposure_nums
                        , ccd_nums = ccd_nums
                        , psf_files = psf_files
                        , weights_files = weights_files
                        , outdir = outdir 
                        , galpos_ras = galpos_ras
                        , galpos_decs = galpos_decs
                        , results_tag = 'pix_1arcsec'
                        , run_time = 120
                        , write_to_file_img_num = 1
                        , SN_counts_guesses = SN_counts_guesses
                        )
    
    #test.run()
    test.plot_pixel_histograms()
    #Check backgrounds
    #write all data and sim files to fits
    #plot flux history

