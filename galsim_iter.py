"""
Dillon Brout
dbrout@physics.upenn.edu

See README.md for To Do List

"""

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

class GalsimKernel:
    """Pixelize the input image and use as a first guess towards a simulated image. 
    Iterate to find original image to be convolved with psf to produce input image.
    
    a. Input image = real image
    
    b. Input Image --> pixelize --> GalSim --> simulated image
    c. Adjust Input image.
    d. Repeat from part b.
    """

    def __init__( self 
                 , real_img # whole image, not stamp
                 , model_img # whole image, not stamp, WITHOUT SN
                 , file_path = ''
                 , galpos_ra = 100
                 , galpos_dec = 100
                 , SN_RA_guess = 0 # arcsec from center of entire image (not stamp)
                 , SN_DEC_guess = 0 # arsec from center of entire image (not stamp)
                 , SN_flux_guess = 0.0
                 , satisfactory = 2 # process is iterated until chisq reaches this value
                 , stamp_RA = 9
                 , stamp_DEC = 9
                 , psf_file = ''
                 , weights_file = ''
                 , outdir = None
                 , trim_edges = 1 # num pixels
                 , coarse_pixel_scale = .5 #arcsec
                 ):


        real_img_file = os.path.join( file_path, real_img )
        weights_file_long = os.path.join( file_path, weights_file )
        model_img_file = os.path.join( file_path, real_img )
        self.DES_PSFEx_file = os.path.join( file_path, psf_file )

        if not os.path.exists(outdir):
            os.makedirs(outdir)
        self.outdir = outdir


        self.real_fits = pf.open( real_img_file )
        self.real_img = self.real_fits[0].data
        self.real_header = self.real_fits[0].header
        self.weights_fits = pf.open( weights_file_long )
        self.weights = self.weights_fits[0].data

        self.model_img = self.real_img

        self.pixel_scale = self.real_header['PIXSCAL1']
        self.coarse_factor = self.pixel_scale/coarse_pixel_scale
        #self.pixel_scale = 0.2634

        self.SN_flux = SN_flux_guess
        self.SN_RA_guess = SN_RA_guess
        self.SN_DEC_guess = SN_DEC_guess

        self.galpos_ra = float( galpos_ra )
        self.galpos_dec = float( galpos_dec )
        self.stamp_RA = float( stamp_RA )
        self.stamp_DEC = float( stamp_DEC )

        self.satisfactory = satisfactory
        self.trim_edges = trim_edges

        # Create supernova point source (ie. very very small gaussian)
        self.sn = galsim.Gaussian( sigma = 1.e-8, flux = self.SN_flux )

        # Shift SN relative to galaxy center
        self.sn = self.sn.shift(galsim.PositionD( self.SN_RA_guess, self.SN_DEC_guess ))#NEED TO DOUBLE CHECK RA vs DEC
        
        #get wcs information from real data file (this doesnt change as model changes)                                                
        self.wcs = galsim.FitsWCS( real_img_file )

        # Get psf over the entire ccd
        self.psf_model = galsim.des.DES_PSFEx( self.DES_PSFEx_file, wcs=self.wcs)

        # position of galaxy in original image. (pixels) (doesnt iterate) NEED TO FIGURE OUT RA VS DEC
        self.image_pos = galsim.PositionD( self.galpos_ra, self.galpos_dec )

        # We just care about psf locally at the image pos
        self.psf = self.psf_model.getPSF( self.image_pos )

        # DEMO 7. Gets rid of FFT Runtime Error
        # http://stackoverflow.com/questions/24966419/fft-runtime-error-in-running-galsim?newreg=fb140d2381ff47cda008d3ade724ed59
        self.big_fft_params = galsim.GSParams( maximum_fft_size = 10240 )

        # Read in real image for comparison to model
        #NEED TO FIGURE OUT HOW/IF THIS NEEDS TO BE PIXELIZED
        full_real_data_image = galsim.fits.read( real_img_file )
        full_weights = galsim.fits.read( weights_file_long )
        #full_real_data_image = pf.open( real_img_file )[0].data

        # Chop out real data stamp NEED TO DOUBLE CHECK RA VS DEC.
        self.real_data_stamp = full_real_data_image[ galsim.BoundsI( int( self.galpos_ra-self.stamp_RA ) 
                                                                    , int( self.galpos_ra+self.stamp_RA )
                                                                    , int( self.galpos_dec-self.stamp_DEC )
                                                                    , int( self.galpos_dec+self.stamp_DEC )
                                                                    ) ]
        self.weights_stamp = full_weights[ galsim.BoundsI( int( self.galpos_ra-self.stamp_RA ) 
                                                                    , int( self.galpos_ra+self.stamp_RA )
                                                                    , int( self.galpos_dec-self.stamp_DEC )
                                                                    , int( self.galpos_dec+self.stamp_DEC )
                                                                    ) ]

        start_model_filename = 'start_model.fits'
        start_model_out = os.path.join(self.outdir,start_model_filename)
        self.real_data_stamp.write(start_model_out)

        weights_stamp_filename = 'weights_stamp_untrimmed_unpix.fits'
        weights_stamp_out = os.path.join(self.outdir,weights_stamp_filename)
        self.weights_stamp.write( weights_stamp_out )

        self.model_img = pf.open(start_model_out)[0].data
        self.real_data_stamp_tobepixelated = pf.open(start_model_out)[0].data
        self.weights_stamp_tobepixelated = pf.open(weights_stamp_out)[0].data

        self.model_img_pix = self.pixelize( self.model_img )
        self.real_data_stamp_pixelated = self.pixelize( self.real_data_stamp_tobepixelated )
        self.weights_stamp_pixelated = self.pixelize( self.weights_stamp_tobepixelated )
        self.model_img_pix = self.model_img

        self.model = np.array( self.model_img_pix, dtype=float )
        self.model = np.ascontiguousarray(np.flipud(np.fliplr(self.model.T)))


        self.real_data_stamp_trimmed = self.real_data_stamp_pixelated[ self.trim_edges:-self.trim_edges
                                                                        , self.trim_edges:-self.trim_edges
                                                                     ]

        self.weights_stamp_trimmed = self.weights_stamp_pixelated[ self.trim_edges:-self.trim_edges
                                                                        , self.trim_edges:-self.trim_edges
                                                                     ]
        #self.real_data_stamp_trimmed = full_real_data_image[ galsim.BoundsI( int( self.galpos_ra-self.stamp_RA+self.trim_edges ) 
        #                                                            ,int( self.galpos_ra+self.stamp_RA-self.trim_edges )
        #                                                            ,int( self.galpos_dec-self.stamp_DEC+self.trim_edges )
        #                                                            ,int( self.galpos_dec+self.stamp_DEC-self.trim_edges )
        #                                                            ) ]

        real_data_filename = 'test_data_out.fits'
        real_data_file_out = os.path.join( self.outdir, real_data_filename )
        os.system('rm '+real_data_file_out)
        pf.writeto(real_data_file_out, self.real_data_stamp_trimmed)
        
        real_data_filename_beforepix = 'test_data_out_before_pix.fits'
        real_data_file_out_beforepix = os.path.join( self.outdir, real_data_filename_beforepix )

        os.system('rm '+real_data_file_out_beforepix )
        pf.writeto(real_data_file_out_beforepix,self.real_data_stamp_tobepixelated)

        self.real_stamp = pf.open( real_data_file_out )[0].data
        self.real_stamp_array = self.real_stamp.ravel()
        self.weights_stamp_array = self.weights_stamp_trimmed.ravel()

        print 'Done Innitting'

    """
    This will manage the iterating process
    """
    def run( self ):
        correlation = 10000
        while correlation > self.satisfactory:
            print 'Press Enter to continue'
            raw_input()
            t1 = time.time()
            self.adjust_sn()
            print 'Done adjusting SN'
            self.adjust_model() 
            print 'Done adjusting model'
            self.kernel()
            print 'Executed Kernel'
            correlation = self.compare_model_and_sim()
            print 'Correlated ' + str( correlation )
            t2 = time.time()
            print t2-t1

    """                                                                                                                                    
    the kernel gets iterated over...                                                                                       
    """
    def kernel( self ):
        t1 = time.time()
        print 'creating galsim image'
        # Convert model to galsim image
        self.im = galsim.Image( array = self.model, scale = self.pixel_scale ) # scale is arcsec/pixel

        t2 = time.time()
        print t2-t1
        print 'creating gal_model'
        # Create interpolated image (can mess around with interp methods...)
        self.big_fft_params = galsim.GSParams(maximum_fft_size=10240)

        self.gal_model = galsim.InterpolatedImage( image = self.im, x_interpolant = 'linear')


        t3 = time.time()
        print t3-t2
        print 'creating total_gal'
        # Combine galaxy model and supernova
        self.total_gal = self.gal_model + self.sn

        t4 = time.time()
        print t4-t3
        print 'convolving'
        # Convolve galaxy+sn model with psf
        self.final_big_fft = galsim.Convolve( [self.total_gal, self.psf], gsparams = self.big_fft_params )

        t5 = time.time()
        print t5-t4
        print 'stamping'
        # create a blank stamp to be used by drawImage
        self.sim_stamp = galsim.ImageF( self.stamp_RA*2 + 1, self.stamp_DEC*2 + 1, 
                                        wcs = self.wcs.local(image_pos=self.image_pos) )
        t6 = time.time()
        print t6-t5
        print 'drawing'

        self.sim_filename = 'test_sim_out.fits'
        self.sim_full_filename = 'test_sim_big_out.fits'
        self.simoutfile = os.path.join(self.outdir,self.sim_filename)


        self.final_out_image = self.final_big_fft.drawImage( image = self.sim_stamp, wcs = self.wcs.local(image_pos=self.image_pos) )
        
        

        #self.final_out_image_trimmed = self.final_out_image[galsim.BoundsI( int( self.trim_edges ) + 1
        #                                                            , int(self.stamp_RA*2+1) - int( self.trim_edges )
        #                                                            , int( self.trim_edges ) + 1
        #                                                            , int(self.stamp_DEC*2+1) - int( self.trim_edges )
        #                                                            )]

        self.final_out_image.write(file_name = self.simoutfile)

        self.model_img = pf.open(self.simoutfile)[0].data
        self.model_img_pix = self.pixelize( self.model_img )

        self.simulated_image = self.model_img_pix[ self.trim_edges:-self.trim_edges
                                                    , self.trim_edges:-self.trim_edges
                                                ]
        t7 = time.time()
        print t7-t6

    """
    Adjusting the guess for the location and flux of the supernova
    """
    def adjust_sn( self, flux_adj = 0.0, ra_adj = 0.0, dec_adj = 0.0 ):        
        self.SN_flux += flux_adj
        self.SN_RA_guess += ra_adj
        self.SN_DEC_guess += dec_adj
        # Shift SN relative to galaxy center                                                                                         
        #NEED TO FIGURE OUT PROPER WAY TO ADJUST (DOES THIS ALWAYS SHIFT FROM CENTER?)
        #NEED TO DOUBLE CHECK RA VS DEC
        if flux_adj == 0.0:
            pass
        else:
            self.sn = galsim.Gaussian( sigma = 1.e-8, flux = self.SN_flux )
        t2 = time.time()

        if (ra_adj == 0.0) & (dec_adj == 0.0):
            pass
        else:
            self.sn = self.sn.shift( galsim.PositionD( self.SN_RA_guess, self.SN_DEC_guess )) # arcsec  
        

    """
    Adjusting the galaxy model pixel values. Completely empirical!
    """
    def adjust_model( self ):
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
    def compare_model_and_sim( self ):
        sim_image_ravel = self.simulated_image.ravel()
        chisq = 0.0
        i = -1
        for sim_pixel in sim_image_ravel:
            i += 1
            chisq += (sim_pixel - self.real_stamp_array[i])**2 * self.weights_stamp_array[i]

        return chisq

    def pixelize( self, img, zoomxfactor=None ):
        if zoomxfactor is None:
            zoomxfactor = self.coarse_factor
        print img.shape
        pix_img = zoom(img,zoomxfactor,order=1)# order 1 is bilinear
        print pix_img.shape
        return pix_img

def read_query( file, image_dir ):
    query = rdcol.read( file, 1, 2, '\t')
    images, exposure_nums, field_ccds = get_all_image_names( image_dir )
    #print np.array( query['Exposure'] )
    #raw_input()
    #print exposure_nums
    #raw_input()
    exposures = np.array( query['Exposure'] )
    ccds = np.array( query['Field-CCD'], dtype= 'string' )


    image_paths = []
    exposures_ = []
    
    query_wheres = {}

    #print exposures
    #print 'hehe'
    #raw_input()
    for exposure in np.unique( exposures ):
        for ccd in np.unique( ccds ):
            this_ccd = int(ccd.split('-')[-1])
            #print 'exposure '+str(int(exposure))
            #print 'ccd ' + str(this_ccd)
            #print 'image path'
            #print images

            #print (exposure_nums == str(int(exposure))) & (field_ccds == this_ccd)
            #print images[10]
            #raw_input()
            #print images[(exposure_nums == str(int(exposure))) & (field_ccds == this_ccd)]
            #print ccds[(exposures == exposure) & (ccds == ccd)]
            query_wheres[ exposure ] = [ (exposures == exposure) & (ccds == ccd) ]

            #print exposure_nums == int(exposure)
            #print exposure_nums[exposure_nums == int(exposure)]
            #print field_ccds[field_ccds == this_ccd]
            image_paths.append( images[ (exposure_nums == str(int(exposure))) & (field_ccds == this_ccd) ]) 

            #print  images[ exposure_nums.index( str(int(exposure) ) ) ]
        
            #exposures_.append( exposure )
            #match up with a .fits file
    
    print 'Made Image Arrays'
    return query, query_wheres, image_paths, exposures


def get_all_image_names( image_dir ):
    images = []
    exposure_nums = []
    field_ccds = []

    print image_dir
    
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
    return np.array(images,dtype='string'), np.array(exposure_nums), np.array(field_ccds,dtype='int')

if __name__=='__main__':
    image_dir = '/global/scratch2/sd/dbrout/'
    #real_img_without_SN = 'SNp1_228717_SN-E1_tile20_g_01.fits'
    #real_image_with_SN = 'SNp1_228717_SN-E1_tile20_g_01+fakeSN.fits'
    #psf_file = 'SNp1_228717_SN-E1_tile20_g_01.psf'
    #psf_file_full = os.path.join(image_dir, psf_file)
    outdir = '/global/u1/d/dbrout/FinalPhot/out'
    query_file = './queries/test.txt'

    print 'Started Reading'
    query, query_wheres, image_paths, exposures = read_query( query_file, image_dir )


    #Start by online looking at one image, one exposure
    #print image_paths
    #raw_input()
    #print exposures

    image_num = 0

    real_img_without_SN = str(image_paths[image_num]).strip('[').strip(']').replace("'",'')
    psf_file = real_img_without_SN.split('.')[0]+'.psf'
    weights_file = real_img_without_SN.split('.')[0]+'.weight.fits'
    this_exposure_and_ccd = query_wheres[exposures[image_num]]
    print real_img_without_SN
    print psf_file
    #Need to double check x and y are correct columns
    galpos_ra = np.array(query['x'])[this_exposure_and_ccd] #in pixels
    galpos_dec = np.array(query['y'])[this_exposure_and_ccd] #in pixels

    #galpos_ra = 33
    #galpos_dec = 33

    print galpos_ra
    print galpos_dec
    print real_img_without_SN
    print psf_file
    #raw_input()

    # Initial guess for model is real img without SN
    test = GalsimKernel( real_img_without_SN, real_img_without_SN
                                            , psf_file = psf_file
                                            , weights_file = weights_file
                                            , outdir = outdir 
                                            , galpos_ra = galpos_ra
                                            , galpos_dec = galpos_dec
                                            )
    
    test.run()

#model_arr = [23,26,44],[12,1,44],[15,161,23]
