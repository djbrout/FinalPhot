import numpy as np
import copy


class run():
	
	def __init__(self
				, model = None
				, stdev = None
				, data = None
				, psfs = None
				, weights = None
				, image_size = None
				, run_time = 100
				):

		if model is None:
            raise AttributeError('Must provide model array!')
        if stdev is None:
            raise AttributeError('Must provide stdev for each model parameter!')
        if data is None:
        	raise AttributeError('Must provide real data for comparison!')
        if psfs is None:
            raise AttributeError('Must provide psfs for each epoch!')
        if weights is None:
        	raise AttributeError('Must provide weights for each epoch!')
        if image_size is None:
        	raise AttributeError('Must provide the size of the image!')

        oktogo = False
        if len( model[ image_size**2+1: ] ) == len( data ):
        	if len( stdev[ image_size**2+1: ] ) == len( data ):
            	if len( data ) == len( psfs ):
            		if len( data ) == len( weights ):
                		oktogo = True
        if not oktogo:
            raise AttributeError('Require that the dimensions of the following all match: \
            \n\tmodel fluxes\n\tdata\n\tpsfs')


        self.model = model
        self.stdev = stdev
        self.data = data
        self.psfs = psfs
        self.run_time = run_time

        self.history = []
        self.sims = []
        [ self.sims.append( np.zeros( image_size, image_size ) ) for epoch in np.arange( np.len( self.model[ image_size**2 + 1 : ] ) ) ]

        self.run()


	def run( self ):
        self.lastchisq = 9999999.9
        self.t1 = time.time()
        self.counter = 0
        self.t2 = time.time()

        while self.t2-self.t1 < self.run_time:
            self.t2 = time.time()
            self.counter += 1
            self.accepted_int += 1
            self.mcmc()

        self.summarize_run()

	def summarize_run( self ):
        self.t2 = time.time()
        print 'Total Time: ' + str( self.t2 - self.t1 )
        print 'Num Iterations: ' + str( self.counter )
        print 'Accepted Percentage: ' + str( self.accepted_history )
        #np.savez(self.results_npz, pixel_history = self.pixel_history
        #                        , simulated_stamps = self.simulated_images
        #                        , data_stamps = self.real_data_stamps_trimmed
        #                        , sn_flux_history  = self.sn_flux_history
        #                        )


	def mcmc( self ):
        
        self.adjust_model()
        [ self.sims[ epoch ] = dans_convolve( self.kicked_model, self.psfs) ]
        self.thischisq = self.compare_sim_and_real()

        #decide whether to accept new values
        accept_bool = self.accept(self.lastchisq,self.thischisq)


        if accept_bool:
            print 'accepted'
            self.lastchisq = self.thischisq
            self.accepted_history = ( self.accepted_history * self.accepted_int + 1.0 ) / ( self.accepted_int + 1 )
            self.copy_adjusted_image_to_model()
            self.update_history()
            self.chisq.append(self.thischisq)
        else:
            self.accepted_history = ( self.accepted_history * self.accepted_int ) / ( self.accepted_int + 1 )
            self.update_unaccepted_history()

    def adjust_model( self ):
    	deltas = []
        [ deltas.append( np.random.normal( scale= self.stdev[i] ) ) for i in np.arange( len( self.model ) ) ]
        self.kicked_model = self.model + deltas
        return

    def compare_sim_and_real( self ):
        chisq = 0.0
        for epoch in np.arange(len(self.simulated_images)):
            this_sim_image_ravel = self.simulated_images[epoch].ravel()
            i = -1
            chisq += np.sum( (this_sim_image_ravel + self.galpos_backgrounds[epoch] - self.real_data_stamps_ravel[epoch])**2 * self.weights_stamps_ravel[epoch])

        return chisq

    def accept( self, last_chisq, this_chisq ):
        alpha = np.exp( last_chisq - this_chisq ) / 2.0
        return_bool = False
        if alpha >= 1:
            return_bool = True
        else:
            if np.random.rand() < alpha:
                return_bool = True
        return return_bool

    def copy_adjusted_image_to_model(self):
        self.model = copy.copy( self.kicked_model )
        return

    def update_history(self):
        self.history.append( self.kicked_model )
        return

    def update_unaccepted_history(self):
        self.history.append( self.model )
        return

    def get_model( self ):
    	return self.model

    def get_model_history( self ):
    	return self.history
