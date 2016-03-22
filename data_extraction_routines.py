
def extract_relevant_data( case_list = [], exceptions = [], y_delta_locs = [],
                         x_2h_locs = [] , plot = False):
    """ This will extract the wall normal data at the spanwise location
    TE at a certain y density
    """

    from os                           import listdir
    from os.path                      import join,split
    from pandas                       import DataFrame, HDFStore, read_pickle
    from boundary_layer_routines      import return_bl_parameters
    from raw_data_processing_routines import find_nearest
    from raw_data_processing_routines import decript_case_name
    from progressbar                  import ProgressBar,Percentage
    from progressbar                  import Bar,ETA,SimpleProgress
    from numpy                        import array, round
    from data_cleaning_routines       import show_surface_from_df

    x_2h_locs    = round( array( x_2h_locs ),    2 )
    y_delta_locs = round( array( y_delta_locs ), 2 )

    # Get the available HDF5 files #############################################
    hdf5_root = '/media/carlos/6E34D2CD34D29783/' +\
                '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'

    if not len(case_list):
        hdf5_files = [f for f in listdir( hdf5_root ) \
                      if f.endswith('.hdf5') \
                      and not f in exceptions ]
    else:
        hdf5_files = [f for f in listdir( hdf5_root ) \
                      if f.endswith('.hdf5') \
                      and f in case_list ]
    # ##########################################################################

    for hf in [join( hdf5_root, f ) for f in hdf5_files]:

        f = split( hf )[1].replace('_AirfoilNormal','')\
                .replace('_Aligned.hdf5','')

        print "   Extracting data from {0}".format(f)
        print "     at the normalized streamwise locations:"
        print "     {0}".format( x_2h_locs )


        hdf_t = HDFStore( hf, 'r' )

        # Get the available coordinates ########################################
        hf_coords = hdf_t.select('data', where = [ 't = 0' ], 
                                 columns = [ 'x', 'y' ] )
        # ######################################################################

        # Turn the non-dim requested locations into physical coords ############
        requested_locations = []
        requested_normalized_locations = []
        #for x,x_norm in zip(x_2h_locs * tooth_length, x_2h_locs):
        #    for y_d in y_delta_locs:
        #        bl_params = return_bl_parameters( f , [x] )
        #        d_99 = bl_params.delta_99.values[0]
        #        #if "STE" in f:
        #        #    d_99 = 9.4
        #        y = y_d * d_99
        #        requested_locations.append( (x,y) )
        #        requested_normalized_locations.append( ( x_norm, y_d ) )

        # Get the normalization locations depending on the case ################
        if 'z00' in f and not 'STE' in f:
            x_bl_loc = 40
        elif 'z05' in f:
            x_bl_loc = 20
        elif 'z10' in f or 'STE' in f:
            x_bl_loc = 0

        bl_params = return_bl_parameters( f , [x_bl_loc] )
        d_99 = bl_params.delta_99.values[0]

        for x,x_norm in zip(x_2h_locs * tooth_length, x_2h_locs):
            for y_d in y_delta_locs:
                y = y_d * d_99
                requested_locations.append( (x,y) )
                requested_normalized_locations.append( ( x_norm, y_d ) )
        print "    Normalizing to a BL thickness of {0:.2f} mm".\
                format(d_99)
        # ######################################################################

        available_xy_locs = []
        for x,y in requested_locations:
            thresh = 0.05

            while hf_coords[ 
                ( hf_coords.x < x + thresh ) & \
                ( hf_coords.x > x - thresh ) 
            ].empty:
                thresh += 0.01
                print "   increasing thresh to {0}".format(thresh)
                if thresh > 0.4: break
            if thresh > 0.4: continue

            # Do an inverse search quadrant around the closest x to find
            # available values of y ########################################
            possible_y_locs = hf_coords[ 
                ( hf_coords.x < x + thresh ) & \
                ( hf_coords.x > x - thresh ) 
            ].y.values

            selected_y = find_nearest( y, possible_y_locs )

            # Based on the found y value closest to that requested, request
            # the closest x value ##########################################
            possible_x_locs = hf_coords[ 
                ( hf_coords.y == selected_y )
            ].x.values

            selected_x = find_nearest( x, possible_x_locs )

            # Append it all into tuples to keep them together ##############
            available_xy_locs.append( ( selected_x, selected_y ) )

        available_xy_locs = list( set( available_xy_locs ) )

        if plot:

            trailing_edge,phi,alpha,U,z = decript_case_name( f )

            if trailing_edge == 'serrated': device = 'Sr20R21'
            elif trailing_edge == 'straight': device = 'STE'
            elif trailing_edge == 'slitted': device = 'Slit20R21'

            case_name = "{0}_phi{1}_alpha{2}_U{3}_loc{4}_tr.dat".format(
                device, phi, alpha, U, z
            )

            df_av = read_pickle( 'averaged_data/' + case_name + '.p' )
            show_surface_from_df( df_av , points = available_xy_locs )


        progress = ProgressBar(
             widgets=[
                 Bar(),' ',
                 Percentage(),' ',
                 ETA(), ' (query bunch  ',
                 SimpleProgress(),')'], 
             maxval = len( available_xy_locs )
             ).start()

        query   = ''
        cnt_all = 0

        cnt = 0
        time_series_hdf = HDFStore( 'ReservedData/' + f + '.hdf5' )
        for (x,y),(xr,yr),(xn,yn) \
                in zip(available_xy_locs, requested_locations, 
                       requested_normalized_locations):

            query = query + " ( x={0:.3f} & y={1:.3f} ) ".\
                    format( x, y )

            df_t = hdf_t.select(
                key   = 'data',
                where = query,
            )

            if not df_t.shape[0]:
                print " Could not find the coordinates"
                print "   x = {0}, y = {1}".format(x,y)

            df_t['x_requested'] = round( xr, 3 )
            df_t['y_requested'] = round( yr, 3 )

            df_t['near_x_2h']    = round( xn, 2 )
            df_t['near_y_delta'] = round( yn, 2 )


            if not cnt:
                time_series_hdf.put( 'data', df_t , 
                                    data_columns = [
                                        'near_x_2h',
                                        'near_y_delta',
                                        't'
                                    ],
                                    format = 't')
            else:
                time_series_hdf.append( 'data', df_t , 
                                       data_columns = [
                                           'near_x_2h',
                                           'near_y_delta',
                                           't'
                                       ],
                               format = 't')

            ## Now get the time series ######################################
            #time_series[-1] = time_series[-1]\
            #        .append( df_t , ignore_index = True )
            ## ##############################################################

            cnt_all += 1
            cnt     += 1

            #if cnt > 20 or cnt_all == len( available_xy_locs )-1:
            #    time_series[-1] = time_series[-1].drop_duplicates()
            #    time_series.append( DataFrame() )
            #    cnt = 0

            query    = ''

            progress.update(cnt_all)

            df_t = DataFrame()


        progress.finish()
        hdf_t.close()
        time_series_hdf.close()

        #all_time_series = DataFrame()
        #for ts in time_series:
        #    all_time_series = all_time_series.append( ts, ignore_index = True )

        #all_time_series.to_pickle( 'ReservedData/' + f + '.p' )


def run_data_extraction():
    from os import listdir
    from numpy import arange

    hdf5_root = '/media/carlos/6E34D2CD34D29783/' + \
                '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'

    case_list = sorted( [
        f for f in listdir( hdf5_root ) \
        if f.endswith('.hdf5') \
        and not 'Slit' in f \
        and 'a0' in f \
        and 'p0' in f \
    ] )

    print "   Going to retrieve selected locations from the following files"
    for c in case_list:
        print "      "+c
    
    #y_delta_locs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    y_delta_locs = arange( 0.1, 1.5, 0.02 )

    for cl in case_list:
        x_2h_locs = []
        if cl.endswith('z00_tr_Aligned.hdf5') and not cl.startswith('STE'):
            #x_2h_locs = [ 0.8, 0.9, 1.0 ]
            x_2h_locs = arange( 0.76, 1.02, 0.01 ) 
        elif cl.endswith('z05_tr_Aligned.hdf5'):
            #x_2h_locs = [ 0.3, 0.4, 0.5 ]
            x_2h_locs = arange( 0.26, 0.52, 0.01 ) 
        elif cl.endswith('_AirfoilNormal.hdf5'):
            #x_2h_locs = [ 0.0, 0.1, 0.2 ]
            x_2h_locs = arange( -0.1, 0.16, 0.01 ) 

        if len( x_2h_locs ):
            df = extract_relevant_data( 
                case_list    = [cl],
                exceptions   = [],
                y_delta_locs = y_delta_locs,
                x_2h_locs    = x_2h_locs,
                plot         = True
            )

    return df

# CONSTANTS ####################################################################

tooth_length = 40

# ##############################################################################
