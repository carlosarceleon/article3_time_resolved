
def extract_relevant_data( case_list = [], exceptions = [], y_delta_locs = [],
                         x_2h_locs = [] ):
    """ This will extract the wall normal data at the spanwise location
    TE at a certain y density
    """

    from os                           import listdir
    from os.path                      import join,split
    from pandas                       import DataFrame, HDFStore
    from boundary_layer_routines      import return_bl_parameters
    from raw_data_processing_routines import find_nearest
    from progressbar                  import ProgressBar,Percentage
    from progressbar                  import Bar,ETA,SimpleProgress
    from numpy                        import array

    x_2h_locs    = array( x_2h_locs )
    y_delta_locs = array( y_delta_locs )

    # Get the available HDF5 files #############################################
    hdf5_root = '/media/carlos/6E34D2CD34D29783/' +\
                '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'

    if not len(case_list):
        hdf5_files = [f for f in listdir( hdf5_root ) \
                      if f.endswith('_Aligned.hdf5') \
                      and not f in exceptions ]
    else:
        hdf5_files = [f for f in listdir( hdf5_root ) \
                      if f.endswith('_Aligned.hdf5') \
                      and f in case_list ]
    # ##########################################################################

    for hf in [join( hdf5_root, f ) for f in hdf5_files]:

        f = split( hf )[1].replace('_Aligned.hdf5','')

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
        for x in x_2h_locs * tooth_length:
            for y_d in y_delta_locs:
                bl_params = return_bl_parameters( f , [x] )
                y = y_d * bl_params.delta_99.values[0]
                requested_locations.append( (x,y) )
        # ######################################################################

        available_xy_locs = []
        for x,y in requested_locations:
            thresh = 0.5

            # Do an inverse search quadrant around the closest x to find
            # available values of y ############################################
            possible_y_locs = hf_coords[ 
                ( hf_coords.x < x + thresh ) & \
                ( hf_coords.x > x - thresh ) 
            ].y.values

            selected_y = find_nearest( y, possible_y_locs )

            # Based on the found y value closest to that requested, request
            # the closest x value ##############################################
            possible_x_locs = hf_coords[ 
                ( hf_coords.y == selected_y )
            ].x.values

            selected_x = find_nearest( x, possible_x_locs )

            # Append it all into tuples to keep them together ##################
            available_xy_locs.append( ( selected_x, selected_y ) )


        progress = ProgressBar(
             widgets=[
                 Bar(),' ',
                 Percentage(),' ',
                 ETA(), ' (query bunch  ',
                 SimpleProgress(),')'], 
             maxval = len( available_xy_locs )
             ).start()

        time_series = DataFrame()

        query   = ''
        cnt_all = 0

        for (x,y),(xr,yr) in zip(available_xy_locs, requested_locations):

            query = query + " ( x={0:.3f} & y={1:.3f} ) ".\
                    format( x, y )

            df_t = hdf_t.select(
                key   = 'data',
                where = query,
            )

            if not df_t.shape[0]:
                print " Could not find the coordinates"
                print "   x = {0}, y = {1}".format(x,y)

            df_t['x_requested'] = xr
            df_t['y_requested'] = yr

            # Now get the time series ######################################
            time_series = time_series.append( df_t , ignore_index = True )
            # ##############################################################

            query    = ''

            cnt_all += 1
            progress.update(cnt_all)


        progress.finish()
        hdf_t.close()

        time_series.to_pickle( 'ReservedData/' + f + '.p' )


def run_data_extraction():
    from os import listdir

    hdf5_root = '/media/carlos/6E34D2CD34D29783/' +\
                '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'

    case_list = sorted( [f for f in listdir( hdf5_root ) \
                  if f.endswith('_Aligned.hdf5') \
                  and not 'Slit' in f] )
    
    y_delta_locs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

    for cl in case_list:
        x_2h_locs = []
        if cl.endswith('z00_tr_Aligned.hdf5') and not cl.startswith('STE'):
            x_2h_locs = [ 0.8, 0.9, 1.0 ]
        elif cl.endswith('z05_tr_Aligned.hdf5'):
            x_2h_locs = [ 0.3, 0.4, 0.5 ]
        elif cl.endswith('z10_tr_Aligned.hdf5') or cl.startswith('STE'):
            x_2h_locs = [ 0.0, 0.1, 0.2 ]

        if len( x_2h_locs ):
            df = extract_relevant_data( 
                case_list    = [cl],
                exceptions   = [],
                y_delta_locs = y_delta_locs,
                x_2h_locs    = x_2h_locs,
            )

    return df

# CONSTANTS ####################################################################

tooth_length = 40

# ##############################################################################
