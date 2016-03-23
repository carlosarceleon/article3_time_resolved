
def extract_relevant_data( case_list = [], exceptions = [], y_delta_locs = [],
                         x_2h_locs = [] , plot = False):
    """ This will extract the wall normal data at the spanwise location
    TE at a certain y density
    """

    from os                           import listdir
    from os.path                      import join,split
    from pandas                       import DataFrame, HDFStore, read_pickle
    from boundary_layer_routines      import return_bl_parameters
    from raw_data_processing_routines import decript_case_name
    from progressbar                  import ProgressBar,Percentage
    from progressbar                  import Bar,ETA,SimpleProgress
    from numpy                        import array, round, linspace
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

        available_xy_locs = hf_coords[
            ( hf_coords.x > min( x_2h_locs ) * 40. ) & \
            ( hf_coords.x < max( x_2h_locs ) * 40. ) & \
            ( hf_coords.y > min( y_delta_locs ) * d_99 ) & \
            ( hf_coords.y < max( y_delta_locs ) * d_99 )
        ][ ['x','y'] ]
              
        available_xy_locs = [tuple(x) for x in available_xy_locs.values]

        if plot:

            trailing_edge,phi,alpha,U,z = decript_case_name( f )

            if trailing_edge == 'serrated': device = 'Sr20R21'
            elif trailing_edge == 'straight': device = 'STE'
            elif trailing_edge == 'slitted': device = 'Slit20R21'

            case_name = "{0}_phi{1}_alpha{2}_U{3}_loc{4}_tr.dat".format(
                device, phi, alpha, U, z
            )

            df_av = read_pickle( 'averaged_data/' + case_name + '.p' )
            show_surface_from_df( df_av , points = available_xy_locs ,
                                plot_name = 'ReservedData/' + f + '.png'
                                )

        query   = ''
        cnt_all = 0

        cnt = 0
        time_series_hdf = HDFStore( 'ReservedData/' + f + '.hdf5' , 'w' )

        vertical_split_blocks = 10

        progress = ProgressBar(
             widgets=[
                 Bar(),' ',
                 Percentage(),' ',
                 ETA(), ' (query bunch  ',
                 SimpleProgress(),')'], 
             maxval = vertical_split_blocks
             ).start()

        # Don't try to get it all at once; split the vertical in 4 pieces
        y_ranges = linspace( 
            min( y_delta_locs ),
            max( y_delta_locs ),
            vertical_split_blocks
        ) * d_99

        xmin = min(x_2h_locs) * 40.
        xmax = max(x_2h_locs) * 40.

        for ymin, ymax in zip( y_ranges[:-1], y_ranges[1:] ):

            query = " x>={0} & x<{1} & y>={2} & y<{3} ".\
                    format( xmin, xmax, ymin, ymax )

            df_t = hdf_t.select(
                key   = 'data',
                where = [ query ],
            )

            df_t['near_x_2h']    = round( df_t.x / 40.,  4 )
            df_t['near_y_delta'] = round( df_t.y / d_99, 4 )

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

            cnt_all += 1
            cnt     += 1

            progress.update(cnt_all)

            df_t = DataFrame()


        progress.finish()
        hdf_t.close()
        time_series_hdf.close()


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
    y_delta_locs = arange( 0.1, 2.0, 0.02 )

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
