def decript_case_name(case_name):
    from re import findall

    try:
        trailing_edge = findall('[litSrTE201R]+', case_name)[0]
    except IndexError:
        print case_name
    if trailing_edge == "STE":     trailing_edge = 'straight'
    if trailing_edge == "Sr20R21": trailing_edge = 'serrated'
    if trailing_edge == "Slit20R21": trailing_edge = 'serrated'

    phi = 0 ; alpha = 0 ; U = 0 ; z = 0

    try:
        phi   = findall('phi[0-9]',case_name)[0].replace('phi','')
        alpha = findall('alpha[0-9][0-9]?',case_name)[0].replace('alpha','')
        z     = findall('loc[0-9][0-9]',case_name)[0].replace('loc','')

        if   z == '00' : z = '0'
        elif z == '05' : z = '0.25'
        elif z == '10' : z = '0.5'

    except:
        # In case it's the DaVis naming convention ############################
        if not phi:
            phi   = findall('p-?[0-9]',case_name)[0].replace('p','')
        if not alpha:
            alpha = findall('a-?[0-9][0-9]?',case_name)[0].replace('a','')
        if not z:
            z     = findall('z[0-9][0-9]',case_name)[0].replace('z','')
        # #####################################################################

    U = findall('U[0-9][0-9]',case_name)[0].replace('U','')

    return trailing_edge,phi,alpha,U,z

def get_angle_between_points(point1, point2):
    from math import atan, pi, degrees

    delta_y = point2[1] - point1[1]
    delta_x = point2[0] - point1[0]
    if delta_x == 0:
        return pi

    return degrees(atan( delta_y / delta_x ))


def read_tecplot_file(
    tecplot_folder,
    tecplot_time_step_file,
    time_step             = 0,
):

    """Reads in a tecplot file, given, and returns a pandas data frame

    Important!
    This data frame that is returned is already CORRECTED and turned to
    the standard coordinate system

    Input: address to tecplot_file

    Output: a data frame containing the 0.25 mm resolved structured grid data

    """
    import pandas as pd 
    from re import findall
    from os.path import split,join

    # Get available variables
    tecplot_file = join( tecplot_folder, tecplot_time_step_file )
    f = open(tecplot_file,'ro')

    # Read the header and get variable and title information from it ###########
    var_flag = False
    dev_flag = False
    for line in f:
        string = findall("^VARIABLES[ _A-Za-z0-9,\"=]+",line)
        if string:
            variables = [v.replace(' ','_').replace("\"","") \
                         for v in string[0].replace("VARIABLES = ",'')\
                         .split(", ")]
            variables = [v for v in variables if len(v)]
            var_flag = True
        string = findall("^TITLE = [ -_A-Za-z0-9,\"=]+",line)
        if string:
            dev_flag = True
        if var_flag and dev_flag:
            break
    f.close()
    ############################################################################

    # Put the data into a data frame ###########################################
    df = pd.read_table(
            tecplot_file,
            skiprows  = 4,
            names     = variables,
            sep       = '[ \t]+',
            index_col = False,
            engine = 'python'
            )
    ############################################################################
    df = df.drop('z',1)
    df = rename_df_columns_from_DaVis_to_standard(df)

    trailing_edge,phi,alpha,U,z = decript_case_name(split(tecplot_folder)[-1])

    if trailing_edge == 'serrated': device = 'Sr20R21'
    elif trailing_edge == 'straight': device = 'STE'
    elif trailing_edge == 'slitted': device = 'Slit20R21'

    df['case_name'] = "{0}_phi{1}_alpha{2}_U{3}_loc{4}_tr.dat".format(
        device, phi, alpha, U, z
    )

    df['tecplot_folder'] = split(tecplot_folder)[-1]

    df[ 't' ] = time_step

    return df

def correct_df_translation_rotation( df ):
    from corrections import correction_dict
    from Masks import Masks as masks
    from numpy import array, round
    
    try:
        mask = array(masks[df.case_name.unique()[0]])
    except KeyError:
        mask = array(masks[df.case_name.unique()[0].replace('_tr','')])

    if df.tecplot_folder.unique()[0] in correction_dict.keys():
        x_corr, y_corr, angle_corr = correction_dict[
        df.tecplot_folder.unique()[0]
        ]
    else:
        print "    Didn't find a correction term for {0}".format( 
            df.tecplot_folder 
        )
        print "    among"
        print "    {0}".format( correction_dict.keys() )
        x_corr, y_corr, angle_corr = (0, 0, 0)

    df['x_orig'] = round( df.x.values, 3 )
    df['y_orig'] = round( df.y.values, 3 )

    df.x = df.x - mask[1,0] + y_corr
    df.y = df.y - mask[1,1] - x_corr

    mask_rotation = get_angle_between_points(
        mask[1], mask[2]
    )

    df = rotate_df( df, mask_rotation + angle_corr)

    df = df.dropna()

    df.x = round(df.x , 3)
    df.y = round(df.y , 3)

    return df


def rotate_df(df,degrees = 0):
    from math import radians
    from numpy import sin, cos
    import pandas as pd

    angle = radians(degrees)

    x = df['x'] 
    y = df['y'] 

    df_rotated = pd.DataFrame()

    df_rotated['x'] =  x*cos(angle) + y*sin(angle) 
    df_rotated['y'] = -x*sin(angle) + y*cos(angle) 

    # Velocity components
    df_rotated['u'] =  df['u']*cos(angle) \
            + df['v']*sin(angle)
    df_rotated['v'] = -df['u']*sin(angle) \
            + df['v']*cos(angle)
    df_rotated['w'] = df['w']

    df_rotated['t'] = df.t
    df_rotated['x_orig'] = df.x_orig
    df_rotated['y_orig'] = df.y_orig


    for col in df.select_dtypes(include=['object']).columns:
        df_rotated[col] = df[col].unique()[0]

    return df_rotated

def rename_df_columns_from_DaVis_to_standard(df):
    DaVis_naming_dict= {
          "x"  : "x",
          "y"  : "y",
          "z"  : "z",
          "Vx" : "u",
          "Vy" : "v",
          "Vz" : "w",
          }

    df.columns = [
        DaVis_naming_dict[col] for col in df.columns
    ]

    return df

def find_nearest(to_point,from_array):
   """ Finds the nearest available value in a array to a given value

   Inputs:
      to_point: value to find the nearest to in the array
      from_array: array of available values of which the nearest has to be 
      found
   Returns:
      The nearest value found in the array
      The difference between the requested and available closest value 
      in the array
   """
   from numpy import ones,argmin
   deltas = ones(len(from_array))*1000
   for v,i in zip(from_array,range(len(from_array))):
      deltas[i] = abs(to_point - v)

   return from_array[argmin(deltas)]

def regrid_df(df,variables = [],resolution=[0]):
    import numpy as np
    from scipy.interpolate import griddata
    import pandas as pd

    string_columns = df.select_dtypes(include=['object'])

    if not len(resolution)==2:
        grid_y, grid_x = np.mgrid[
                    df['y'].min() : df['y'].max() : resolution[0],
                    df['x'].min() : df['x'].max() : resolution[0],
                    ]
    else:
        grid_y, grid_x = np.mgrid[
                    df['y'].min() : df['y'].max() : resolution[1]*1j,
                    df['x'].min() : df['x'].max() : resolution[0]*1j,
                    ]
        
    df_interpolated = pd.DataFrame({
            'x' : grid_x.ravel(),
            'y' : grid_y.ravel(),
            })

    if not len(variables):
        variables = [f for f in df.columns if not f == 'x' and not f == 'y']
    for v in variables:
        grid_var = griddata(
                ( df['x'].values , df['y'].values) , 
                df[v].values, 
                (grid_x,grid_y),
                method='linear'
                )

        df_interpolated[v] = grid_var.ravel()
        df_interpolated = df_interpolated.fillna(0)
        
    # Re-center the array to the TE location at (0,0) ##########################
    df_interpolated.y = df_interpolated.y - \
            find_nearest(0,df_interpolated.y.values)
    df_interpolated.x = df_interpolated.x - \
            find_nearest(0,df_interpolated.x.values)
    ############################################################################

    for col in string_columns.columns:
        df_interpolated[col] = df[col].unique()[0]

    return df_interpolated

def read_raw_tecplot_folder_and_write_pandas_hdf5(
    case_folder,
    root                  = 0,
    output_file           = 0,
    output_root           = 0,
    overwrite             = False,
):
    from os.path import isfile,join,splitext
    from os import listdir
    from progressbar import ProgressBar,Percentage,Bar,ETA,SimpleProgress
    from pandas import DataFrame

    # File related things ######################################################
    if not output_file:
        output_file = case_folder+"_Aligned.hdf5"

    if not output_root:
        output_root = '/media/carlos/6E34D2CD34D29783/' +\
                '2015-02_SerrationPIV/TR_Data_Location_Calibrated_Article3'

    if not output_file.endswith('_Aligned.hdf5'):
        output_file = output_file.replace("_Aligned.hdf5","")+"_Aligned.hdf5"

    if isfile(join( output_root, output_file )) and not overwrite:
        print "  Exiting; file exists:\n{0}".format(output_file)
        return 0
    else:
        print "  Writing\n{0}".format(output_file)

    # ##########################################################################


    time_step_files = sorted(
        [join(root,case_folder,f) for f in listdir(join( root, case_folder )) \
         if splitext(f)[1] == '.dat']
    )

    progress = ProgressBar(
         widgets=[
             Bar(),' ',
             Percentage(),' ',
             ETA(), ' (file ',
             SimpleProgress(),')'], 
         maxval=len(time_step_files)
         ).start()

    cnt = 0
    for f,t in zip(time_step_files,range(len(time_step_files))):

       df_t = read_tecplot_file(
           tecplot_folder         = join( root, case_folder ),
           tecplot_time_step_file = f,
           time_step              = t,
       )

       if cnt == 0:
           df = df_t.copy()
       else:
           df = df.append( df_t, ignore_index = True)

           try:
               x_cnt = df.x.value_counts().max() 
           except AttributeError:
               print df
               raise
           if not x_cnt.max() == x_cnt.min():
               print "  There's something wrong, counted {0} instances of x"\
                       .format(x_cnt.max())
               return 0

       if cnt == 50:

           df = correct_df_translation_rotation( df )\
                   [['x','y','t','u','v','w']]

           df = df.sort_values( by = ['x','y','t'] )

           df.set_index( ['x','y'], inplace = True)

           if t == cnt: append = False; mode = 'w'
           else: append = True; mode = 'a'

           df.to_hdf( join( output_root, output_file ), 
                     'data', format = 't', append = append ,
                     mode = mode)

           cnt = 0

           df = DataFrame()

       cnt += 1

       progress.update(t)

    progress.finish()

    return 1

def run_process( root = 0 , overwrite = False ):
    from os import listdir
    from os.path import isdir,join

    if not root:
        root = '/media/carlos/6E34D2CD34D29783/2015-02_SerrationPIV/'\
                +'TR_Data_NewProcessing'

    case_folders = [f for f in listdir(root) \
                    if isdir( join( root , f ) ) \
                    and f.endswith('_tr')
                   ]

    for cf in case_folders:

        read_raw_tecplot_folder_and_write_pandas_hdf5(
            case_folder           = cf,
            root                  = root,
            output_file           = 0,
            overwrite             = overwrite,
        )
