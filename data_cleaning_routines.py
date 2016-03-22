DaVis_orientaion_to_standard_dict = {
    'x':                                 'y',
    'y':                                 'x',
    'v':                                 'u',
    'u':                                 'v',
    'w':                                 'w',
    'u_rms':                             'v_rms',
    'v_rms':                             'u_rms',
    'w_rms':                             'w_rms',
    "Reynold_stress_uv":                 'Reynold_stress_uv',
    "Reynold_stress_uw":                 "Reynold_stress_vw",
    "Reynold_stress_vw":                 "Reynold_stress_uw",
    "Reynold_stress_uu":                 "Reynold_stress_vv",
    "Reynold_stress_vv":                 'Reynold_stress_uu',
    "Reynold_stress_ww":                 'Reynold_stress_ww',
    "Length of Avg V":                   "Length of Avg V",
    "Length of Standard deviation of V": "Length of RMS V",
    "Length of RMS V":                   "Length of RMS V",
    "case_name":                         "case_name",
}

def read_davis_tecplot_folder_and_rotate_to_serration_surface(tecplot_folder,
                                                             plot = False):
    """Reads in a tecplot folder, given, and returns a pandas data frame
    that is rotated to the airfoil surface as specified in the mask

    """
    import pandas as pd 
    import re
    from os import listdir
    from os.path import split,join
    from Masks import Masks as masks
    from numpy import array
    from corrections import correction_dict

    tecplot_files = [f for f in listdir(tecplot_folder) \
                     if f.endswith('.dat')]

    for tecplot_file in tecplot_files:

        # Get available variables
        f = open(join(tecplot_folder,tecplot_file),'ro')
        var_flag = False
        dev_flag = False
        for line in f:
            string = re.findall("^VARIABLES[ _A-Za-z0-9,\"=]+",line)
            if string:
                variables = [v.replace(' ','_').replace("\"","") \
                             for v in string[0].replace("VARIABLES = ",'')\
                             .split(", ")]
                variables = [v for v in variables if len(v)]
                var_flag = True
            string = re.findall("^TITLE = [ -_A-Za-z0-9,\"=]+",line)
            if string:
                dev_flag = True
            if var_flag and dev_flag:
                break
        f.close()

        # Put the data into a data frame #######################################
        data = pd.read_table(
                join(tecplot_folder,tecplot_file),
                skiprows  = 4,
                names     = variables,
                sep       = '[ \t]+',
                index_col = False,
                engine    = 'python'
                )
        if tecplot_file == tecplot_files[0]:
            df = data.copy()
        else:
            variable = data.drop(['x','y'],1).columns[0]
            df[variable] = data[variable]
        # ######################################################################

    df = rename_df_columns_from_DaVis_to_standard(df)

    trailing_edge,phi,alpha,U,z = decript_case_name(split(tecplot_folder)[-1])

    if trailing_edge == 'serrated': device = 'Sr20R21'
    elif trailing_edge == 'straight': device = 'STE'
    elif trailing_edge == 'slitted': device = 'Slit20R21'
    if 'tr' in tecplot_folder: time_resolved = True
    else: time_resolved = False

    df['case_name'] = "{0}_phi{1}_alpha{2}_U{3}_loc{4}.dat".format(
        device, phi, alpha, U, z
    )

    if time_resolved:
        df.case_name = df.case_name.unique()[0].replace('.dat','_tr.dat')

    try:
        mask = array(masks[df.case_name.unique()[0]])
    except KeyError:
        mask = array(masks[df.case_name.unique()[0].replace('_tr','')])

    if split( tecplot_folder )[1] in correction_dict.keys():
        x_corr, y_corr, angle_corr = correction_dict[
        split( tecplot_folder )[1]
        ]
    else:
        print "    Didn't find a correction term for {0}".format( 
            split( tecplot_folder )[1] 
        )
        print "    among"
        print "    {0}".format( correction_dict.keys() )
        x_corr, y_corr, angle_corr = (0, 0, 0)

    df.x = df.x - mask[1,0] + y_corr 
    df.y = df.y - mask[1,1] - x_corr 

    data_shape = (
        len(df.y.unique()),len(df.x.unique())
    )

    mask_rotation = get_angle_between_points(
        mask[1], mask[2]
    )

    if 'STE' in df.case_name.unique()[0] or 'loc10' in df.case_name.unique()[0]:
        airfoil_angle_correction = -11
    else:
        airfoil_angle_correction = 0
    df = rotate_df( df, mask_rotation + angle_corr + airfoil_angle_correction )

    # Regrid ###################################################################
    df = regrid_df(
        df, 
        variables = df.columns,
        resolution = data_shape
    )
    ######################################################################

    # Rotate variables to the standard coordinate system ###################
    #rotated_variables = [ DaVis_orientaion_to_standard_dict[v] \
    #                     for v in df.columns]
    #df.columns = rotated_variables
    # ######################################################################

    if plot:
        show_sample_bls_from_df(df)

    df = get_vorticity( df )

    df.to_pickle("averaged_data/{0}.p".format(df.case_name.unique()[0]))


def show_sample_bls_from_df( df ):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib              import rc
    from numpy                   import meshgrid
    from boundary_layer_routines import return_bl_parameters

    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=1.6,style='ticks',
            rc={"axes.axisbelow": False})
    rc('font',family='serif', serif='cm10')

    # At 0% 2h:
    if 'STE' in df.case_name.unique()[0] or 'loc10' in df.case_name.unique()[0]:
        pct0 = find_nearest( -0.05*40, df.x  )
    else:
        pct0 = find_nearest(  0.05*40, df.x  )
    # At 10% 2h:
    pct80 = find_nearest( 0.75*40., df.x )

    trailing_edge,phi,alpha,U,z = \
            decript_case_name( df.case_name.unique()[0] )

    if trailing_edge == 'serrated': device = 'Sr20R21'
    elif trailing_edge == 'straight': device = 'STE'
    elif trailing_edge == 'slitted': device = 'Slit20R21'

    case_name = "{0}_a{1}_p{2}_U20_z{3:02.0f}_tr".\
            format( device, alpha, phi, float(z)*20 )

    bl_params = return_bl_parameters( 
        case_name,
        [pct0] )

    X,Y = meshgrid( df.x.unique() , df.y.unique() )
    U   = df['u'].reshape(X.shape)

    V   = df['v'].reshape(X.shape)

    fig,axes = plt.subplots(1,2, figsize = (20,7))

    cf = axes[0].contourf(X,Y,U)

    axes[0].plot( 
        df[ ( df.x == pct0  ) ].x, 
        df[ ( df.x == pct0  ) ].y,
        lw = 3, c = 'k' 
    )
    axes[0].plot( 
        df[ ( df.x == pct80  ) ].x, 
        df[ ( df.x == pct80  ) ].y,
        lw = 3, c = 'k' 
    )

    d_99 = bl_params.delta_99.values[0]
    axes[1].axhline( d_99 , lw = 2, c = 'r' )

    axes[0].axhline( 0 , lw = 3 , c = 'k' )

    stride = 10

    axes[0].quiver(
        X[::stride, ::stride],
        Y[::stride, ::stride],
        U[::stride, ::stride],
        V[::stride, ::stride],
        scale=300
    )

    axes[0].set_aspect('equal')

    axes[1].plot( 
        df[ ( df.x == pct0 ) & ( df.u > 0.1 ) ].u, 
        df[ ( df.x == pct0 ) & ( df.u > 0.1 ) ].y
    )
    axes[1].plot( 
        df[ ( df.x == pct80 ) & ( df.u > 0.1 ) ].u, 
        df[ ( df.x == pct80 ) & ( df.u > 0.1 ) ].y
    )

    axes[1].set_xlim(0,20)
    axes[1].axhline( 0 , lw = 3 , c = 'k' )

    plt.colorbar(cf)
    plt.show()
    plt.close(fig)


def get_vorticity(df):
    from numpy import shape,zeros
    from sys import exit

    if "vorticity_xy" in df.columns:
        # Do nothing and return the same DF
        return df

    # Get shape of 2D plane
    nx = len(df['x'].unique())
    ny = len(df['y'].unique())
    Ux = df['u'].values.reshape((ny,nx))
    Uy = df['v'].values.reshape((ny,nx))
    ax = df['x'].values.reshape((ny,nx))/1000. # [mm] -> [m]
    ay = df['y'].values.reshape((ny,nx))/1000. # [mm] -> [m]

    i,j = shape(Ux)

    # Sanity check:
    if i != shape(Uy)[0] or i != shape(ax)[0] or i != shape(ay)[0]:
        exit("   The shape of the matrices while getting the "+\
             "vorticity is not the same!")

    duy_dax = zeros((i,j))
    dux_day = zeros((i,j))
    for ii in range(1,i-1):
        for jj in range(1,j-1):
            duy_dax[ii,jj] = (Uy[ii,jj+1]-Uy[ii,jj-1])\
                    /(ax[ii,jj+1]-ax[ii,jj-1])
    for ii in range(1,i-1):
        for jj in range(1,j-1):
            dux_day[ii,jj] = (Ux[ii+1,jj]-Ux[ii-1,jj])\
                    /(ay[ii+1,jj]-ay[ii-1,jj])

    vorticity = duy_dax - dux_day

    df['vorticity_xy'] = vorticity.ravel()

    return df


def rename_df_columns_from_DaVis_to_standard(df):
    DaVis_naming_dict= {
          "x"                                 : "x",
          "y"                                 : "y",
          "Avg_Vx"                            : "u",
          "Avg_Vy"                            : "v",
          "Avg_Vz"                            : "w",
          "RMS_Vx"                            : "u_rms",
          "RMS_Vy"                            : "v_rms",
          "RMS_Vz"                            : "w_rms",
          "Reynold_stress_XY"                 : "Reynold_stress_uv",
          "Reynold_stress_XZ"                 : "Reynold_stress_uw",
          "Reynold_stress_YZ"                 : "Reynold_stress_vw",
          "Reynold_stress_XX"                 : "Reynold_stress_uu",
          "Reynold_stress_YY"                 : "Reynold_stress_vv",
          "Reynold_stress_ZZ"                 : "Reynold_stress_ww",
          "Length_of_Avg_V"                   : "Length of Avg V",
          "Length_of_Standard_deviation_of_V" : "Length of RMS V",
          "Length_of_RMS_V"                   : "Length of RMS V",
          }
    DaVis_naming_dict2= {
          "x"                                 : "x",
          "y"                                 : "y",
          "Avg Vx"                            : "u",
          "Avg Vy"                            : "v",
          "Avg Vz"                            : "w",
          "RMS Vx"                            : "u_rms",
          "RMS Vy"                            : "v_rms",
          "RMS Vz"                            : "w_rms",
          "Reynold stress_XY"                 : "Reynold_stress_uv",
          "Reynold stress_XZ"                 : "Reynold_stress_uw",
          "Reynold stress_YZ"                 : "Reynold_stress_vw",
          "Reynold stress_XX"                 : "Reynold_stress_uu",
          "Reynold stress_YY"                 : "Reynold_stress_vv",
          "Reynold stress_ZZ"                 : "Reynold_stress_ww",
          "Length of_Avg_V"                   : "Length of Avg V",
          "Length of_Standard_deviation_of_V" : "Length of RMS V",
          "Length of_RMS_V"                   : "Length of RMS V",
          }

    if not "Avg_Vx" in df.columns:
        return df

    df.rename( columns = DaVis_naming_dict  , inplace = True )
    df.rename( columns = DaVis_naming_dict2 , inplace = True )


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

def regrid_df(df,variables,resolution=[0]):
    import numpy as np
    from scipy.interpolate import griddata
    import pandas as pd

    string_columns = df.select_dtypes(include=['object'])

    if not len(resolution)==2:
        grid_y, grid_x = np.mgrid[
                    int(df['y'].min()):int(df['y'].max()):0.1,
                    int(df['x'].min()):int(df['x'].max()):0.1,
                    ]
    else:
        grid_y, grid_x = np.mgrid[
                    int(df['y'].min())-1:int(df['y'].max()+1):resolution[1]*1j,
                    int(df['x'].min())-1:int(df['x'].max()+1):resolution[0]*1j,
                    ]

    df_interpolated = pd.DataFrame({
            'x' : grid_x.ravel(),
            'y' : grid_y.ravel(),
            })
    
    variables = [v for v in variables if not v in string_columns \
                 and not v == 'x' and not v == 'y']
    for v in variables:
        grid_var = griddata(
                (df['x'].values,df['y'].values), 
                df[v].values, 
                (grid_x,grid_y),
                method='cubic'
                )
        df_interpolated[v] = grid_var.ravel()


    # Re-center the array to the TE location at (0,0)
    df_interpolated.y = df_interpolated.y - \
            find_nearest(0,df_interpolated.y.values)
    df_interpolated.x = df_interpolated.x - \
            find_nearest(0,df_interpolated.x.values)

    if len(string_columns):
        for col in string_columns.columns:
            df_interpolated[col] = df[col].unique()[0]

    if len(string_columns):
        for sc in string_columns:
            df_interpolated[sc] = df[sc].unique()[0]

    df_interpolated = df_interpolated.fillna(0)
    return df_interpolated

def rotate_polygon(polygon,theta):
    """Rotates the given polygon which consists of corners 
    represented as (x,y),
    around the ORIGIN, clock-wise, theta degrees
    
    """
    
    import math
    from numpy import array
    theta = math.radians(-theta)
    rotatedPolygon = []

    for corner in polygon.T :
        rotatedPolygon.append(
            ( 
                corner[0]*math.cos(theta)-corner[1]*math.sin(theta) , 
                corner[0]*math.sin(theta)+corner[1]*math.cos(theta) 
            ) 
        )

    return array(rotatedPolygon).T

def mask_orientation_from_DaVis_to_standard(mask):
    from numpy import array

    mask_x_rot =  mask[1,:]
    mask_y_rot = -mask[0,:]

    return array([mask_x_rot , mask_y_rot])

def return_mask(case_name,zero_at_TE = True, coordinates = 'standard'):
    from Masks import Masks as masks
    from numpy import array

    mask = array(masks[case_name])

    if zero_at_TE:
        mask = mask - mask[1]
    mask = mask.T
    if not coordinates == 'DaVis':
        mask = mask_orientation_from_DaVis_to_standard(mask)

    return mask

def get_angle_between_points(point1, point2):
    from math import atan, pi, degrees

    delta_y = point2[1] - point1[1]
    delta_x = point2[0] - point1[0]
    if delta_x == 0:
        return pi

    return degrees(atan( delta_y / delta_x ))



def show_surface_from_df(df, variable='u', points = [], mask = []):
    import matplotlib.pyplot as plt
    from numpy import meshgrid,linspace

    if len(df.x.unique())>1000:

        df = regrid_df( df , variables = [variable]+['u','v'], 
                       resolution = [0.1])

    X,Y = meshgrid( df.x.unique() , df.y.unique() )
    Z   = df[variable].reshape(X.shape)

    U   = df['u'].reshape(X.shape)
    V   = df['v'].reshape(X.shape)

    fig,axes = plt.subplots(1,1)

    if variable == 'flow_angle':
        df.flow_angle.fillna(90)
        levels = linspace(-15,15)
        cf = axes.contourf(X,Y,Z,levels=levels)

    else:
        cf = axes.contourf(X,Y,Z)

    stride = 20
    axes.quiver(
        X[::stride, ::stride],
        Y[::stride, ::stride],
        U[::stride, ::stride],
        V[::stride, ::stride],
        scale=300
    )

    if len(mask):
        axes.plot(mask[0,:],mask[1,:])

    if len(points):
        for p in points:
            px = find_nearest( p[0], df.x.unique() )
            py = find_nearest( p[1], df.y.unique() )
            axes.scatter( px, py )

    plt.colorbar(cf)
    axes.set_aspect('equal')

    plt.show()
    plt.close(fig)

def show_streamlined_surface_from_df(df, variable='u', 
                                     points = [], mask = [],
                                     height_correction     = 0,
                                     streamwise_correction = 0,
                                     angle_correction      = 0,
                                     x_max                 = 0,
                                     x_min                 = 0,
                                     y_max                 = 0,
                                     y_min                 = 0,
                                     plot_name             = 'SurfacePlot.png',
                                     airfoil_rotate        = 0,
                                    ):
    import matplotlib.pyplot as plt
    from numpy import meshgrid,linspace,array,append
    from matplotlib import rc
    import seaborn as sns

    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=1.6,style='whitegrid')
    rc('font',family='serif', serif='cm10')


    data_shape = (
        len(df[ (df.x >= x_min) & (df.x <= x_max) ].x.unique()),
        len(df[ (df.y >= y_min) & (df.y <= y_max) ].y.unique())
    )

    df = correct_flow_plane_df(
        df, 
        rotation_angle        = angle_correction,
        height_correction     = height_correction,
        streamwise_correction = streamwise_correction,
    )

    if x_min or x_max:
        df = df[ (df.x >= x_min) & (df.x <= x_max) ]
    if y_min or y_max:
        df = df[ (df.y >= y_min) & (df.y <= y_max) ]

    df = regrid_df(df,variables=list(set(['u','v',variable])),
                   resolution = data_shape)

    X,Y = meshgrid( df.x.unique() , df.y.unique() )
    Z   = df[variable].values.reshape(X.shape)
    U   = df['u'].values.reshape(X.shape)
    V   = df['v'].values.reshape(X.shape)

    fig,axes = plt.subplots(1,1)

    if variable == 'flow_angle':
        df.flow_angle.fillna(90)
        levels = linspace(-15,15)
        cf = axes.contourf(X,Y,Z,levels=levels)

    else:
        cf = axes.contourf(X,Y,Z,
                          cmap = plt.get_cmap('Reds'))

    if len(mask):
        axes.plot(mask[0,:],mask[1,:])

    if len(points):
        axes.scatter(points[0],points[1])

    axes.streamplot(
        X, Y, U, V, 
        density      = 0.5,
        zorder       = 1,
    )

    naca0018_profile = NACA0018_trailing_edge_profile()

    axes.fill_between(
        append(array(naca0018_profile[0]),[df.x.max()]),
        append(array(naca0018_profile[1]),[naca0018_profile[1][-1]]),
        array([df.y.min()]*(len(naca0018_profile[1])+1)),
        facecolor='k',
        zorder = 10
    )

    axes.set_xlabel('$x$ [mm]')
    axes.set_ylabel('$y$ [mm]')

    clb = plt.colorbar(cf)
    clb.set_label('$u$ [m/s]')

    axes.set_aspect('equal')

    if x_min or x_max:
        plt.xlim(x_min, x_max) 
    if y_min or y_max:
        plt.ylim(y_min, y_max)
    
    if not plot_name:
        plt.show()
    else:
        plt.savefig(plot_name,bbox_inches='tight')
    plt.close(fig)

def rotate_df(df,degrees = 0):
    from math import radians
    from numpy import sin, cos, sqrt
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

    # RMS components
    df_rotated['Reynold_stress_uu'] = \
            df['Reynold_stress_uu']*cos(angle)**2\
            + df['Reynold_stress_vv']*sin(angle)**2\
            + 2.*cos(angle)*sin(angle)*df['Reynold_stress_uv']

    df_rotated['Reynold_stress_vv'] = \
            df['Reynold_stress_uu']*sin(angle)**2\
            + df['Reynold_stress_vv']*cos(angle)**2\
            - 2.*cos(angle)*sin(angle)*df['Reynold_stress_uv']

    #df_rotated['Reynold_stress_uv'] = \
    #        df['Reynold_stress_uv']*cos(2*angle) \
    #        + (cos(angle)*sin(angle)) \
    #        * ( df['Reynold_stress_uu']**2 * -df['Reynold_stress_vv']**2 )
    df_rotated['Reynold_stress_uv'] = \
            df['Reynold_stress_uv']*cos(2*angle) \
            + (cos(angle)*sin(angle)) \
            * ( df['Reynold_stress_uu'] - df['Reynold_stress_vv'] )

    df_rotated['u_rms'] = sqrt(
        df_rotated['Reynold_stress_uu'].fillna(0).values
    )
    df_rotated['v_rms'] = sqrt(
        df_rotated['Reynold_stress_vv'].fillna(0).values
    )
    df_rotated['u_rms'] = df_rotated['u_rms'].fillna(0)
    df_rotated['v_rms'] = df_rotated['v_rms'].fillna(0)

    ###

    for col in df.select_dtypes(include=['object']).columns:
        df_rotated[col] = df[col].unique()[0]

    return df_rotated

def decript_case_name(case_name):
    from re import findall

    try:
        trailing_edge = findall('[litSrTE201R]+', case_name)[0]
    except IndexError:
        print case_name
    if trailing_edge == "STE":       trailing_edge = 'straight'
    if trailing_edge == "Sr20R21":   trailing_edge = 'serrated'
    if trailing_edge == "Slit20R21": trailing_edge = 'slitted'

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

def correct_flow_plane_df( df, rotation_angle = 0, 
                          height_correction = 0, 
                          streamwise_correction = 0
                         ):
    # Do all the plane correcions ########################################
    if rotation_angle:
        df = rotate_df( df, rotation_angle )
    ######################################################################

    # Correct height and ignore below minimum trustable y ################
    df.y = df.y+height_correction
    ######################################################################

    # Correct streanwise translation #####################################
    df.x = df.x-streamwise_correction
    ######################################################################

    return df


def get_Reynolds_number(U, C = 0.2):
    rho = 1.193
    mu  = 1.813e-5

    return rho*U*C/mu

