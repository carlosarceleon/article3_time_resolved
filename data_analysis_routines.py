def remove_angle_jumps(df):
    from numpy import sign
    from math import pi

    for var in ['u','v','w']:
        df['phi_'+var].loc[df['phi_'+var]<0] = \
                df['phi_'+var].loc[df['phi_'+var]<0] + pi

        for ix in range(len(df))[:-2]:
            dif = df['phi_'+var].ix[ix+1] - df['phi_'+var].ix[ix]
            if abs(dif) > pi*0.4:
                df['phi_'+var].ix[ix+1] = df['phi_'+var].ix[ix+1] \
                        - sign(dif) * pi

        df['phi_'+var].loc[df['phi_'+var]<0] = \
                df['phi_'+var].loc[df['phi_'+var]<0] + pi

        for ix in range(len(df))[:-2]:
            dif = df['phi_'+var].ix[ix+1] - df['phi_'+var].ix[ix]
            if abs(dif) > pi*0.4:
                df['phi_'+var].ix[ix+1] = df['phi_'+var].ix[ix+1] \
                        - sign(dif) * pi

    return df

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    import numpy as np
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def calculate_Uc_from_St(df,delta_x):
    from math import pi
    from scipy.stats import linregress
    import pandas as pd
    from numpy import nan

    for var in ['u','v']:
        df = df.append(
            pd.DataFrame( data = {
                'St_'+var:   0,
                'phi_'+var: 0
            }, index = [0]
            ), ignore_index = True
        )

        r_value      = 0
        truncated_df = df[['St_'+var, 'phi_'+var]].copy()
        truncated_df = truncated_df[ truncated_df['phi_'+var] != 0 ] 
        truncated_df = truncated_df.sort_values( by = [ 'St_' + var ] )
        truncated_df = truncated_df.dropna()
        truncated_df = truncated_df.reset_index( drop = True )
        consider     = len(truncated_df)-1

        while r_value**2 < 0.97:

            truncated_df = truncated_df.sort_values( by = ['St_' + var] ).\
                    ix[:consider].\
                    reset_index(drop=True)

            slope, intercept, r_value, p_value, std_err = linregress(
                truncated_df['phi_'+var],
                truncated_df['St_'+var],
            )

            consider -= 1
            if len(truncated_df) < 3:
                slope = 0
                intercept = 0
                break
        
        Ue    = df.Ue.unique()[0]
        delta = df.delta_99.unique()[0]
    
        if var == 'u' and slope:
            Uc_u = Ue*2*pi*slope*(delta_x/1000.)/(delta/1000.)
            slope_u = slope
            if Uc_u > 25: Uc_u = 0 
        elif not slope:
            Uc_u = nan
            slope_u = nan

        if var == 'v' and slope:
            Uc_v = Ue*2*pi*slope*(delta_x/1000.)/(delta/1000.)
            slope_v = slope
        elif not slope:
            Uc_v = nan
            slope_v = nan

    return Uc_u, Uc_v, df, intercept, slope_u, slope_v

def do_the_coherence_analysis(df_upstream, df_downstream):
    import pandas as pd
    from raw_data_processing_routines import find_nearest

    coherence_df = pd.DataFrame()
    # Check how many cases there are ###########################################
    cases = df_upstream.case_name.unique()

    for c in cases:

        for near_x_down in df_downstream[df_downstream.case_name == c].near_x\
                           .unique():

            x_down = df_downstream[
                ( df_downstream.case_name == c ) &\
                ( df_downstream.near_x == near_x_down)
            ].x.min()

            x_up = find_nearest( 
                x_down, 
                df_upstream[ df_upstream.case_name == c ].x.unique() 
            )

            y_up = df_upstream[ 
                ( df_upstream.x == x_up ) & \
                ( df_upstream.case_name == c )
            ].y.unique()[0] 

            y_down = find_nearest( 
                y_up, 
                df_downstream[ 
                    ( df_downstream.x == x_down ) & \
                    ( df_downstream.case_name == c )  
                ].y.unique()
            )

            print c, y_up, y_down, x_down

            coherence_df = coherence_df.append(
                get_Uc_phi_and_coherence( 
                    df_upstream[ 
                        ( df_upstream.case_name == c ) & \
                        ( df_upstream.x == x_up )      & \
                        ( df_upstream.y == y_up ) 
                    ] , 

                    df_downstream[ 
                        ( df_downstream.case_name == c ) & \
                        ( df_downstream.x == x_down )    & \
                        ( df_downstream.y == y_down ) 
                    ], 
                ),
                ignore_index = True
            )

    return coherence_df


def get_Uc_phi_and_coherence( signal1_df, signal2_df, case_name ):
    import pandas as pd
    from scipy.signal import csd
    from numpy import abs,sqrt,angle
    from boundary_layer_routines import return_bl_parameters

    max_lag          = 10000
    nperseg          = 2**6
    fs               = 10000
    x_1              = signal1_df.x.unique()
    x_2              = signal2_df.x.unique()
    y_1              = signal1_df.y.unique()
    y_2              = signal2_df.y.unique()

    delta_x = sqrt( abs(x_1[0] - x_2[0])**2 + abs(y_1[0] - y_2[0])**2 )

    df = pd.DataFrame()

    if not signal1_df.u.mean() or not signal1_df.v.mean()\
       or not signal1_df.w.mean():
        return df
    if not signal2_df.u.mean() or not signal1_df.v.mean()\
       or not signal2_df.w.mean():
        return df

    for var in ['u', 'v','w']:

        # Get the perturbations ################################################
        s1 = signal1_df[var].values[0:max_lag] 

        s2 = signal2_df[var].values[0:max_lag] 
        # ######################################################################

        f,Pxy = csd(
            s2,s1,
            nperseg = nperseg,
            fs      = fs,
        )

        f,Pxx = csd(
            s1,s1,
            nperseg = nperseg,
            fs      = fs,
        )

        f,Pyy = csd(
            s2,s2,
            nperseg = nperseg,
            fs      = fs,
        )

        gamma_squared = abs( Pxy )**2 / ( Pxx * Pyy )

        gamma = sqrt(gamma_squared)

        #Phi = arctan( Pxy.imag / Pxy.real )
        Phi = angle( Pxy )

        data = pd.DataFrame()
        data['f_'+var]     = f[:-1]
        data['phi_'+var]   = Phi[:-1]
        data['gamma_'+var] = gamma[:-1]
        data['mean_'+var]  = signal2_df[var].values[0:max_lag].mean()
        data['std_'+var]   = signal2_df[var].values[0:max_lag].std()

        #data = data[
        #    data['f_'+var] >= freq_lower_limit
        #].reset_index( drop = True )

        df = pd.concat( [ df, data ], axis = 1 )

        bl_data = return_bl_parameters( case_name , [ x_2 ] )

        df['delta_99'] = bl_data.delta_99.unique()[0]
        df['Ue']       = bl_data.Ue.unique()[0]

        df['St_'+var] = get_Strouhal( 
            df['f_'+var], bl_data.delta_99.values[0], 
            bl_data.Ue.values[0] 
        )


    df = remove_angle_jumps(df)
    df = df.drop_duplicates()

    df.loc[ is_outlier( df[ 'phi_v' ] , thresh = 2.0) , 'phi_v' ] = 0
    df.loc[ is_outlier( df[ 'phi_u' ] , thresh = 2.0) , 'phi_u' ] = 0

    Uc_u, Uc_v,data,intercept,slope_u,slope_v = calculate_Uc_from_St(
        df,
        delta_x = delta_x
    )

    df['Uc_u']                  = Uc_u
    df['Uc_v']                  = Uc_v
    df['trusted_f_v']           = data.f_v
    df['trusted_f_u']           = data.f_u
    df['slope_u']               = slope_u**-1
    df['slope_v']               = slope_v**-1
    df['x_upwind']              = signal1_df.x.unique()[0]
    df['near_x_upwind']         = signal1_df.near_x.unique()[0]
    df['y_upwind']              = signal1_df.y.unique()[0]
    df['near_y_upwind']         = signal1_df.near_y.unique()[0]
    df['near_y_delta_upwind']   = signal1_df.near_y_delta.unique()[0]
    df['x_downwind']            = signal2_df.x.unique()[0]
    df['near_x_downwind']       = signal2_df.near_x.unique()[0]
    df['y_downwind']            = signal2_df.y.unique()[0]
    df['near_y_downwind']       = signal2_df.near_y.unique()[0]
    df['near_y_delta_downwind'] = signal2_df.near_y_delta.unique()[0]
    df['case']                  = signal1_df.case_name.unique()[0]
    df['delta_x']               = delta_x

    # Kill the obvious outliers ################################################
    df.loc[ ( df.Uc_u > 25 ) , 'Uc_u' ] = 0
    df.loc[ ( df.Uc_u < 15 ) & ( df.y_upwind > 5 ) , 'Uc_u' ] = 0
    df.loc[ ( df.Uc_u > 12 ) & ( df.y_upwind < 3 ) , 'Uc_u' ] = 0
    #df.loc[ df.gamma_u < 0.3 , 'gamma_u' ] = nan
    #df.loc[ df.gamma_v < 0.3 , 'gamma_v' ] = nan
    # ##########################################################################

    return df

def get_color_and_marker(case_name):

    if "STE" in case_name:
        color  = (0.0, 0.4470588235294118, 0.6980392156862745)
        marker = 'x'
        cmap = 'Blues'
    elif 'z00' in case_name:
        color = (0.0, 0.6196078431372549, 0.45098039215686275)
        marker = '2'
        cmap = 'Greens'
    elif 'z05' in case_name:
        color = (0.8352941176470589, 0.3686274509803922, 0.0)
        marker = '+'
        cmap = 'Oranges'
    elif 'z10' in case_name:
        color = (0.8, 0.4745098039215686, 0.6549019607843137)
        marker = 'o'
        cmap = 'RdPu'

    else: print case_name; return 0,0
    return color,marker,cmap
                        
def do_the_reynolds_stress_quadrant_analysis(cases_df,y_delta, plot_name = ''):

    from matplotlib import rc
    import seaborn as sns
    import matplotlib.pyplot as plt
    from boundary_layer_routines import return_bl_parameters

    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=5.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    cases = sorted(cases_df.file.unique(),reverse=True)
    fig,axes = plt.subplots(
        1,len(cases), 
        figsize = (figsize[0]*len(cases), figsize[1]), 
        sharex=True,
        sharey=True, 
    )

    if not len(cases)>2:
        return 0
        
    for case_file, case_i, ax in zip(cases, range(len(cases)), axes):

        case = cases_df[cases_df.file == case_file]\
                .sort_values( by = ['t'] )

        bl_data = return_bl_parameters( case_file , [case.x.unique()[0]])

        color, marker, cmap = get_color_and_marker( 
            case_file
        )

        if not color and not marker:
            break

        case["uprime"] = ( case.u - case.u.mean() ) / bl_data.Ue.values[0]
        case["vprime"] = ( case.v - case.v.mean() ) / bl_data.Ue.values[0]
        
        ax.plot(
            case['uprime'].values[::100],
            case['vprime'].values[::100],
            ls              = '',
            marker          = marker,
            markeredgewidth = markeredgewidth,
            markerfacecolor = markerfacecolor,
            markeredgecolor = 'k',
            markersize      = markersize*1.8,
            mew             = mew,
            color           = 'k',
            alpha           = 0.4
        )

        kde = sns.kdeplot(case.uprime, case.vprime,
                    cmap         = cmap,
                    ax           = ax,
                    shade        = True,
                    shade_lowers = False,
                    gridsize     = 30
                   )
        
        ax.set_aspect('equal')
        kde.collections[0].set_alpha(0)

        ax.set_xlabel('')
        ax.set_ylabel('')

        ax.set_xlim( -0.3 , 0.3 )
        ax.set_ylim( -0.3 , 0.3 )
        ax.set_xticks([-0.2,0,0.2])
        ax.set_yticks([-0.2,0,0.2])

        ax.axhline( 0, ls = '--', lw=3 , c = 'k')
        ax.axvline( 0, ls = '--', lw=3 , c = 'k')

        if y_delta == 0.1:
            ax.set_xlabel(r"$u'/u_e$")
            ax.grid(False)

    axes[0].set_ylabel(r"$v'/u_e$")

    t = axes[0].text( -0.25, 0.2, r'$y/\delta_{{99}} = {0}$'.format(y_delta))
    t.set_bbox(dict(color='white', alpha=0.7, edgecolor='white'))


    plot_name = 'Results/ReynoldsQuadrant_{0}_ydelta{1:.2f}.png'\
        .format( plot_name,y_delta )

    fig.savefig( 
        plot_name.replace('.','_').replace('_png','.png'),
        bbox_inches = 'tight'
    )
    plt.cla()

def get_Strouhal(f,delta,U):
    delta = delta/1000.
    return f*delta/U


def do_the_frequency_analysis(cases_df, y, plot_name = '', schematic = ''):
    
    from scipy.signal import welch
    from numpy import log10
    from matplotlib import rc
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.cbook import get_sample_data
    from boundary_layer_routines import return_bl_parameters
    
    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=3.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    if not len(cases_df):
        print "   No cases were passed to process"
        return 0

    for var in ['u','v','w']:

        figsize = (8,5)
        fig,ax = plt.subplots(1,1, figsize = figsize)
        
        kolmogorov_law_curve = [[0],[0]]

        if len(cases_df.file.unique()) <= 1:
            return 0

        for case_file in cases_df.file.unique():

            case = cases_df[cases_df.file == case_file]\
                    .sort_values( by = ['t'] )

            freq, Pxx = welch(
                x       = case[var].values,
                nperseg = nperseg,
                fs      = fs,
                scaling = 'spectrum',
            )

            case['case'] = case.file.unique()[0]
            bl_data = return_bl_parameters( case.file.unique()[0] , 
                                           [case.x.unique()[0]])

            color, marker, cmap = get_color_and_marker( 
                case_file
            )

            if not color and not marker:
                break

            res = dict( zip( 
                get_Strouhal( freq, bl_data.delta_99.values[0], 
                                 bl_data.Ue.values[0] ),
                10 * log10( Pxx ),
            ))
            
            ax.plot(
                get_Strouhal( freq, bl_data.delta_99.values[0], 
                                 bl_data.Ue.values[0] )[:-1],
                10 * log10( Pxx )[:-1],
                marker          = marker,
                markeredgewidth = markeredgewidth,
                markerfacecolor = markerfacecolor,
                #markeredgecolor = color,
                markersize      = markersize*1.5,
                mew             = mew,
                #color           = color
                label = case_file.replace("_"," ")
            )

            #if "z05" in case_file and not y < 0.25 and var == 'u':
            #    pass

            #    case['edge_normal'] = \
            #            case.u * cos( deg2rad( serration_angle ) ) +\
            #            case.w * sin( deg2rad( serration_angle ) )

            #    freq, Pxx = welch(
            #        x       = case.edge_normal.values,
            #        nperseg = nperseg,
            #        fs      = fs,
            #        scaling = 'spectrum',
            #    )

            #    color = "#3498db"

            #    ax.plot(
            #        get_Strouhal( freq, bl_data.delta_99.values[0], 
            #                         bl_data.Ue.values[0] )[:-1],
            #        10 * log10( Pxx )[:-1],
            #        marker          = 's',
            #        markeredgewidth = markeredgewidth,
            #        markerfacecolor = markerfacecolor,
            #        markeredgecolor = color,
            #        markersize      = markersize*1.5,
            #        mew             = mew,
            #        color           = color,
            #    )


            k_lims = ( sorted(res.keys())[8], sorted(res.keys())[20] )

            if not any(kolmogorov_law_curve[1]) or \
               res[k_lims[0]] > kolmogorov_law_curve[1][0]:

                kolmogorov_law_curve = get_kolmogorov_law_curve(x_lim = k_lims)
                kolmogorov_law_curve[1] = kolmogorov_law_curve[1] \
                        + res[k_lims[0]] + 3

        ax.plot( 
            kolmogorov_law_curve[0] , 
            kolmogorov_law_curve[1] , 
            '--',
            color = 'k' ,
            lw    = 3,
        )

        if not k_lims == (0,0):
            t = ax.text(
                k_lims[0], #+(k_lims[1]+k_lims[0])/2.,
                kolmogorov_law_curve[1][0]-3,
                r"$\textrm{St}_\delta^{-5/3}$",
                ha = 'center',
            )
            t.set_bbox(dict(color='white', alpha=0.7, edgecolor='white'))

        ax.set_xscale('log')

        ax.legend( bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=1, mode="expand", borderaxespad=0.)

        ax.set_xlim( 0.09 , 2.2 )
        ax.set_ylim( -25 , 5 )
        if y == 0.1:
            ax.set_xlabel(r"$\textrm{St}_\delta = f\delta_{99}/u_e$")
        ax.set_ylabel(
            r"$10\log_{10}\left(\Phi_{"+\
            str(var)+r"}\right)$ [dB]"
        )
        #t = ax.text(
        #    1.8, -23,
        #    r'$y/\delta_{{99}} = {0}$'.format(y),
        #    ha = 'right'
        #)
        #t.set_bbox(dict(color='white', alpha=0.7, edgecolor='white'))

        plt.grid(True, which='both')

        if schematic and y == 0.1:
            if var == 'v':
                schematic = schematic.replace('_with_edge_normal',
                                              '_noSTE')
            im = plt.imread( get_sample_data( schematic  ) )
            newax = fig.add_axes([0.175, 0.175, 0.4, 0.4], anchor = 'SW', 
                                         zorder=100)
            newax.imshow(im)
            newax.axis('off')

        plot_composed_name = 'Results/FreqSpectra_{0}_ydelta{1:.2f}_{2}.png'\
            .format( plot_name, y, var )
        print " Going to save\n   {0}".format( plot_composed_name )

        fig.savefig( 
            plot_composed_name.replace('.','_').replace('_png','.png'),
            bbox_inches = 'tight'
        )
        plt.cla()

def get_kolmogorov_law_curve( x_lim = (0.5,1.5) ):
    from numpy import linspace, log10

    slope_x = linspace( x_lim[0], x_lim[1] ,100 )

    slope_y = 10*log10(slope_x**(-5/3.))

    return [ slope_x, slope_y ]

def plot_mean_and_std( df , plot_name = '' ):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import seaborn as sns
    from numpy import array
    from boundary_layer_routines import return_bl_parameters

    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=4.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    for var in ['u','v','w']:
        fig_std,   ax_std   = plt.subplots(1, 1, figsize = figsize)
        fig_mean,  ax_mean  = plt.subplots(1, 1, figsize = figsize)
        for case in df.case_name.unique():
            case_df = df[df.case_name == case]
            bl_data = return_bl_parameters( case , [df.x.unique()[0]])
            y_locs = sorted( case_df.near_y_delta.unique() )
            Um  = []
            std = []
            real_y_locs = []
            for y_loc in y_locs:
                y_df = case_df[ 
                    (case_df.near_y_delta == y_loc) &\
                    (case_df.case_name == case) 
                ]
                Um.append(
                    y_df[var].mean()
                )
                std.append(
                    y_df[var].std()
                )
                real_y_locs.append(
                    case_df[ case_df.near_y_delta == y_loc ]\
                    .y.unique()[0]
                )

            color, marker, cmap = get_color_and_marker( case )

            plot_config = {
                'marker'          : marker,
                'markeredgewidth' : markeredgewidth,
                'markerfacecolor' : markerfacecolor,
                'markeredgecolor' : color,
                'markersize'      : markersize,
                'mew'             : mew,
                'color'           : color,
            }

            ax_std.plot(
                array(std) / bl_data.Ue.values[0], 
                real_y_locs / bl_data.delta_99.values[0],
                ls='',
                **plot_config
            )

            ax_mean.plot(
                array(Um) / bl_data.Ue.values[0], 
                real_y_locs / bl_data.delta_99.values[0],
                ls='',
                **plot_config
            )


        fig_std.savefig(
            "Results/std_{0}_{1}.png".format(
                plot_name.replace('.','_').replace('_png','.png'),
                var,
            ), bbox_inches = 'tight'
        )
        fig_mean.savefig(
            "Results/mean_{0}_{1}.png".format(
                plot_name.replace('.','_').replace('_png','.png'),
                var,
            ), bbox_inches = 'tight'
        )

def plot_coherence_Uc_phi( coherence_df , plot_name = '', schematic = ''):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    import seaborn as sns
    from numpy import array,linspace,append#,sqrt,diag
    from math import pi
    from matplotlib.cbook import get_sample_data
    from scipy.optimize import curve_fit
    from boundary_layer_routines import return_bl_parameters

    def log_law( y, y0, alpha ):
        from numpy import log

        return alpha * log( y / y0 )

    def exponential_eta( phi , eta ):
        from numpy import exp

        return exp( - eta * phi)


    rc('text',usetex=True)
    rc('font',weight='normal')

    sns.set_context('paper')
    sns.set(font='serif',font_scale=3.0,style='whitegrid')
    rc('font',family='serif', serif='Linux Libertine')

    for var in ['u']:#['u','v','w']:

        fig_Uc,    ax_Uc    = plt.subplots(1, 1, figsize = (8,5))

        for case in coherence_df.case.unique():

            case_df = coherence_df[coherence_df.case == case]
        
            bl_data = return_bl_parameters( case , [case_df.x.unique()[0]])

            Uc  = []
            y_locs = sorted( case_df.near_y_delta_downwind.unique() )

            real_y_locs = []
            for y_loc in y_locs:
                y_df = case_df[ 
                    (case_df.near_y_delta_downwind == y_loc) &\
                    (case_df.case == case) 
                ]

                if not var == 'w':
                    Uc.append(
                        y_df['Uc_'+var].unique()[0]
                    )
                real_y_locs.append(
                    case_df[ case_df.near_y_delta_downwind == y_loc ]\
                    .y_downwind.unique()[0]
                )

            real_y_locs = array(real_y_locs)

            color, marker, cmap = get_color_and_marker( case )

            plot_config = {
                'marker'          : marker,
                'markeredgewidth' : markeredgewidth,
                'markerfacecolor' : markerfacecolor,
                'markeredgecolor' : color,
                'markersize'      : markersize,
                'mew'             : mew,
                'color'           : color,
            }

            if not var == 'w':

                clean_Uc = array(Uc)
                clean_y  = array(real_y_locs)

                ax_Uc.plot(
                    clean_Uc[clean_Uc > 0] / bl_data.Ue.values[0], 
                    clean_y[clean_Uc > 0] / bl_data.delta_99.values[0],
                    ls='',
                    **plot_config
                )

                try:
                    popt, pcov = curve_fit( 
                        log_law, 
                        append(
                            clean_y[ ( clean_Uc > 0 ) & ( clean_Uc < 20 ) ],
                            array([0.95,1.0,1.1])*bl_data.delta_99.values[0]
                        ), 
                        append(
                            clean_Uc[ ( clean_Uc > 0 ) & ( clean_Uc < 20 ) ],
                            [1.0*bl_data.Ue.values[0]]*3
                        ),
                    )
                    plot_fit = True
                    #alpha, y0 = sqrt(diag(pcov))
                    #print case
                    #print 100*alpha/ popt[0]
                    #print 100*y0/     popt[1]
                except:
                    plot_fit = False
                    pass

            if plot_fit:

                ax_Uc.plot(
                    log_law( 
                        real_y_locs[1:], 
                        popt[0], 
                        popt[1] 
                    ) / bl_data.Ue.values[0], 
                    real_y_locs[1:] / bl_data.delta_99.values[0],
                    lw = 5,
                    color = 'w',
                )
                ax_Uc.plot(
                    log_law( 
                        real_y_locs[1:], popt[0], popt[1] 
                    ) / bl_data.Ue.values[0], 
                    real_y_locs[1:] / bl_data.delta_99.values[0],
                    lw=2.5,
                    color = color,
                )

        ax_Uc.set_xlabel( r"$u_c/u_e$" )
        ax_Uc.set_ylabel( r"$y/\delta_{{99}}$" )
        ax_Uc.set_xlim(0.5,1.25)
        ax_Uc.set_ylim(0,1.0)
        ax_Uc.set_xticks([0.5, 0.75, 1, 1.25])

        if schematic:
            im = plt.imread( get_sample_data( schematic  ) )
            newax = fig_Uc.add_axes([0.60, 0.12, 0.4, 0.4], anchor = 'SW', 
                                         zorder=100)
            newax.imshow(im)
            newax.axis('off')

        fig_Uc.savefig(
            "Results/Uc_{0}_{1}_xi{2:.2f}.png".format(
                plot_name,
                var,
                case_df.delta_x.unique()[0]
            ).replace('.','_').replace('_png','.png'), 
            bbox_inches = 'tight'
        )

    for y_loc in coherence_df.near_y_delta_downwind.unique():
        print y_loc

        y_df = coherence_df[ coherence_df.near_y_delta_downwind == y_loc ]

        for var in ['u','v']:

            fig_phi,  ax_phi  = plt.subplots(1, 1, figsize = figsize)
            fig_coh,  ax_coh  = plt.subplots(1, 1, figsize = figsize)
            fig_quad, ax_quad = plt.subplots(1, 1, figsize = figsize)

            case_cnt = 0
            for case in y_df.case.unique():

                case_df = y_df[ y_df.case == case ]
                
                bl_data = return_bl_parameters( case , 
                                               [case_df.x.unique()[0]])
                
                color, marker, cmap = get_color_and_marker( case )

                plot_config = {
                    'marker'          : marker,
                    'markeredgewidth' : markeredgewidth,
                    'markerfacecolor' : markerfacecolor,
                    'markeredgecolor' : color,
                    'markersize'      : markersize*2,
                    'mew'             : mew,
                    'color'           : color,
                }

                ax_phi.plot(
                    get_Strouhal( case_df['f_'+var], bl_data.delta_99.values[0],
                                 bl_data.Ue.values[0]),
                    case_df['phi_'+var],
                    ls = '',
                    **plot_config
                )
                #ax_phi.plot(
                #    get_Strouhal( case_df['f_'+var], bl_data.delta_99.values[0],
                #                 bl_data.Ue.values[0] ),
                #    case_df['f_'+var]*case_df['slope_'+var],
                #    '--',
                #    lw = 3,
                #    color = color,
                #)

                ax_coh.plot(
                    case_df['phi_'+var],
                    case_df['gamma_'+var],
                    ls = '',
                    **plot_config
                )

                plot_fit = False
                try:
                    popt, pcov = curve_fit( 
                        exponential_eta, 
                        case_df[
                            ( case_df[ 'gamma_'+var ] > 0.4   ) &\
                            ( case_df[ 'phi_'+var ]   > pi/6. ) 
                        ]['phi_'+var],
                        case_df[
                            ( case_df[ 'gamma_'+var ] > 0.4   ) &\
                            ( case_df[ 'phi_'+var ]   > pi/6. ) 
                        ]['gamma_'+var]
                    )
                    plot_fit = True
                except:
                    plot_fit = False
                    pass

                eta = popt[0]
                if plot_fit:
                    ax_coh.plot(
                        linspace(0,2*pi,30),
                        exponential_eta( linspace(0,2*pi,30), eta ),
                        lw = 5,
                        color = 'w',
                    )
                    ax_coh.plot(
                        linspace(0,2*pi,30),
                        exponential_eta( linspace(0,2*pi,30), eta ),
                        lw = 2.5,
                        color = color,
                    )

                    annotate_ix = case_cnt * 3 + 5
                    ax_coh.annotate(
                        r'$\eta = {0:.2f}$'.format(eta),
                        xy = (
                            linspace(0,2*pi,30)[ annotate_ix + 4],
                            exponential_eta( 
                                linspace(0,2*pi,30), eta 
                            )[ annotate_ix + 4],
                        ),
                        xytext = (
                            linspace(0,2*pi,30)[ annotate_ix + 7],
                            linspace(0,1,30)[::-1][ annotate_ix ],
                        ),
                        arrowprops=dict(
                            arrowstyle="-",
                            connectionstyle="angle,angleA=0,angleB=90,rad=20",
                            lw = 3,
                        ),
                        
                    )




                case_cnt += 1
            # Configure the phi plot ###########################################
            ax_phi.set_yticks(array(
                [0,1/2.,1,3/2.,2,2.5,3 ]
            )*pi)
            ax_phi.set_yticklabels(
                ['$0$','$\\pi/2$','$\\pi$',
                 '$3\\pi/2$','$2\\pi$','$5\\pi/2$','$3\\pi$'
                ]
            )
            ax_phi.set_xlim(0,St_max)
            ax_phi.set_ylim(0,3*pi)

            ax_phi.set_xlabel(r"$\textrm{{St}}_\delta=f\delta_{{99}}/u_e$")
            ax_phi.set_ylabel(
                r"$\phi_{0}$ [rad]"\
                .format(var)
            )
            t = ax_phi.text(
                0.15,5*pi/2., 
                r'$y/\delta_{{99}} = {0}$'.format(
                    case_df.near_y_delta_downwind.unique()[0],
                ))
            t.set_bbox(dict(color='white', alpha=0.7, edgecolor='white'))
            # Save the Phi plot ################################################
            phi_plot_name = "Results/Phi_{0}_y{1:.2f}_{2}_xi{3:.2f}.png".format(
                plot_name.replace('.','_').replace('_png','.png'),
                case_df.near_y_delta_downwind.unique()[0],
                var,
                case_df.delta_x.unique()[0]
            )
            fig_phi.savefig(
                phi_plot_name.replace('.','_').replace('_png','.png'),
                bbox_inches = 'tight'
            )
            # ##################################################################

            # Configure the coherence plot #####################################
            ax_coh.set_xticks(array(
                [0,1/2.,1,3/2.,2]
            )*pi)
            ax_coh.set_xticklabels(
                ['$0$', '$\\pi/2$', '$\\pi$', '$3\\pi/2$', '$2\\pi$' ]
            )
            ax_coh.set_ylim(0,1)
            ax_coh.set_xlim(0,2*pi)
            ax_coh.set_xlabel(
                r"$\phi = \mu_{{x}}\xi,\;\xi \approx {0:.2f}\cdot 2h$".\
                format(case_df.delta_x.unique()[0]/40.)
            )
            ax_coh.set_ylabel(r"$\gamma_{{{0}}}$".format(var))

            t = ax_coh.text(
                pi/8.,0.05,
                r'$y/\delta_{{99}} = {0}$'.format(
                    case_df.near_y_delta_downwind.unique()[0],
                ))

            t.set_bbox(dict(color='white', alpha=0.7, edgecolor='white'))

            if schematic and case_df.near_y_delta_downwind.unique()[0] == 0.5:
                im = plt.imread( get_sample_data( schematic  ) )
                newax = fig_coh.add_axes([0.1, 0.1, 0.3, 0.3], 
                                         anchor = 'SW', zorder=100)
                newax.imshow(im)
                newax.axis('off')

            # Save it ##########################################################
            plot_name_coherence = \
                    "Results/Coherence_{0}_y{1:.2f}_{2}_xi{3:.2f}.png".format(
                        plot_name,
                        case_df.near_y_delta_downwind.unique()[0],
                        var,
                        case_df.delta_x.unique()[0]
                    )
            fig_coh.savefig(
                plot_name_coherence.replace('.','_').replace('_png','.png'),
                bbox_inches = 'tight'
            )
            plt.close(fig_coh)
            plt.close(fig_phi)
            plt.close(fig_Uc)

def do_the_time_resolved_analysis():
    import pandas as pd
    from os.path import split,isfile

    def do_the_frequency_plot(df,plot_name, schematic = ''):
        for y in df.near_y_delta.unique():

            df_y_cases = df[ df.near_y_delta == y ]

            do_the_frequency_analysis( 
                df_y_cases,
                y = y,
                plot_name = plot_name,
                schematic = schematic
            )

    def do_the_Reynolds_quadrant_analysis(df, plot_name):
        for y in df.near_y_delta.unique():

            df_y_cases = df[ df.near_y_delta == y ]

            do_the_reynolds_stress_quadrant_analysis(
                df_y_cases,
                y_delta = y,
                plot_name = plot_name,
            )

    def do_the_coherence_analysis(df_upstream, df_downstream
                                  ,plot_name,schematic = ''):
        coherence_df = pd.DataFrame()
        for y_up, y_down in zip(
            sorted(df_upstream.near_y_delta.unique()),
            sorted(df_downstream.near_y_delta.unique())
        ):

            partial_coherence_df = do_the_coherence_analysis(
                df_upstream[ df_upstream.near_y_delta == y_up ],
                df_downstream[ df_downstream.near_y_delta == y_down ],
            )

            if not partial_coherence_df.empty:
                coherence_df = coherence_df.append( 
                    partial_coherence_df, ignore_index = True 
                )

        plot_coherence_Uc_phi( coherence_df ,
                                 plot_name = plot_name,
                                 schematic = schematic)



    z05_dfs = pd.DataFrame()
    z05_comparison_cases = [
        './ReservedData/Sr20R21_a-12_p0_U20_z05_tr.p',
        './ReservedData/Sr20R21_a0_p0_U20_z05_tr.p',
        #'./ReservedData/Sr20R21_a-12_p6_U20_z05_tr.p',
        #'./ReservedData/Sr20R21_a12_p0_U20_z05_tr.p',
        './ReservedData/Sr20R21_a12_p0_U20_z05_tr.p',
        './ReservedData/Sr20R21_a12_p6_U20_z05_tr.p',
    ]

    for case in z05_comparison_cases:
        if not isfile(case):
            print "   This file doesn't exist:"
            print "     {0}".format(case)
            return 0
        else:
            print "   Appending file:"
            print "     {0}".format(case)
            
        df = pd.read_pickle( case )
        df['file'] = split( case )[1].replace('.p','')
        z05_dfs = z05_dfs.append(
            df[ df.near_x_2h == 0.4] , ignore_index = True
        )

    do_the_frequency_plot( z05_dfs, 'z05',  schematic = '')

    #schematic_TE = '/home/carlos/Documents/PhD/Articles/Article_2/'+\
    #        'Figures/measurement_locations_TE_m2_noSTE.png'

    #TE_cases = TE_cases[ TE_cases.case_name != 'STE_a0_p0_U20_z00_tr' ]

    #do_the_Reynolds_quadrant_analysis( TE_cases, 'TE' )

    #do_the_coherence_analysis( TE_cases_upstream, TE_cases, 'TE' , 
    #                          schematic = schematic_TE)
    #plot_mean_and_std( TE_cases )

    #do_the_frequency_plot( x0_cases, 'x0', schematic = schematic_x0 )
    #do_the_Reynolds_quadrant_analysis( x0_cases, 'x0' )
    #do_the_coherence_analysis( x0_coherence_cases , "x0",
    #                         schematic = schematic_x0)

# Constants ####################################################################

nperseg         = 2**6
fs              = 10000

St_min          = 0.225
St_max          = 2.4

markeredgewidth = 3
markerfacecolor = 'none'
markersize      = 12
mew             = 4 # Marker edge width

figsize         = (8,7)

serration_angle = 76.

# ##############################################################################

