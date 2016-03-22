from data_cleaning_routines       import \
        read_davis_tecplot_folder_and_rotate_to_serration_surface
from checklist                    import checklist
from raw_data_processing_routines import run_raw_data_collection
from boundary_layer_routines      import run_bl_analysis
from data_extraction_routines     import run_data_extraction
from data_analysis_routines       import do_the_time_resolved_analysis

def average_tecplot_dat_to_pandas_p(folder, plot = False, 
                                    ignore_checklist = False):
    import os
    for tecplot_folder in [f for f in os.listdir( folder )]:
        if not checklist[tecplot_folder] or ignore_checklist:
            print "  Processing {0}".format(tecplot_folder)
            read_davis_tecplot_folder_and_rotate_to_serration_surface(
                tecplot_folder = os.path.join(folder,tecplot_folder), 
                plot = plot
            )

            if not ignore_checklist:
                return

# 1) ###########################################################################
# Extract the averaged time-resolved results from the TECPLOT files from DaVis
# to a pickled pandas dataframe pickle #########################################
#
#average_tecplot_dat_to_pandas_p(
#    '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/time_resolved'+\
#    '/tecplot_data_avg',
#    plot = True,
#    ignore_checklist = False
#)
# ##############################################################################

# 2) ###########################################################################
# Do the boundary layer analysis on the averaged time-resolved pickles and put
# the boundary layer parameters in a pickle of its own #########################
#
#run_bl_analysis()
#
# ##############################################################################

# 3) ###########################################################################
# The pre-processing of the data, from raw TECPLOT to an aligned data frame ####
#
#run_raw_data_collection( overwrite = True , 
#                        only_for = [
#                            'STE_a0_p0_U20_z00_tr',
#                            'Sr20R21_a0_p0_U20_z10_tr'
#                        ]
#                       )
#
# ##############################################################################

# 4) ###########################################################################
# The data extraction of the interesting coordinate timeseries from the aligned
# data frame created in the previous step ######################################
#
run_data_extraction()
#
# ##############################################################################

# 5) ###########################################################################
# Do the data analysis #########################################################
#
#do_the_time_resolved_analysis()
#
# ##############################################################################

