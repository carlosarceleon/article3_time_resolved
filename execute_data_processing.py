import data_cleaning_routines as dcr
from checklist import checklist
from raw_data_processing_routines import run_process
from boundary_layer_routines import run_bl_analysis

def average_tecplot_dat_to_pandas_p(folder, plot = False, 
                                    ignore_checklist = False):
    import os
    for tecplot_folder in [f for f in os.listdir( folder )]:
        if not checklist[tecplot_folder] or ignore_checklist:
            print "  Processing {0}".format(tecplot_folder)
            dcr.read_davis_tecplot_folder_and_rotate_to_serration_surface(
                tecplot_folder = os.path.join(folder,tecplot_folder), 
                plot = plot
            )

            if not ignore_checklist:
                return

#average_tecplot_dat_to_pandas_p(
#    '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/time_resolved'+\
#    '/tecplot_data_avg',
#    plot = False,
#    ignore_checklist = True
#)

# The pre-processing of the data, from raw TECPLOT to an aligned data frame ####
#run_process( overwrite = False )
# ##############################################################################

run_bl_analysis()
