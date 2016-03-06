import data_cleaning_routines as dcr
from checklist import checklist

def tecplot_dat_to_pandas_p(folder, plot = False):
    import os
    for tecplot_folder in [f for f in os.listdir( folder )]:
        if not checklist[tecplot_folder]:
            print "  Processing {0}".format(tecplot_folder)
            dcr.read_davis_tecplot_folder_and_rotate_to_serration_surface(
                tecplot_folder = os.path.join(folder,tecplot_folder), 
                plot = plot
            )

            return

tecplot_dat_to_pandas_p(
    '/home/carlos/Documents/PhD/Articles/Article_3/Scripts/time_resolved'+\
    '/tecplot_data_avg',
    plot = True
)
