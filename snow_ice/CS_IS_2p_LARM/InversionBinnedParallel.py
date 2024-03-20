
'''
The objective of this program is to provide the entire pipeline from directly consuming 
Carmen's interpolated data 

This time including support for multi-threading.
'''
import numpy as np
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
from itertools import chain
from mpl_toolkits.basemap import Basemap
import os
import warnings
import netCDF4 as nc



def draw_map(m, scale=0.2):
    # draw a shaded-relief image
    m.shadedrelief(scale=scale)
    
    # lats and longs are returned as a dictionary
    lats = m.drawparallels(np.linspace(-90, 90, 13))
    lons = m.drawmeridians(np.linspace(-180, 180, 13))

    # keys contain the plt.Line2D instances
    lat_lines = chain(*(tup[1][0] for tup in lats.items()))
    lon_lines = chain(*(tup[1][0] for tup in lons.items()))
    all_lines = chain(lat_lines, lon_lines)
    
    # cycle through these lines and set the desired style
    for line in all_lines:
        line.set(linestyle='-', alpha=0.3, color='w')


def mask_observations(observations, inversion):
    new_inversion = np.copy(inversion)
    for i in range(len(observations)):
        for j in range(len(observations[0])):
            if np.isnan(observations[i][j]):
                new_inversion[i][j] = np.nan
    return new_inversion

'''
FROM DAILY DATA: 
fb_path1 = "../carmen/daily_numpys/AK_CPOM/FB_interp_2016-2017_25km_20170309.npy"
fb_path2 = "../carmen/daily_numpys/CS2_CPOM/FB_interp_2016-2017_25km_20170309.npy"

'''
def main(
    fb_path1 = "../carmen/daily_numpys/IS2/FB_interp_2018_25km_20181015.npy",
    fb_path2 = "../carmen/daily_numpys/CS2_LARM/FB_interp_2018-2019_25km_20181015.npy",
    verbose=False,
    minlat = -5000000.0/1000000, maxlat = 5000000.0/1000000,
    minlon = -5000000.0/1000000, maxlon = 5000000.0/1000000,
    parametrization = 1,
    initial_cells = 1,
    iterations_number = 250000,
    verbosity = 5000,
    independent_chains = 4,
    temperature_levels = 1,
    maximum_temperature = 2.0,
    iterations_between_tempering_attempts = 10,
    skipping = 100000,
    thinning = 5,
    render_observations=False,
         ):

    if verbose:
        print("Starting inversion")
    
    '''
    Step 1: Data cleaning and adapting to the TransTessellate standard
    '''
    is2 = np.load(fb_path1)
    cs2 = np.load(fb_path2)
    

    fb1_filename = os.path.split(fb_path1)[1]
    fb2_filename = os.path.split(fb_path2)[1]


    grid_x = np.load("../../new_x_25km.npy")
    grid_y = np.load("../../new_y_25km.npy")


    is2_observations = []
    cs2_observations = []

    for i in range (360):
        for j in range (360):
            if not np.isnan(is2[i][j]):
                is2_observations.append([grid_x[i][j]/1000000, grid_y[i][j]/1000000, 1, is2[i][j], 0.01])
            if not np.isnan(cs2[i][j]):
                cs2_observations.append([grid_x[i][j]/1000000, grid_y[i][j]/1000000, 0, cs2[i][j], 0.01])


    cs2_observations.extend(is2_observations)
    data = pd.DataFrame(cs2_observations, columns=["Longitude", "Latitude", "Type", "Value", "StdDev"])


    if render_observations:
        fig, ax = plt.subplots(1, 2, figsize=(15, 12))
        
        data[data['Type']==0].plot(kind='scatter', x='Longitude', y='Latitude', c='Value', cmap='seismic', ax=ax[0], title=f'{fb1_filename}')
        data[data['Type']==1].plot(kind='scatter', x='Longitude', y='Latitude', c='Value', cmap='seismic', ax=ax[1], title=f'{fb2_filename}')

        plt.show()
    
    observations_matrix_subset = data.values

    np.savetxt("observations.txt", observations_matrix_subset, '%5.1f %5.1f %d %5.5f %5.5f')

    # Add the total number of observations at the top of the file
    with open('observations.txt', 'r') as original: data = original.read()
    with open('observations.txt', 'w') as modified: modified.write(f"{observations_matrix_subset.shape[0]}\n" + data)

    '''
    Step 2: Performing inversion
    '''
    # Hyperparameters
    # number_of_processes = 1
    #parametrization = 0 # 0 for Voronoi, 1 for Delaunay linear, 2 for Delaunay Clough-Tocher
    #iterations_number = 100000
    #verbosity = 50000
    #
    # Run inversion
    subprocess.run([
                "mpirun", "-np", str(independent_chains * temperature_levels),
                "./snow_icept", 
                "-i", "observations.txt", 
                "-o", "results/", 
                "-P", "priors/prior_snow.txt",
                "-P", "priors/prior_ice.txt", 
                "-M", "priors/positionprior_snow.txt", 
                "-M", "priors/positionprior_ice.txt",
                "-H", "priors/hierarchical_snow.txt", 
                "-H", "priors/hierarchical_ice.txt", 
                "-C", str(initial_cells),
                "-x", str(minlon), "-X", str(maxlon),
                "-y", str(minlat), "-Y", str(maxlat),
                "-A", str(parametrization), "-A", str(parametrization),
                "-t", str(iterations_number), 
                "-v", str(verbosity),
                "-c", str(independent_chains),    # Independent chains to run at each temperature
                "-K", str(temperature_levels),    # Number of temperature levels for parallel tempering
                "-m", str(maximum_temperature),  # Maximum temperature for the parallel tempering log temperature
                "-e", str(iterations_between_tempering_attempts)    # Number of iterations between parallel tempering exchange attempts
                ])

    # Step 3: Compute means 
    parameter_W = 360
    parameter_H = 360

    file_snow = f"images/{fb1_filename}_snow"
    file_ice = f"images/{fb1_filename}_ice"


    subprocess.run([
                "mpirun", "-np", str(independent_chains),
                "./post_mean_mpi", "-i", 
                "results/ch.dat", "-o", file_snow,
                "-x", str(minlon), "-X", str(maxlon),
                "-y", str(minlat), "-Y", str(maxlat),
                "-s", str(skipping),
                "-t", str(thinning),
                "-A", str(parametrization), "-A", str(parametrization),
                "-W", str(parameter_W), "-H", str(parameter_H),
                "-D", str(file_snow + "_stddev"),

                "-I", str(0)])

    subprocess.run([
                "mpirun", "-np", str(independent_chains),
                "./post_mean_mpi", "-i", 
                "results/ch.dat", "-o", file_ice,
                "-x", str(minlon), "-X", str(maxlon),
                "-y", str(minlat), "-Y", str(maxlat),
                "-s", str(skipping),
                "-t", str(thinning),
                "-A", str(parametrization), "-A", str(parametrization),
                "-W", str(parameter_W), "-H", str(parameter_H),
                "-D", str(file_ice + "_stddev"),
                "-I", str(1)])            
       


if __name__ == "__main__":
    main()
