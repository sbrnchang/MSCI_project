'''
The objective of this program is to provide the entire pipeline from directly consuming CPOM 
data and directly producing plots with the results. 

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


def main(
    observation_points = 600,
    model_1 = "SyntheticSnow",
    model_2 = "SyntheticIce",
    noise = 0.01,
    verbose=False,
    minlat = -0.1, maxlat = 7950000.0/1000000,
    minlon = -0.1, maxlon = 7950000.0/1000000,
    parametrization = 2,
    iterations_number = 100000,
    verbosity = 5000,
    independent_chains = 4,
    temperature_levels = 1,
    maximum_temperature = 2.0,
    iterations_between_tempering_attempts = 10,
    skipping = 10000,
    thinning = 2,
    render_map=False,
    render_matrix=False,
    render_models = False,
    render_median = False,
    render_stddev = False,
    render_histogram = False
    ):

    if verbose:
        print("Starting inversion")
    

    if render_models:
        imgA = np.loadtxt("synthetic/syntheticobs_franke.img.A")
        imgB = np.loadtxt("synthetic/syntheticobs_franke.img.B")

        fig, ax = plt.subplots(1, 2, figsize=(15, 12))
        
        img = ax[0].imshow(imgA, cmap='seismic', aspect='auto', interpolation='None', origin='lower')
        ax[0].set_title('Synthetic snow depth (m)')
        plt.colorbar(img, ax=ax[0])

        img = ax[1].imshow(imgB, cmap='seismic', aspect='auto', interpolation='None', origin='lower')
        ax[1].set_title('Synthetic ice thickness (m)')
        plt.colorbar(img, ax=ax[1])

        plt.show()


    '''
    Step 2: Performing inversion
    '''
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
    parameter_W = 160
    parameter_H = 160

    file_snow = f"images/snow"
    file_ice = f"images/ice"


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
                "-m", str(file_snow + "_median"),
#                "-g", str(file_snow + "_histogram"),
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
                "-m", str(file_ice + "_median"),
#                "-g", str(file_ice + "_histogram"),
                "-I", str(1)])            
       
    '''
    Step 4: Produce and save plots
    '''
    snow_mat = np.loadtxt(file_snow)
    ice_mat = np.loadtxt(file_ice)


    lon = np.linspace(minlon, maxlon, 160)
    lat = np.linspace(minlat, maxlat, 160)
    lon_g, lat_g = np.meshgrid(lon, lat)

    lon_g = np.load("will_lons.npy")[:-1, :-1]
    lat_g = np.load("will_lats.npy")[:-1, :-1]

    extent = [minlon, maxlon, minlat, maxlat]




if __name__ == "__main__":
    main()
