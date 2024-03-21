# TransTesselate2D for snow depth and sea ice thickness retrievals
This code is the basis for the MSci Physics Project.

Original code is by: Hawkins R., Bodin T., Sambridge M., Choblet G. and Husson L., "Trans-dimensional surface reconstruction with different classes of parameterization", Geochemistry, Geophysics, Geosystems, 2019

See documentation/manual.tex Latex file for the original tutorial/manual on how to customize the code.

# Compilation

This software is written in C++ and requires

- GNU g++ Version 6.x or greater
- GNU fortran Version 6.x or greater
- GNU Make version 4.x or greater
- GNU Scientific Library (GSL) version 1 or 2
- OpenMPI version 1.10

In the root of the directory:
```
> make 
```

## Snow and ice retrieval

Inside of the `snow_ice` folder there are a few models that can be used. In order to run them they first need to be compiled calling:
```
> make
```

If error during make try 
``` 
> make clean
> make
```
Each model is set with Voronoi cell, for 4 chains. 2p models run for roughly one hour. 

To run the inversion:
```
> python InversionBinnedParallel.py
```

After inversion:
Validation plots with Operation Ice Bridge are in Jupyter Notebooks of each model labelled `OIB.ipynb`
Mapping of inversion plots are in each Jupyter Notebooks of each model `mapping.ipynb`

## Synthetic simulations

In `synthetic`:

`acceptance_tests` and `initial_prior_test` are used for simple synthetic case to test how sensitive the software is with different parameters.

`realistic_simulation` is for synthetic case at the Arctic.

in each folder, run tests in `SyntheticTest.ipynb` and analysis, including acceptance rate in `Analysis.ipynb`.


