# How to run this code:

## About

This code runs a series of models of increasing complexity. The models are fitting a distribution, or distributions, to a set of body size data. The models start but fitting to simulated data, as a test of concept, before fitting to real data. The models are also fit to a subset of the species before fitting to all species.

The file structure must remain the same in order for the code to run (i.e. raw data remaining in the `input/data/` folder).

## Input

Observed fish body size data can be found in the `input/data/` folder. The stan code for fitting the bayesian models can be found in the `input/models/` folder.

## Output

# To-do list

-   run the population models for longer, with tighter priors on eps parameters (4000 runs, with mean of eps at -6 and sigma of 0.5, upper limit of -2)

    -   **158pops**

    -   1223pops

-   run for ecoregions

    -   188ecos

    -   822ecos
