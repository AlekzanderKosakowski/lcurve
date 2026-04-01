Modified copy of Tom Marsh's source code used for "lcurve". For the full set of original files, see the source repository:

    https://github.com/trmrsh


Modifications are described below.


### 1) "visualise" colors:
visualise.cc now plots colors based on surface element temperature, normalized by the [min_temperature, max_temperature] of all of the surface elements across star1 and star2 together.

This has only been tested for stellar components, not for accretion discs or accretion stream-disc impact spots.

### 2) Starspot irradiation:
The original code treats star1 as a point source for irradiation onto star2. This worked well, but meant ignoring contribution from star1's starspots. I've adjusted the source code to allow the user to enable a flag (finite_irr12) in their parameters.txt file to include contribution from all of star1's surface elements as irradiation contributors towards star2.

This was done via a nested loop, so the runtime increases significantly when enabled (I observed a runtime increase from roughly 0.5s to 1.5s per calculation). Try to keep the stellar grid resolution parameters (nlat1c, nlat2c, nlat1f, nlat2f) low since set_star_continuum() runtime grows as $\mathcal{O}(N^2)$ with this flag enabled.

![Starspot irradiation at latitude 0 degrees](figures/starspot_irradiation_0deg.gif)
![Starspot irradiation at latitude 60 degrees](figures/starspot_irradiation_60deg.gif)


