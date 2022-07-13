[![Travis-CI](https://travis-ci.org/CICE-Consortium/CICE.svg?branch=master)](https://travis-ci.org/CICE-Consortium/CICE)
[![GHActions](https://github.com/CICE-Consortium/CICE/workflows/GHActions/badge.svg)](https://github.com/CICE-Consortium/CICE/actions)
[![Documentation Status](https://readthedocs.org/projects/cice-consortium-cice/badge/?version=master)](http://cice-consortium-cice.readthedocs.io/en/master/?badge=master)
[![lcov](https://img.shields.io/endpoint?url=https://apcraig.github.io/coverage.json)](https://apcraig.github.io)

<!--- [![codecov](https://codecov.io/gh/apcraig/Test_CICE_Icepack/branch/master/graph/badge.svg)](https://codecov.io/gh/apcraig/Test_CICE_Icepack) --->

## The CICE Consortium sea-ice model with WIM

CICE6 includes the capabilities for measuring the impacts of waves on sea ice across a global scale. In initial testing of CICE6 it was apparent that there is no code for the propagation and attenuation of wind waves and swell through sea ice. Without waves no ice breakup will occur and all new floes are assumed to have the geometry of pancake ice, rather than using the theory of Shen et al. (2001). The Waves-in-Ice Module (WIM) module is a way to propagate waves through the ice pack in CICE6, such that preliminary runs can be completed using CICE's new FSD routine. This module is from the code used in Bennetts et al. (2017) ([https://bitbucket.org/puotila/cicewithwaves/src/master/](https://bitbucket.org/puotila/cicewithwaves/src/master/)), although some modifications have been made to make it compatible with the latest version of CICE. It was decided that CICE's ice breakup routine (Roach et al., 2018; Horvat et al., 2015) would be used within this model to keep comparisons as simple as possible. The short-coming of this decision is that wave-sea ice interactions are not coupled, we are currently investigating the effects of this.

The code primarily aims to achieve: 

1. The definition of the edge of the ice pack (determined when the ice concentration < puny)
2. Read in data for these respective cells from a WW3 hindcast dataset ([https://data.csiro.au/collections/collection/CI6616v010](https://data.csiro.au/collections/collection/CI6616v010)).
3. Use a Bretscheider spectrum to define the wave spectrum for each cell along the edge of the wave mask.
4. Propagate waves given the significant wave height, peak period, and mean wave direction into the ice pack.
5. Feed this data into the icepack thermodynamic routine.

Attenuation is computed as per the observation of Meylan et al. (2014) currently, there is the code in place to implement the attenuation model presented in Williams et al. (2013a,2013b), but this has not been tested within CICE6 yet.
    
Method:
1. Define wave mask edge (from areal ice concentration) 
2. Define the Bretschneider spectrum along the mask edge
3. Calculate the attenuation coefficients from Meylan et al. (2014)
4. Propagate waves through the sea ice
5. Calculate the peak period and signficant wave height for each cell 
 
For more details see https://noahday.notion.site/README-af9c753fb9ab47f2a0220bda3cfa256d.        

## Dependencies

The code is initiated within `ice_step_mod`, in the subroutine `step_therm2` 

### `CICE_RUN`

### `ice_forcing`

### `ice_step_mod`

## Reference list:

Bennetts, L., O’Farrell, S. and Uotila, P. (2017) ‘Brief communication: impacts of ocean-wave-induced breakup of Antarctic sea ice via thermodynamics in a stand-alone version of the CICE sea-ice model’, The Cryosphere

Horvat, C. and Tziperman, E. (2015) ‘A prognostic model of the sea-ice floe size and thickness distribution’, The Cryosphere, 9(6), pp. 2119–2134. doi:10.5194/tc-9-2119-2015.

Meylan, M.H., Bennetts, L.G. and Kohout, A.L. (2014) ‘In situ measurements and analysis of ocean waves in the Antarctic marginal ice zone.’, Geophysical Research Letters, 41(14), pp. 5046–5051. doi:10.1002/2014GL060809.)

Roach, L.A. et al. (2018) ‘An Emergent Sea Ice Floe Size Distribution in a Global Coupled Ocean-Sea Ice Model’, Journal of Geophysical Research: Oceans, 123(6), pp. 4322–4337. doi:10.1029/2017JC013692.. doi:10.5194/tc-11-1035-2017.

Shen, H.H., Ackley, S.F., Hopkins, M.A., 2001. A conceptual model for pancake-ice formation in a wave field. Ann. Glaciol. 33, 361–367. https://doi.org/10.3189/172756401781818239

Williams, T.D. et al. (2013) ‘Wave–ice interactions in the marginal ice zone. Part 1: Theoretical foundations’, Ocean Modelling, 71, pp. 81–91. doi:10.1016/j.ocemod.2013.05.010.

Williams, T.D. et al. (2013) ‘Wave–ice interactions in the marginal ice zone. Part 2: Numerical implementation and sensitivity studies along 1D transects of the ocean surface’, Ocean Modelling, 71, pp. 92–101. doi:10.1016/j.ocemod.2013.05.011.

## ❓FAQ

- What dimension is the WIM model?

    Currently, the WIM is limited to 1D (along longitudinal lines). There is some infrastructure for moving to a 2D wave propagation scheme, although Alberto Alberello and Luke Bennetts ran into problems regarding the transfer of wave energy between blocks in CICE.
    
- What is the wave spectrum made up of?
    
    Using significant wave height and a peak period a Bretschneider spectrum is defined. This spectrum is then propagated either North or South according to it's mean wave direction.
    
- Where are waves propagated from?
    
    Yes, from the closest ice-free cell. It uses a 1D wave scattering advection scheme, so it reads in wave direction but theses directions don't evolve over time.
    
- What output does this module produce?
    
    Wave spectrum, peak period, significant wave height for each cell.
    
- What FSD do they use?
    
    They use the FSD routine in CICE6 (Roach et al., 2017). However it could be possible to implement a split power law to describe the FSD following Toyota et al. 2017 (as seen in Bennetts et al., 2017) for the breakup routine.
    
- Does it use wave direction?
    
    They do use wave direction from the closest ice free cell. The wave spectrum is 'advected' through the ice. I say 'advected' as it is really only attenuation rather than full advection. And the wave direction does not change as the waves collide with floes and other waves.
    
- Relevant variable dimensions
    
    ```fortran
    real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
             ifd     , & ! ice floe diameter (m)
             swh     , & ! significant wave height (m)
             mwd     , & ! mean wave direction (Rads)
             ppd         ! wave peak period (s)
    ```
    
- What are the attenuation settings?
    
    There are two kinds of attenuation present in this module, the simplest is based off of observations from Meylan et al. (2014). The more complex attenuation model has not currently been tested in CICE6 and is based off of the work by Williams et al. (2013a,2013b).
    

### Extra techinical details
    
`increment_floe` is the main script in the module, it calls all other files in the `wave_ice_code` and exists within CICE code. This code first initalises: frequency min/max, angular frequencies, frequency, wave direction, initial energy spec, wavelength, wave number, dummy spectral moments, wave spectrum tolerance, spectrum passed through ice-covered ocean, sine mean wave direction for rows, values required for calculations of mean wave direction, WIM termination flags, attenuation parameters, and points in frequency and angular domains.
    
This code proceeds by defining the wave mask edge and propagate waves into ice, define attenuation coefficient coefficients. Then it calculates the wavelengths and wave number by `lambda = gT^2/2pi`, `kappa = 2pi/lambda` as well as direction `th_in = -pi/2 + (lp_i-1)pi/(nth_in-1)`. Then the Bretschneider SDF function is called to initialise the wave spectrum. 
    
This gets the code to the point where it can update the wave spectrum and the floe sizes. This can work in both the coupled and uncoupled set ups by calling `sub_Balance` or `sub_Uncoupled` respectively. `sub_Balance` first calculates the `mass,`and  `fr`, where `mass` is the scaled mass and `fr` is the flexural ridging (note that `dom` and `dth` are `d*omega` and `d(th)`). It then implements [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule) for `dth` and calculates `lambda_ice`. Its advection scheme differs from `sub_Uncoupled` as it uses the functions `fn_IntAttn` and `fn_AvAttn` to calculate the integrated attenuation coefficient with respect to floe diameter and the averaged attenuation coefficient with respect to floe diameter respectively. This calculation takes the maximum floe size, `lambda_ice(lp_i)` and `om(lp_i)` as input.  The wave attenuation is handled in the same exponential (+ cosine) scheme as in the previous subroutine (although with different coefficient values).
