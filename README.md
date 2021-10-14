[![Travis-CI](https://travis-ci.org/CICE-Consortium/CICE.svg?branch=master)](https://travis-ci.org/CICE-Consortium/CICE)
[![GHActions](https://github.com/CICE-Consortium/CICE/workflows/GHActions/badge.svg)](https://github.com/CICE-Consortium/CICE/actions)
[![Documentation Status](https://readthedocs.org/projects/cice-consortium-cice/badge/?version=master)](http://cice-consortium-cice.readthedocs.io/en/master/?badge=master)
[![lcov](https://img.shields.io/endpoint?url=https://apcraig.github.io/coverage.json)](https://apcraig.github.io)

<!--- [![codecov](https://codecov.io/gh/apcraig/Test_CICE_Icepack/branch/master/graph/badge.svg)](https://codecov.io/gh/apcraig/Test_CICE_Icepack) --->

## The CICE Consortium sea-ice model
# README

The Waves-in-Ice Module (WIM) module is a way to propagate waves through the ice pack in CICE6, such that preliminary runs can be completed using CICE's new FSD routine. This module is from the code used in Bennetts et al. (2017) ([https://bitbucket.org/puotila/cicewithwaves/src/master/](https://bitbucket.org/puotila/cicewithwaves/src/master/)), although some modifications have been made to make it compatible with the latest version of CICE. It was decided that CICE's ice breakup routine (Roach et al., 2018; Horvat et al., 2015 check!!) would be used within this model to keep comparisons as simple as possible. 

The code primarily aims to achieve: 

1. The definition of the edge of the ice pack (determined when the ice concentration < puny)
2. Read in data for these respective cells from a WW3 hindcast dataset ([https://data.csiro.au/collections/collection/CI6616v010](https://data.csiro.au/collections/collection/CI6616v010)).
3. Use a Bretscheider spectrum to define the wave spectrum for each cell along the edge of the wave mask.
4. Propagate waves given the significant wave height, peak period, and mean wave direction into the ice pack.
5. Feed this data into the icepack thermodynamic routine.

Attenuation is computed as per the observation of Meylan et al. (2014) currently, there is the code in place to implement the attenuation model presented in Williams et al. (2013a,2013b), but this has not been tested within CICE6 yet.

# WIM

- Explanation of code
    
    `increment_floe` is the main script in the module, it calls all other files in the `wave_ice_code` and exists within CICE code. This code first initalises: frequency min/max, angular frequencies, frequency, wave direction, initial energy spec, wavelength, wave number, dummy spectral moments, wave spectrum tolerance, spectrum passed through ice-covered ocean, sine mean wave direction for rows, values required for calculations of mean wave direction, WIM termination flags, attenuation parameters, and points in frequency and angular domains.
    
    This code proceeds by defining the wave mask edge and propagate waves into ice, define attenuation coefficient coefficients. Then it calculates the wavelengths and wave number by $\lambda = \frac{gT^2}{2\pi}$, $\kappa = \frac{2\pi}{\lambda}$ as well as direction $th_{in} = -\frac{\pi}{2} + \frac{(lp_{i}-1)\pi}{nth_{in}-1}$. Then the Bretschneider SDF function is called to initialise the wave spectrum. 
    
    This gets the code to the point where it can update the wave spectrum and the floe sizes. This can work in both the coupled and uncoupled set ups by calling `sub_Balance` or `sub_Uncoupled` respectively. `sub_Balance` first calculates the `mass,`and  `fr`, where `mass` is the scaled mass and `fr` is the flexural ridging (note that `dom` and `dth` are $d\omega$ and $d(th)$). It then implements [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule) for `dth` and calculates $\lambda_{ice}$. Its advection scheme differs from `sub_Uncoupled` as it uses the functions `fn_IntAttn` and `fn_AvAttn` to calculate the integrated attenuation coefficient with respect to floe diameter and the averaged attenuation coefficient with respect to floe diameter respectively. This calculation takes the maximum floe size, $\lambda_{ice}(lp\_i)$ and $om(lp\_i)$ as input.  The wave attenuation is handled in the same exponential (+ cosine) scheme as in the previous subroutine (although with different coefficient values).
    

        Method:
        
        1. Define wave mask edge and propagate waves into ice
        2. Define attenuation coefficient coefficients
        3. Calculate wave numbers and wavelengths
            
 
        
- Summary of `wave_ice_code`
    
    [m_waveice documentation](https://www.notion.so/m_waveice-documentation-0938744c39e9488b9e01d2ac499a44e2)
    
    There are two attenuation routines.
    
    1. `fn_Attn_MBK` uses the attenuation model derived from [Meylan et al. 2014](https://www.notion.so/Meylan-et-al-2014-2c785c506c784e4c85fa00270fff63ff) (Geophys Res Lett). Which is called in `sub_Uncoupled` and does not depend on the floe size.
    2. `fn_Attn_WIM_v1` uses the attenuation model from [Williams et al. 2013a,b](https://www.notion.so/Williams-et-al-2013a-b-5eae5120f7164a1ba164154ec1189e9f) (Ocean Model.) Which is called in `sub_Balance` and does depend on floe size
    
    There is also a [Chebyshev interpolation](https://www.notion.so/Chebyshev-interpolation-1ebf5c9254504928a367e0d391175f5d) function (in 2D)
    
    `lp_i` runs from `1,...,nw` uses a Bretschneider spectrum.
    
    - `sub_Uncoupled`, attenuation does not depend on floe size
        
        This subroutine uses the Meylan et al. 2014 attenuation coefficient scheme. The floe size is not accounted for in the attenuation. Interestingly, there is a slight difference in a coefficient between *just wave advection* and *potential for fracture* cases (1/0.7 and 1/0.75 respectively), the attenuation coefficient is also a product of ice concentration. The wave spectrum takes an exponential attenuation.
        
        [sub_Uncoupled](https://www.notion.so/d2b5b8129794403886474f8c879342d3)
        
        L: MIZ length
        
        - More details
            
            ### Nothing but wave advection, line 288
            
            1. Calculation of attenuation coefficient, $\alpha_{lp\_i} = \text{conc}*\text{fn\_Attn\_MBK}(\text{om(lp\_i)}\frac{1}{0.7})$
            2. Calculation of attenuated wave spectrum, 
            
            $$S_{\text{attn}}(\text{lp\_i+nw*(lp\_j-1)}) = S_{\text{init}}(\text{lp\_i+nw*(lp\_j-1)})*\exp(-\alpha_{lp\_i}*\cos(\text{th(lp\_j))*L})$$
            
            ### Potential for fracture, line 312
            
            1. Calculation of attenuation coefficient, $\alpha_{lp\_i} = \text{conc}*\text{fn\_Attn\_MBK}(\text{om(lp\_i)}\frac{1}{0.75}$
            2. Calculation of attenuated wave spectrum, 
            
            $$S_{\text{attn}}(\text{lp\_i+nw*(lp\_j-1)}) = S_{\text{init}}(\text{lp\_i+nw*(lp\_j-1)})*\exp(-\alpha_{lp\_i}*\cos(\text{th(lp\_j))*L})$$
            
            1. call `sub_StrainSpec(S_attn,Es)`, `Es` is the strain
            2. ! Set attenuated spectrum (as not floe size dependent)
            
            ### Floes already smaller than breaking size
            
            1. `floe_sz_max = floe_sz_init`
            
            ### Waves strong enough to break ice at end of cell
            
            1. `floe_sz_max = floe_sz_emp`
            
            else search for a point at which waves can break ice
            
            ```bash
            Lc = zero(i0,Lcell,toli,toli,fn_proplng_uncoupled)
            L           = Lc
            floe_sz_max = (Lc*floe_sz_brk_emp + (Lcell-Lc)*floe_sz_init)/Lcell
            tmt         = 1   ! waves attenuated (in all liklihood)
            ```
            
    - `sub_Balance`,  no description provided
        
        This subroutine first calculates the `mass`, `fr`, and `dom`, where `mass` is the scaled mass and `fr` is the flexural ridging. It then implements [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule) for `dth` and calculates $\lambda_{ice}$. Its advection scheme differs from `sub_Uncoupled` as it uses the functions `fn_IntAttn` and `fn_AvAttn` to calculate the attenuation coefficient $\alpha$. This calculation takes the maximum floe size, $\lambda_{ice}(lp\_i)$ and $om(lp\_i)$ as input. But the wave attenuation is handled in the same exponential scheme as in the previous subroutine (although with different coefficient values).
        
    - `fn_AvAttn`, the averaged attenuation coefficient with respect to floe diameter
        
        This function uses the [Williams et al. 2013a,b](https://www.notion.so/Williams-et-al-2013a-b-5eae5120f7164a1ba164154ec1189e9f) formulation for attenuation.
        
        1. Find the average floe size
        2. Calculate the attenuation coefficient, $\alpha$
        3. Set the attenuation to be $= \frac{\text{concentration}\times \alpha}{\text{Average floe size}}$ 
        
        From Williams et al. (2013a), 
        
        $$\hat{\alpha} = \frac{\alpha c}{\lang D \rang}$$
        
        where $\alpha$ is the non-dimensional attenuation coefficient, i.e. the (average) amount of attenuation per individual floe, which is a function of ice thickness and wave period.
        
        ```fortran
        ! D average
         
         Dav = gamma*(floe_sz_max*((floe_sz_min/floe_sz_max)**gamma)-floe_sz_min)/(1d0-gamma)
         
         ! non-dim attenuation
         
         if (ATTEN_MODEL.eq.0) then
          alpha_nd = fn_apxAttn(Dav,lam,1d0)
         elseif (ATTEN_MODEL.eq.1) then 
          alpha_nd = fn_Attn_WIM_v1(dum_om,hice,dum_idl)
         end if
         
         ! Dimensionalise
         
         fn_AvAttn = conc*alpha_nd/Dav
        ```
        
    - `fn_IntAttn`, the integrated attenuation coefficient with respect to floe diameter
        
        This function find the midpoints of the floe size (`dD`), calculates the FSD parameters $\beta_0$ and $\beta_1$  (see [Bennetts et al. 2017](https://www.notion.so/Bennetts-et-al-2017-b3546893eb5844a3a4b42f5b1b3c664c)). From which, it calculated $\mathbb{P}_0$ and then goes on to calculate $\alpha\_vec$ according to whether the code specified for `fn_apxAttn` or `fn_Attn_WIM_v1` attenuation models. $p\_vec$ is also calculated using `fn_splitPDF`. Finally, `fn_IntAttn` is set to `dD*conc*dot_product(alpha_vec,p_vec/D_vec)`, for which it is outputted.
        
        From Bennetts et al. (2017)
        
        $$\beta_0 = \frac{1}{D_{mn}^{-\gamma_0} - D_{cr}^{-\gamma_0}}$$
        
        $$\beta_1 = D_{cr}^{\gamma_1}$$
        
        $$\mathbb{P}_0 = 1-q \left( \frac{D_{pr}}{D_{cr}} \right) ^{\gamma_1}$$
        
        where $q = 0.05$, and $D_{pr} = \lambda/2$.
        
   
    - `fn_proplng`, no description provided
        
        From Williams et al. (2013a):
        
        Equation 14
        
        The criterion $\mathbb{P}_\epsilon > \mathbb{P}_c$ can be written as
        
        $$E_s > E_c = \epsilon_c \sqrt{-2/\log{\mathbb{P}_c)}}$$
        
        ```fortran
        L = dum_L
         
         dum_D1 = zero(i0,i1,toli,toli,fn_wlng)
         
         do loop_w=1,nw
          if (ATTEN_METH.eq.1) then
           alpha(loop_w) = fn_IntAttn(dum_D1, lam_ice(loop_w), om(loop_w))
          else
           alpha(loop_w) = fn_AvAttn(dum_D1, lam_ice(loop_w), om(loop_w))
          end if
          S_attn(loop_w) = S_init(loop_w)*exp(-alpha(loop_w)*L)
         end do 
         
         call sub_StrainSpec(S_attn, dum_Es)
         
         fn_proplng = dum_Es - epsc*sqrt(-2/log(Pc))
        ```
        
    - `fn_proplng_uncoupled`, `fn_proplng` for problem when attn/floe length are uncoupled.
        
        Uses the same method as `fn_proplng` just with the Meylan et al. (2014) attenuation rather than Williams et al. (2013a).
        
        ```fortran
        L = dum_L
         
         do loop_w=1,nw
          alpha(loop_w)  = conc*fn_Attn_MBK(om(loop_w))/0.7d0 
          S_attn(loop_w) = S_init(loop_w)*exp(-alpha(loop_w)*L)
         end do 
         
         call sub_StrainSpec(S_attn, dum_Es)
         
         fn_proplng_uncoupled = dum_Es - epsc*sqrt(-2/log(Pc))
        ```
        
    - `fn_SpecMoment`, no description provided
        
        
        Set `dum_simp` to constants
        
        ```fortran
        dom_local = om_in(2)-om_in(1)
         
         if (nth_in.eq.1) then
          dth_local = 3d0         ! for Simpson's rule
         else 
          dth_local = th_in(2)-th_in(1)
         end if
        
        do loop_w=1,nw_in
          do loop_th=1,nth_in
           wt_int(loop_w+nw_in*(loop_th-1)) = (dom_local/3d0)*(dth_local/3d0)*wt_simp(loop_w+nw_in*(loop_th-1))
           dum_v(loop_w+nw_in*(loop_th-1))  = dum_S(loop_w+nw_in*(loop_th-1))*(om_in(loop_w)**mom)
          end do
         end do
        
        fn_SpecMoment   = dot_product(wt_int,dum_v)
        ```
        
        $$\text{fn\_SpecMoment} = \text{wt\_int} \cdot \text{dum\_v} $$
        
    - `sub_WavelenSpec`, calculates a representative wavelength for the spectrum
        
        `kappa` := frequency parameter. 
        
        See [wave propagation](https://www.notion.so/Sea-Ice-Wiki-4a4ee6b89ca74373a7792d6bfa17c42d).
        
        The gravity wave angular frequency can be expressed as $\omega = \sqrt{gk}$ (using shallow water relation, [source](https://en.wikipedia.org/wiki/Gravity_wave)), which can be rearranged to $k = \frac{\omega^2}{g}$.
        
        `om_crit` is calculated following [this](https://www.notion.so/Theory-Calculation-of-significant-wave-height-in-the-ice-pack-6191341a58e442b885f59981633fa986),  $\omega_{crit} = 1/T_{02} = \sqrt{\frac{m_2}{m_0}}$, where $T_{02}$ is the mean wave period ([source](http://www.coastalwiki.org/wiki/Statistical_description_of_wave_parameters)).
        
        `wlng_crest` := representative wavelength for the spectrum.
        
        ```fortran
        mom0     = dot_product(wt_int,dum_u) 
        mom2     = dot_product(wt_int,dum_v)
        om_crit = sqrt(mom2/mom0)
        kappa = om_crit**2d0/gravity
        wlng_crest = 2d0*pi/zero(0d0,max(kappa/tanh(kappa),sqrt(sqrt(kappa*mass/fr))),toli,toli,fn_DispRel_ice_inf)
        ```
        
    - `sub_StrainSpec`, calculates the strain imposed on the ice by the wave spectrum
        
        From Williams et al. (2013a) The mean square value for the strain is $\lang \epsilon^2 \rang = m_0[\epsilon]$ where
        
        $m_0[\epsilon] = \int_0^\infty S(\omega) E^2(\omega)d\omega , E(\omega) = \frac{h}{2}k_{ice}^2W(\omega)$ 
        
        and the significant strain amplitude is defined as $E_s = 2\sqrt{m_o[\epsilon]}$, which is two standard deviatios in strain.
        
        ```fortran
        mom0_eps = dot_product(wt_int,dum_vec)
        
         Es = 2d0*sqrt(mom0_eps)
        ```
        
    - `fn_DispRel_ice_inf`, infinite depth ice-coupled dispersion relation
        
        `mass` := scaled mass
        
        `kappa` := frequency parameter
        
        `fr` := flexural ridging
        
        `dum_k` := 
        
        ```fortran
        fn_DispRel_ice_inf = (1d0-(mass*kappa)+(fr*(dum_k**4d0)))*dum_k - kappa
        ```
        
        $$dispRel = (1-mass\times \kappa + fr\times dum\_k^4)dum\_k - \kappa$$
        
    - `fn_wlng`, calculates wavelength (`fn_wlng = lam/2 - D1`)
        
        
        Input: `dum_D1`
        
        Output: `fn_wlng` = `lam_attn/2 - dum_D1`
        
        `fn_wlg` := wavelength
        
        ```fortran
        do loop_w=1,nw
          if (ATTEN_METH.eq.1) then
           alpha(loop_w) = fn_IntAttn(dum_D1, lam_ice(loop_w), om(loop_w))
          else
           alpha(loop_w) = fn_AvAttn(dum_D1, lam_ice(loop_w), om(loop_w))
          end if
          do loop_th=1,nth
           S_attn(loop_w+nw*(loop_th-1)) = S_init(loop_w+nw*(loop_th-1))* &
             exp(-alpha(loop_w)*cos(th(loop_th))*L)
          end do 
         end do 
        
         call sub_WavelenSpec(S_attn, lam_attn)
         
         fn_wlng = lam_attn/2d0 - dum_D1
        ```
        
    - `SDF_Bretschneider`, [Bretschneider wave spectrum](https://www.notion.so/USNA-waves-summary-significant-wave-height-wave-spectrum-2719a3d538c241fab7b8ab23af8a56ab)
        
        Input: $\omega, moment\_no, H_s, T_m$ ($H_s$ is significant wave height and $T_m$ is peak period.
        
        $$om\_m = \frac{2\pi}{T_m}$$
        
        $$\tau = \frac{2\pi}{\omega}$$
        
        $$Bretschneider = \frac{5}{16}H_s^2om\_m^4 \omega^{moment\_no-5}e^{-1.245(\tau/T_m)^4}$$
        
        $$S_\zeta(\omega) = \frac{1.25}{4} \left( \frac{\omega_0}{\omega} \right)^4 \frac{\bar{H}^2_{1/3}}{\omega} e^{-1.25(\omega_0/\omega)^4}$$
        
        Note that $om\_m = \omega_0$ as $T_m = 2\pi/\omega_0$, so these equations are equivalent.
        
        This includes another way of writing the Bretschneider spectrum: 
        
        [](https://ww2.eagle.org/content/dam/eagle/rules-and-guides/current/offshore/238_Guidance_Notes_on_Selecting_Design_Wave_by_Long_Term_Stochastic_Method/Long_Term_Design_Wave_GN_e.pdf)
        
        From [USNA waves summary, significant wave height, wave spectrum](https://www.notion.so/USNA-waves-summary-significant-wave-height-wave-spectrum-2719a3d538c241fab7b8ab23af8a56ab) $\omega_0 = \frac{w\pi}{T_0}$ where $T_0$ is the modal wave period.
        
- Summary `m_waveattn`
    
    Author: Luke Bennetts
    
    ```fortran
    use m_prams_waveice
    ```
    
    - `fn_Attn_MBK`, lines 59-85
        
        The attenuation model derived by Meylan et al. (2014), Geophys Res Lett
        
        $$\alpha = \beta_0 + \beta_1 \times om^2$$
        
        $$fn\_Attn\_MBK = attn\_fac \times( \beta_0(dum\_om^2) + \beta_1(dum\_om^4))$$
        
        From Meylan et al. (2014)
        
        The attenuation coefficient and a nonlinear fit of the median values are of the form
        
        $$\alpha(T) = \frac{a}{T^2} + \frac{b}{T^4}$$
        
        with $a = 2.12 \times 10^{-3} (s^2/m)$ and $b = 4.59\times 10^{-2} (s^4/m)$, which were found using least squares regression, and $T$ denotes wave period.
        
        In `m_prams_waveice.F90` line 61:
        
        ```bash
        attn_fac=1.0d0 ! adjust attenuation rate
        ```
        
        Line 72:
        
        ```fortran
        real(kind=8) :: fn_Attn_MBK
        
        !
        !EOP
        !
         
         real(kind=8), parameter :: beta0 = 5.376168295200780E-005, &
             beta1 = 2.947870279251530E-005
                 
         fn_Attn_MBK = beta0*(dum_om**2) + beta1*(dum_om**4)
         
         fn_Attn_MBK = attn_fac*fn_Attn_MBK
        ```
        
    - `fn_Attn_WIM_v1`, lines 91-185
        
        The attenuation model used for Williams et al (2013a,b) Ocean Model.
        
        Input: `dum_om`, `dum_h`, `dum_idl`
        
        `dum_om` := angular frequency
        
        `dum_h` := thickness
        
        Output: `fn_Attn_WIM_v1`
        
        `Ncheb_f` := degree of polynomial in period/frequency (set to 26)
        
        `Ncheb_h` := degree of polynomial in thickness (set to 26)
        
        Line 162:
        
        ```fortran
        ! resetting the minumum/maximum ice thickness, and the current ice thickness
        if (dum_h.lt.hmin) then
          new_h=hmin 
         elseif (dum_h.gt.hmax) then 
          new_h=hmax 
         else
          new_h=dum_h
         end if
         
         !print*, 'hnew=',new_h
         
         Limits=(/ xmin,xmax,hmin,hmax /) ! limits of the cell
         
         !print*, 'Limits',Limits
         
         fn_Attn_WIM_v1=fn_OP_chebinterp2d(dum_om,new_h,Ncheb_f,Ncheb_h,alp_coeffs,Limits)
          
         fn_Attn_WIM_v1 = attn_fac*fn_Attn_WIM_v1
        ```
        
        From Williams et al. (2013b)
        
        > We note that the quantities provided by the wave and sea ice models to the WIM are likely to require interpolation onto the high resolution grid it uses.
        > 
    - `fn_OP_chebinterp2d`, lines 185-258
        
        Chebyshev interpolation of a function of 2 variables.
        
        [Resources](https://www.notion.so/Chebyshev-interpolation-1ebf5c9254504928a367e0d391175f5d).
        
        ```fortran
        dx       = x_max - x_min
        tt1      = -1d0 + 2d0*(xx-x_min)/dx
        TnVals1  = fn_OP_chebinterp1d(tt1,Ncheb1)
        
        dy       = y_max - y_min 
        tt2      = -1d0 + 2d0*(yy-y_min)/dy
        TnVals2  = fn_OP_chebinterp1d(tt2,Ncheb2)
        
        fn_OP_chebinterp2d = dot_product(TnVals1,matmul(F_coeffs,TnVals2))
        ```
        
    - `fn_OP_chebinterp1d`, lines 261-326
        
        Chebyshev interpolation of a function of 1 variable.
        
        Input:
        
        `tt` := 
        
        `An` := 
        
        ```fortran
        hn(1) = pi
        do lp=1,An
          hn(lp+1)=pi/2
        end do
        
        C0=1d0
        C1=tt 
        f(1)=C0
        f(2)=C1
        
        do lp=2,An
         Cn      = 2*tt*C1-C0
         f(lp+1) = Cn
         C0      = C1
         C1      = Cn
        end do
        
        fn_OP_chebinterp1d = f
        ```
        
- Wave Watch III data reading
    
    In subroutine `CICE_Run`
    
    Line 136:
    
    ```fortran
    !--------------------------------------------------------------------
    !  wave forcing: dummy at the moment
    !--------------------------------------------------------------------
          nmth=1
          if (WAVE_METH.eq.1) then
             call sub_WW3_dataread(nmth,N_tm,N_lat,N_lon)
             allocate(ww3_swh(N_lon,N_lat))
             allocate(ww3_fp(N_lon,N_lat))
             allocate(ww3_dir(N_lon,N_lat))
          endif
                  
          nww3=1-nww3_dt
          
          timeLoop: do
           
             !!!!! LUKE'S WAVE STUFF !!!!!
           
             nww3=nww3+nww3_dt
             if (WAVE_METH.eq.1) then
                print*, 'LB: istep,nww3,N_tm=', istep, nww3, N_tm
                if (nww3.le.N_tm) then
                    ww3_swh(:,:) = ww3_swh_full(:,:,nww3)
                    ww3_fp(:,:)  = ww3_fp_full(:,:,nww3)
                    ww3_dir(:,:) = ww3_dir_full(:,:,nww3)
                else
                    ww3_swh(:,:) = ww3_swh_full(:,:,N_tm)
                    ww3_fp(:,:)  = ww3_fp_full(:,:,N_tm)
                    ww3_dir(:,:) = ww3_dir_full(:,:,N_tm)
                    deallocate(ww3_swh_full)
                    deallocate(ww3_fp_full)
                    deallocate(ww3_dir_full)
                    deallocate(ww3_lat)
                    deallocate(ww3_lon)
                    deallocate(ww3_tm)
                endif
             endif
             nmth = nmth+1
             !print*, '    month=', nmth
             if (WAVE_METH.eq.1) call sub_WW3_dataread(nmth,N_tm,N_lat,N_lon)
             nww3=1-nww3_dt
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ```
    
    - `sub_WW3_dataread`
        
        This checks if the data has been read in correctly, the function `check` is found in `ice_forcing`. This code expects the data comes in monthly netcdf files.
        
        Input: `mth`, month 1-12
        
        Output: `N_tm`, `N_lat`, `N_lon`
        
        ```fortran
        if (mth.eq.1) then
               call check( nf90_open(waveicedatadir // fname_ww3 // '01_full.nc', NF90_NOWRITE, ncid) )
              elseif (mth.eq.2) then
        ...
        ```
        
        - `waveicedatadir`
            
            In `m_prams_waveice`
            
            ```fortran
            ! WAVES
            
             integer, parameter                    :: WAVE_METH = 0   ! inc waves (0=user,1=ww3)
             real(kind=8), allocatable             :: ww3_lat(:,:), ww3_lon(:,:), ww3_tm(:,:)
             real(kind=8), allocatable             :: ww3_swh(:,:), ww3_fp(:,:), ww3_dir(:,:), &
             				ww3_swh_full(:,:,:), ww3_fp_full(:,:,:), ww3_dir_full(:,:,:)
             integer, parameter                    :: nww3_dt = 2 ! ww3 time step relative to cice
            
            waveicedatadir='../../waveice_data/'
            character(10)            :: fname_alp='alp_coeffs'
             !character(24)            :: fname_ww3='waves/ww3.197803_full.nc'
             character(14)            :: fname_ww3='waves/ww3.1978'
            integer, parameter       :: OVERWRITE_DIRS = 0   ! overwrite wave directions with usr set ones (0=no,1=yes)
            ```
            
            Need to change this to
            
            ```fortran
            waveicedatadir='../../waveice_data/'
            character(10)            :: fname_alp='alp_coeffs'
             !character(24)            :: fname_ww3='waves/ww3.197803_full.nc'
             character(14)            :: fname_ww3='waves/ww3.1978'
            integer, parameter       :: OVERWRITE_DIRS = 0   ! overwrite wave directions with usr set ones (0=no,1=yes)
            ```
            

# Dependencies

The code is initiated within `ice_step_mod`, in the subroutine `step_therm2` 

## `CICE_RUN`

## `ice_forcing`

## `ice_step_mod`

# Reference list:

Bennetts, L., O’Farrell, S. and Uotila, P. (2017) ‘Brief communication: impacts of ocean-wave-induced breakup of Antarctic sea ice via thermodynamics in a stand-alone version of the CICE sea-ice model’, The Cryosphere

Horvat, C. and Tziperman, E. (2015) ‘A prognostic model of the sea-ice floe size and thickness distribution’, The Cryosphere, 9(6), pp. 2119–2134. doi:10.5194/tc-9-2119-2015.

Meylan, M.H., Bennetts, L.G. and Kohout, A.L. (2014) ‘In situ measurements and analysis of ocean waves in the Antarctic marginal ice zone.’, Geophysical Research Letters, 41(14), pp. 5046–5051. doi:10.1002/2014GL060809.)

Roach, L.A. et al. (2018) ‘An Emergent Sea Ice Floe Size Distribution in a Global Coupled Ocean-Sea Ice Model’, Journal of Geophysical Research: Oceans, 123(6), pp. 4322–4337. doi:10.1029/2017JC013692.. doi:10.5194/tc-11-1035-2017.

Williams, T.D. et al. (2013) ‘Wave–ice interactions in the marginal ice zone. Part 1: Theoretical foundations’, Ocean Modelling, 71, pp. 81–91. doi:10.1016/j.ocemod.2013.05.010.

Williams, T.D. et al. (2013) ‘Wave–ice interactions in the marginal ice zone. Part 2: Numerical implementation and sensitivity studies along 1D transects of the ocean surface’, Ocean Modelling, 71, pp. 92–101. doi:10.1016/j.ocemod.2013.05.011.

# ❓FAQ

- Is this model in 2D?
    
    When the wave direction override is turned on it is a 1D model, as waves only propagate directly southward. However, when the wave direction data is used the model does become 2D.
    
    It uses 2D [Chebyshev interpolation](https://www.notion.so/Chebyshev-interpolation-1ebf5c9254504928a367e0d391175f5d), still have to check with [Meylan et al. 2014](https://www.notion.so/Meylan-et-al-2014-2c785c506c784e4c85fa00270fff63ff) and [Williams et al. 2013a,b](https://www.notion.so/Williams-et-al-2013a-b-5eae5120f7164a1ba164154ec1189e9f) what dimensions these schemes are in.
    
- What is the wave spectrum made up of?
    
    They construct a 2D Bretschneider spectrum, defined by a significant wave height and a peak period, which is propagated in the mean wave direction.
    
- Are waves propagated?
    
    Yes, from the closest ice-free cell. It uses a 1D wave scattering advection scheme, so it takes in wave direction but theses directions don't evolve over time.
    
- What output does this module produce?
    
    Change in floe size.
    
- What FSD do they use?
    
    They use the FSD routine in CICE6 (Roach et al., 2017). However it would be possible to implement a split power law to describe the FSD following Toyota et al. 2017 (as seen in Bennetts et al., 2017).
    
- Does it use wave direction?
    
    They do use wave direction from the closest ice free cell. The wave spectrum is 'advected' through the ice. I say 'advected' as it is really only attenuation rather than full advection. And the wave direction does not change as the waves collide with floes and other waves.
    
- Variable dimensions
    
    ```fortran
    real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
             ifd     , & ! ice floe diameter (m)
             swh     , & ! significant wave height (m)
             mwd     , & ! mean wave direction (Rads)
             ppd         ! wave peak period (s)
    ```
    
- What are the attenuation settings?
    
    There are two kinds of attenuation present in this module, the simplest is based off of observations from Meylan et al. (2014). The more complex attenuation model has not currently been tested in CICE6 and is based off of the work by Williams et al. (2013a,2013b).
    
