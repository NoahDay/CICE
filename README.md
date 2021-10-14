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
    
- Summary of `ice_floe`
    
    Found in `source_floe_test`
    
    Authors: Luke Bennetts, Siobhan O'Farrell
    
    - `init_floe`, initialise ice floe tracer (call prior to reading restart data)
        
        ```fortran
        if (trim(runtype) == 'continue') restart_floe = .true.
        
              if (restart_floe) then
                 call read_restart_floe
              else
                 call init_floe_0
        			endif
        ```
        
    - `init_floe_0`, initialise ice floe tracer
        
        
        ```fortran
        trcrn(:,:,nt_ifloe,:,:) = c0
        ifd(:,:,:) = c0
        	 do j = 1, ny_block
        	  do i = 1, nx_block
        	   trcrn(i,j,nt_ifloe,:,:) = max_floediam !c300
        	   ifd(i,j,:)              = max_floediam !c300
        	  enddo
        	 enddo
        ```
        
    - `increment_floe`, increase ice floe tracer by scaled timestep length, wave spectrum tolerance, spectrum passed through an ice-covered ocean, sine mean wave direction for row, mean wave directions, WIM termination flags, attenuation parameters, points in frequency and angular domains.
        - Variables:
            - m_prams_waveice
            - m_waveice
            - m_waveatn
            
            Parameters:
            
            - freqency min/max (`fmin = 1/16`, `fmax = 1/6` )
            - angular frequencies ($om_1 = 2\pi*fmin, om_2 = 2\pi \times fmax$)
            - $om_0 = \frac{om_2 - om_1}{nw_{in}-1}$
            - frequency, `om_in`, `T`
            - direction (not used yet), `th_in`
            - initial energy spec, `S_init_in`
            - wavelength/number, `lam_wtr_in`, `k_wtr_in`
            - dummy spectral moments (`dum_sm0`, `dum_sm2`)
            - wave spectrum tolerance, `ws_tol` = 0.1
            - spectrum passed through ice-covered ocean, `wspec_row`
            - sine mean wave direction for row, `mwd_row`
            - array to hold values for calculation of mean wave directions, `mwd_hld`
            - WIM termination flag, `tmt`, `tmt_hld`
            - attenuation, `idd_alp`, `dum_alp`, `fname_alp`
            - points in frequency and angular domains, `nw_in` = 31, `nth_in` = 1
        
        Method:
        
        1. Define wave mask edge and propagate waves into ice
        2. Define attenuation coefficient coefficients
        3. Calculate wave numbers and wavelengths
            
            Is direction used?
            
        
        ```fortran
        om_in(lp_i) = om1 + (lp_i-1)*om_0 ! angular frequency
        T(lp_i) = 2d0*pi/om_in(lp_i) ! frequency
        lam_wtr_in(lp_i) = gravity*(T(lp_i)**2d0)/2d0/pi ! wavelength
        k_wtr_in(lp_i) = 2d0*pi/lam_wtr_in(lp_i) !wavenumber
        th_in(lp_i) = -pi/2 + (lp_i-1)*pi/(nth_in-1) ! direction
        ```
        
        1. Initialise the wave spectrum using Bretschneider
        
        ```fortran
        ! if there are waves
        S_init_in(lp_i) = SDF_Bretschneider(om_in(lp_i),0,loc_swh(i,j),loc_ppd(i,j))
        ```
        
        1. A1. Use WIM to update wave spectrum and floe sizes
        
        ```fortran
        if (do_coupled.ne.0) then
        	call sub_Balance(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j), nw_in,nth_in,om_in,th_in,k_wtr_in,S_init_in,wspec_row_hld(i,:),tmt(i),nu_diag)
        		 else
        	call sub_Uncoupled(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j), nw_in,nth_in,om_in,th_in,k_wtr_in,S_init_in,wspec_row_hld(i,:),tmt(i),nu_diag)
        endif
        ```
        
        1. A2. redistribute energy according to mean directions
            
            a) For the LH boundary
            
        
        ```fortran
        ! Calculate a dummy swh to weight directions:
        dum_sm0 = fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
        
        dum_sm0 = 4d0*(dum_sm0**0.5d0)
        ```
        
        1. Find the direction the waves need to be advected in
        
         b) do the same for the RH boundary
        
        1. A3. Calculate the mean parameters
        
        ```fortran
        dum_sm0 = fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
        dum_sm2 = fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,2,nu_diag)
        loc_swh(i,j)   = 4d0*(dum_sm0**0.5d0)
        loc_ppd(i,j)   = 2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
        ```
        
        1. B1. Use WIM to update wave spectrum and floe size
        
        ```fortran
        if (do_coupled.ne.0) then
          call sub_Balance(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j),nw_in, nth_in,om_in,th_in,k_wtr_in,wspec_row(i,:),wspec_row_hld(i,:),tmt(i),nu_diag)
          else 
           call sub_Uncoupled(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j),nw_in,nth_in,om_in,th_in,k_wtr_in,wspec_row(i,:),wspec_row_hld(i,:),tmt(i),nu_diag)
        endif ! ENDIF (do_coupled.ne.0)
        ```
        
        continues same as the A routine.
        
    - Series of subroutines that write/read Fortran unformatted data files
        
        these subroutines w
        
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
        
        Also calculates the length of the MIZ.
        
        - More details
            
            If there is no ice, then set the maximum floe size to the floe size of pancake ice.
            
            If there is some ice, $mass  = reldens*hice$, 
            
            $fr    = Y*(hice**3)/12/water\_density/gravity/(1-poisson**2)$
            
            - Not sure what this quantity is or what formula it follows
                
                ```fortran
                dom = om(2)-om(1)
                 if (nth.eq.1) then
                  dth = 3d0         ! for Simpson's rule
                 else 
                  dth = th(2)-th(1)
                 end if
                ```
                
                ```fortran
                do lp_i=1,nw
                  kappa           = om(lp_i)**2d0/gravity
                  k_ice(lp_i)     = zero(0d0,max(kappa/tanh(kappa),sqrt(sqrt(kappa*mass/fr))), &
                                         toli,toli,fn_DispRel_ice_inf)
                  lam_ice(lp_i)   = 2d0*pi/k_ice(lp_i)
                 end do
                ```
                
                Then writes "Wave number spectra: k_wtr(1) → k_wrt(nw), lam_wtr = 2pi/k_wtr(1) → 2pi/k_wtr(nw), k_ice(1) → k_ice(nw), lam_ice(1)→lam_ice(nw)"
                
                ```bash
                L = 0  ! for test of incident spectrum
                
                call sub_StrainSpec(S_init, Es)
                call sub_WavelenSpec(S_init, lam_init)
                ```
                
                Then writes the incident spectrum output
                
                ```fortran
                write(idl,*) '>>>>>>> INCIDENT SPECTRUM:'
                  write(idl,*) 'S_init   = ', S_init(1), '->', S_init(nw*nth)
                  call sub_StrainSpec(S_init, Es_init)
                  write(idl,*) 'lam_init = ', lam_init
                  write(idl,*) 'Es_init  = ', Es_init
                ```
                
                If Inc strain < fracture threshold then, nothing to do but advect waves
                
                ```fortran
                if (ATTEN_METH.eq.1) then ! for every lp_i
                   alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                  else
                   alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                  end if
                !  S_attn(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
                  do lp_j=1,nth
                   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
                     exp(-alpha(lp_i)*cos(th(lp_j))*L)
                  end do
                ```
                
                else, there is potential for fracture
                
                ```fortran
                L = Lcell ! initialise propagation length
                 
                 call sub_WavelenSpec(S_init, lam_init)
                 
                 floe_sz_max = zero(i0,i1,toli,toli,fn_wlng)
                
                if (ATTEN_METH.eq.1) then ! for every lp_i
                   alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                  else
                   alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                  end if
                !  S_attn(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
                  do lp_j=1,nth
                   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
                     exp(-alpha(lp_i)*cos(th(lp_j))*L)
                  end do
                
                call sub_Strain(S_attn,Es)
                ```
                
                Waves strong enough to break ice at end of cell
                
                ```fortran
                do lp_i=1,nw*nth
                   S_attn_out(lp_i) = S_attn(lp_i)
                  end do
                ```
                
                Else, search for a point at which waves can break ice
                
                ```fortran
                Lc = zero(i0,Lcell,toli,toli,fn_proplng)
                
                if (idl.ne.0.and.cmt.ne.0) then
                   write(idl,*) '>>>>>>> RESULT:'
                   write(idl,*) '                      --> Breakage occurs within cell ', &
                   Lc/1000d0,'km' !, '(zero chk: ', fn_proplng(Lc),')'
                !   write(idl,*) '<<<---------------------------------------------<<<'
                  end if
                  L = Lc
                  floe_sz_max = zero(i0,i1,toli,toli,fn_wlng)
                
                   if (ATTEN_METH.eq.1) then ! at every lp_i
                    alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                   else
                    alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                   end if
                !   S_attn_out(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
                   do lp_j=1,nth
                    S_attn_out(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
                     exp(-alpha(lp_i)*cos(th(lp_j))*L)
                   end do      
                
                tmt = 1
                ```
                
    - `fn_apxAttn`, approximate the attenuation coefficient
        
        ```fortran
        Rsq = 0.002       ! Reflected energy (10s period, 1m thickness)
         A = -2*log(1-Rsq) ! long-floe limit
        
         B = 5/(lambda/2)**gamma ! tanh(5) = 1
        
         fn_apxAttn = A*tanh(B*(fl_len**gamma))
        ```
        
    - `fn_splitPDF`, the split PDF of [Toyota et al. 2017](https://www.notion.so/Toyota-et-al-2011-fd26dea41caa4564b03a4dcf47c37af0)
        
        ```fortran
        if (D.ge.floe_sz_crit) then
            fn_splitPDF = (1-P0)*beta1*gamma1/D**(gamma1+1d0)
         elseif (D.lt.floe_sz_crit.and.D.ge.floe_sz_min) then    
            fn_splitPDF = P0*beta0*gamma0/D**(gamma0+1d0)
         else
            fn_splitPDF=0
         end if
        ```
        
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
        
        ```fortran
        !mid-points
         dD = (floe_sz_max-floe_sz_min)/(res+1)
         do loop_w=1,res+1
          dum_D_vec(loop_w) = floe_sz_min + (loop_w-1)*dD
         end do
         do loop_w=1,res
          D_vec(loop_w) = (dum_D_vec(loop_w+1)+dum_D_vec(loop_w))/2d0
         end do
         
         ! FSD params
         
         beta0 = 1/(floe_sz_min**(-gamma0) - floe_sz_max**(-gamma1))
         beta1 = floe_sz_crit**gamma1
         
         if (floe_sz_max.ge.floe_sz_crit) then
          P0 = 1d0-0.05*((floe_sz_max/floe_sz_crit)**gamma1)
         else
          P0=0d0
         end  if
          
         if (P0.ge.1.0) then
          P0=1d0
         elseif (P0.le.0.0) then
          P0=0d0
         end if
         
         do loop_w=1,res
          if (ATTEN_MODEL.eq.0) then
           alpha_vec(loop_w) = fn_apxAttn(D_vec(loop_w),lam,1d0)
          elseif (ATTEN_MODEL.eq.1) then
           alpha_vec(loop_w) = fn_Attn_WIM_v1(dum_om, hice, dum_idl)
          end if 
          p_vec(loop_w) = fn_splitPDF(D_vec(loop_w),beta0,beta1,P0)
         end do 
         
         fn_IntAttn = dD*conc*dot_product(alpha_vec,p_vec/D_vec)
        ```
        
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
    

# 🏗️ Structure

- `cice_with_waves`
    - `bld`
        
        Make files.
        
    - `csm_share`
        
        This module was borrowed from the CCSM shared code library and modified to work within CICE without other CCSM shared code.  All such changes are marked with `!echmod`  -- E.C. Hunke, 2007
        
    - `doc`
    - `drivers`
    - `input_templates`
    - `mpi`
        
        No cases of `wave`.
        
    - `serial`
        
        No cases of `wave` found.
        
    - `source_floe_test`
    - `wave_ice_code`
        - `m_fzero.F90`
        - `m_prams_waveice.F90`
            
            Author: Luke Bennetts
            
            Line 82:
            
            ```fortran
            ! WAVES
             
             integer, parameter                    :: WAVE_METH = 0   ! inc waves (0=user,1=ww3) 
             real(kind=8), allocatable             :: ww3_lat(:,:), ww3_lon(:,:), ww3_tm(:,:)
             real(kind=8), allocatable             :: ww3_swh(:,:), ww3_fp(:,:), ww3_dir(:,:), &
             				ww3_swh_full(:,:,:), ww3_fp_full(:,:,:), ww3_dir_full(:,:,:)
             integer, parameter                    :: nww3_dt = 2 ! ww3 time step relative to cice				
              
             ! DATA
             
             character(19)            :: waveicedatadir='../../waveice_data/'
             character(10)            :: fname_alp='alp_coeffs'
             !character(24)            :: fname_ww3='waves/ww3.197803_full.nc'
             character(14)            :: fname_ww3='waves/ww3.1978'
             integer, parameter       :: OVERWRITE_DIRS = 0   ! overwrite wave directions with usr set ones (0=no,1=yes)
            ```
            
        - `m_waveattn.F90`
            
            Author: Luke Bennetts
            
            ```fortran
            use m_prams_waveice
            ```
            
            - `fn_Attn_MBK`, lines 59-85
                
                The attenuation model derived by Meylan et al. (2014), Geophys Res Lett
                
                $$\alpha = \beta_0 + \beta_1 \times om^2$$
                
                $$fn\_Attn\_MBK = attn\_fac \times \beta_0(dum\_om^2) + \beta_1(dum\_om^4)$$
                
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
                
        - `m_waveice.F90`
            
            Author: Luke Bennetts
            
            ```fortran
            use m_fzero
            use m_prams_waveice
            use m_waveattn
            ```
            
            - Routine: `sub_Uncoupled`
                
                Attenuation does not depend on floe sizes.
                
                Line 143: `! No ice = nothing to do`
                
                ```fortran
                if (conc.lt.tolice.or.hice.lt.tolh) then ! conc = dum_conc (cell concentration)
                  floe_sz_max = floe_sz_pancake
                  do lp_i=1,nw
                   S_attn_out(lp_i) = dum_S_init(lp_i) ! S_attn is attn spectrum, l_pi/l_pj are numerical params
                  end do
                  
                  if (idl.ne.0.and.cmt.ne.0) then ! idl
                   write(idl,*) '>>>>>>> RESULT:'
                   if (hice.eq.0d0) then
                    write(idl,*) '                       --> no ice in this cell: h<', tolh  
                   else
                    write(idl,*) '                       --> no ice in this cell: c<', tolice
                   end if
                  write(idl,*) '                           tmt =', tmt   
                  write(idl,*) '                   floe_sz_max =', floe_sz_max 
                !   write(idl,*) '<<<---------------------------------------------<<<'
                  end if
                ```
                
                `else ! there is some ice`
                
                ```fortran
                mass  = reldens*hice 
                 fr    = Y*(hice**3)/12/water_density/gravity/(1-poisson**2)
                ```
                
                - Write
                    
                    ```bash
                     if (idl.ne.0.and.cmt.ne.0) then !
                      write(idl,*) '>>>>>>>>> PARAMS:'
                      write(idl,*) 'mass, fr = ', mass, fr
                      write(idl,*) 'E_crit=epsc*sqrt(-2/log(Pc))      = ', epsc*sqrt(-2/log(Pc))
                      write(idl,*) 'nw, nth', nw, nth
                      write(idl,*) 'om=',dum_om(1),'->',dum_om(nw)
                      write(idl,*) 'th=',dum_th(1),'->',dum_th(nth)
                     endif
                    ```
                    
                
                Line 278, `Ice strain < fracture threshold = nothing to do but advect waves`
                
                ```fortran
                if (Es.lt.epsc*sqrt(-2/log(Pc))) then 
                
                 L = Lcell ! initialise propagation length
                 
                 floe_sz_max = floe_sz_init ! initialise the maximum floe size
                  
                 do lp_i=1,nw
                  alpha(lp_i)  = conc*fn_Attn_MBK(om(lp_i))/0.7d0 ! attn coeff
                  do lp_j=1,nth
                   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
                     exp(-alpha(lp_i)*cos(th(lp_j))*L) ! attn spectrum calculation
                  end do   
                 end do
                 
                 do lp_i=1,nw*nth
                  S_attn_out(lp_i) = S_attn(lp_i)
                 end do
                ```
                
                - Write
                    
                    ```bash
                     if (idl.ne.0.and.cmt.ne.0) then
                      write(idl,*) '>>>>>>> RESULT:'
                      write(idl,*) '                      --> S_init isnt strong enough to break ice: ', &
                       Es, epsc*sqrt(-2/log(Pc))
                      write(idl,*) '                          tmt =', tmt
                      write(idl,*) '                  floe_sz_max =', floe_sz_max  
                    !   write(idl,*) '<<<---------------------------------------------<<<'
                     end if
                     
                     !tmt = 0      ! Es proportional to hice -> so can't just switch off (LB Apr 14)
                    ```
                    
                
                ```fortran
                else ! potential for fracture
                L = Lcell ! initialise propagation length
                  
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) '>>>>>>> AT END OF CELL:'
                  !write(idl,*) 'D1 (init) = ', floe_sz_init, ': empirical breaking size = ', floe_sz_brk_emp
                 end if
                    
                 do lp_i=1,nw
                  alpha(lp_i)  = conc*fn_Attn_MBK(om(lp_i))/0.75d0
                  do lp_j=1,nth
                   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
                     exp(-alpha(lp_i)*cos(th(lp_j))*L)
                  end do     
                 end do
                 
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) 'alpha = ', alpha(1), '->', alpha(nw)
                  write(idl,*) 'S_attn = ', S_attn(1), '->', S_attn(nw)
                 end if
                 
                 call sub_StrainSpec(S_attn, Es)
                  
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) 'Es = ', Es
                 end if
                 
                 ! Set attenuated spectrum (as not floe size dependent)
                 do lp_i=1,nw*nth
                  S_attn_out(lp_i) = S_attn(lp_i)
                 end do
                ```
                
                ```fortran
                if (floe_sz_brk_emp.gt.floe_sz_init) then  !!! floes already smaller than breakage size
                
                floe_sz_max = floe_sz_init
                  !tmt = 0   
                  if (idl.ne.0.and.cmt.ne.0) then
                   write(idl,*) '>>>>>>> RESULT:'
                   write(idl,*) '                       --> Breakage occurs: location inconsequential'
                   write(idl,*) '                           initial floe size', floe_sz_init, '<', &
                      											floe_sz_brk_emp
                   write(idl,*) '                           tmt =', tmt 
                   write(idl,*) '                   floe_sz_max =', floe_sz_max     											
                !   write(idl,*) '<<<---------------------------------------------<<<'
                  end if
                ```
                
                ```fortran
                elseif (Es.ge.epsc*sqrt(-2/log(Pc))) then   !!! waves strong enough to break ice at end of cell
                
                  floe_sz_max = floe_sz_brk_emp
                  !tmt = 0    
                  if (idl.ne.0.and.cmt.ne.0) then
                   write(idl,*) '>>>>>>> RESULT:'
                   write(idl,*) '                       --> Breakage occurs at end of cell',L/1000d0,'km'
                   write(idl,*) '                           tmt =', tmt    
                   write(idl,*) '                   floe_sz_max =', floe_sz_max
                !   write(idl,*) '<<<---------------------------------------------<<<'
                  end if
                ```
                
                ```fortran
                else        !!! search for point at which waves can break ice
                  Lc = zero(i0,Lcell,toli,toli,fn_proplng_uncoupled)
                  L           = Lc
                  floe_sz_max = (Lc*floe_sz_brk_emp + (Lcell-Lc)*floe_sz_init)/Lcell
                  tmt         = 1   ! waves attenuated (in all liklihood)
                  if (idl.ne.0.and.cmt.ne.0) then
                   write(idl,*) '>>>>>>> RESULT:'
                   write(idl,*) '                      --> Breakage occurs within cell ', &
                   Lc/1000d0,'km' !, '(zero chk: ', fn_proplng_uncoupled(Lc),')'
                   write(idl,*) '                          tmt =', tmt   
                   write(idl,*) '                  floe_sz_max =', floe_sz_max
                !   write(idl,*) '<<<---------------------------------------------<<<'
                  end if
                  
                 end if ! (Es.ge.epsc*sqrt(-2/log(Pc)))
                
                 end if ! (Esinit.lt.Ec)
                ```
                
                ```fortran
                Deallocate things
                
                end if  ! (conc.lt.tolice.or.hice.eq.0d0)
                
                 end subroutine sub_Uncoupled
                ```
                
            - Routine: `sub_Balance`, lines 427-771.
                
                Line 507:
                
                ```fortran
                ! No ice = nothing to do
                if (conc.lt.tolice.or.hice.lt.tolh) then
                  floe_sz_max = 5d0
                  do lp_i=1,nw
                   S_attn_out(lp_i) = dum_S_init(lp_i)
                  end do
                  
                  if (idl.ne.0.and.cmt.ne.0) then
                   write(idl,*) '>>>>>>> RESULT:'
                   if (hice.eq.0d0) then
                    write(idl,*) '                       --> no ice in this cell: h<', tolh  
                   else
                    write(idl,*) '                       --> no ice in this cell: c<', tolice
                   end if
                !   write(idl,*) '<<<---------------------------------------------<<<'
                  end if
                ```
                
                ```fortran
                else ! there is some ice
                
                mass  = reldens*hice 
                 fr    = Y*(hice**3)/12/water_density/gravity/(1-poisson**2)
                   
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) '>>>>>>>>> PARAMS:'
                  write(idl,*) 'mass, fr = ', mass, fr
                  write(idl,*) 'E_crit=epsc*sqrt(-2/log(Pc))      = ', epsc*sqrt(-2/log(Pc))
                  write(idl,*) 'nw, nth', nw, nth
                  write(idl,*) 'om=',dum_om(1),'->',dum_om(nw)
                  write(idl,*) 'th=',dum_th(1),'->',dum_th(nth)
                 endif 
                  
                 allocate(om(1:nw))
                 allocate(th(1:nth))
                 
                 allocate(k_wtr(1:nw))
                 
                 allocate(k_ice(1:nw))
                 allocate(lam_ice(1:nw))
                 
                 allocate(S_init(1:nth*nw))
                 allocate(S_attn(1:nw*nth))
                 
                 allocate(alpha(1:nw))
                 
                 do lp_i=1,nw
                  om(lp_i)     = dum_om(lp_i)
                  k_wtr(lp_i)  = dum_k_wtr(lp_i)
                 end do
                 
                 do lp_i=1,nth*nw 
                  S_init(lp_i) = dum_S_init(lp_i)
                 end do
                 
                 do lp_i=1,nth
                  th(lp_i)     = dum_th(lp_i)
                 end do 
                 
                 dom = om(2)-om(1)
                 if (nth.eq.1) then
                  dth = 3d0         ! for Simpson's rule
                 else 
                  dth = th(2)-th(1)
                 end if 
                 
                 do lp_i=1,nw
                  kappa           = om(lp_i)**2d0/gravity
                  k_ice(lp_i)     = zero(0d0,max(kappa/tanh(kappa),sqrt(sqrt(kappa*mass/fr))), &
                                         toli,toli,fn_DispRel_ice_inf)
                  lam_ice(lp_i)   = 2d0*pi/k_ice(lp_i)
                 end do
                 
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) '>>>>>>> WAVE NUMBER SPECTRA:'
                  write(idl,*) 'k_wtr   = ', k_wtr(1), '->', k_wtr(nw)
                  write(idl,*) 'lam_wtr = ', 2*pi/k_wtr(1), '->', 2*pi/k_wtr(nw)
                  write(idl,*) 'k_ice   = ', k_ice(1), '->', k_ice(nw)
                  write(idl,*) 'lam_ice = ', lam_ice(1), '->', lam_ice(nw)
                 endif 
                 
                 L = 0  ! for test of incident spectrum
                 
                 call sub_StrainSpec(S_init, Es)
                 
                 call sub_WavelenSpec(S_init, lam_init)
                 
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) '>>>>>>> INCIDENT SPECTRUM:'
                  write(idl,*) 'S_init   = ', S_init(1), '->', S_init(nw*nth)
                  call sub_StrainSpec(S_init, Es_init)
                  write(idl,*) 'lam_init = ', lam_init
                  write(idl,*) 'Es_init  = ', Es_init
                 end if
                ```
                
                ```fortran
                ! Inc strain<fracture threshold =  nothing to do but advect waves
                
                if (Es.lt.epsc*sqrt(-2/log(Pc))) then
                
                L = Lcell ! initialise propagation length
                 
                 floe_sz_max = floe_sz_init
                 
                 do lp_i=1,nw
                  if (ATTEN_METH.eq.1) then
                   alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                  else
                   alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                  end if
                !  S_attn(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
                  do lp_j=1,nth
                   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
                     exp(-alpha(lp_i)*cos(th(lp_j))*L)
                  end do   
                 end do
                 
                 do lp_i=1,nw*nth
                  S_attn_out(lp_i) = S_attn(lp_i)
                 end do
                 
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) '>>>>>>> RESULT:'
                  write(idl,*) '                      --> S_init isnt strong enough to break ice: ', &
                   Es, epsc*sqrt(-2/log(Pc))
                !   write(idl,*) '<<<---------------------------------------------<<<'
                 end if
                 
                 tmt = 0      ! Es proportional to hice -> so can't just switch off (LB Apr 14)
                ```
                
                ```fortran
                else ! potential for fracture
                L = Lcell ! initialise propagation length
                 
                 call sub_WavelenSpec(S_init, lam_init)
                 
                 floe_sz_max = zero(i0,i1,toli,toli,fn_wlng)
                  
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) '>>>>>>> AT END OF CELL:'
                  write(idl,*) 'D1 = ', floe_sz_max, ': f(z) = ', fn_wlng(floe_sz_max)
                 end if
                   
                 do lp_i=1,nw
                  if (ATTEN_METH.eq.1) then
                   alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                  else
                   alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                  end if
                !  S_attn(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
                  do lp_j=1,nth
                   S_attn(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
                     exp(-alpha(lp_i)*cos(th(lp_j))*L)
                  end do     
                 end do
                 
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) 'alpha = ', alpha(1), '->', alpha(nw)
                  write(idl,*) 'S_attn = ', S_attn(1), '->', S_attn(nw)
                 end if
                 
                 call sub_StrainSpec(S_attn, Es)
                  
                 if (idl.ne.0.and.cmt.ne.0) then
                  write(idl,*) 'Es = ', Es
                 end if
                ```
                
                ```fortran
                if (Es.ge.epsc*sqrt(-2/log(Pc))) then   !!! waves strong 
                 !!!                   enough to break ice at end of cell
                if (idl.ne.0.and.cmt.ne.0) then
                   write(idl,*) '>>>>>>> RESULT:'
                   write(idl,*) '                       --> Breakage occurs at end of cell',L/1000d0,'km'
                !   write(idl,*) '<<<---------------------------------------------<<<'
                  end if
                  do lp_i=1,nw*nth
                   S_attn_out(lp_i) = S_attn(lp_i)
                  end do
                ```
                
                ```fortran
                else        !!! search for point at which waves can break ice
                  Lc = zero(i0,Lcell,toli,toli,fn_proplng)
                  if (idl.ne.0.and.cmt.ne.0) then
                   write(idl,*) '>>>>>>> RESULT:'
                   write(idl,*) '                      --> Breakage occurs within cell ', &
                   Lc/1000d0,'km' !, '(zero chk: ', fn_proplng(Lc),')'
                !   write(idl,*) '<<<---------------------------------------------<<<'
                  end if
                  L = Lc
                  floe_sz_max = zero(i0,i1,toli,toli,fn_wlng)
                  do lp_i=1,nw
                   if (ATTEN_METH.eq.1) then
                    alpha(lp_i)  = fn_IntAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                   else
                    alpha(lp_i)  = fn_AvAttn(floe_sz_max, lam_ice(lp_i), om(lp_i))
                   end if
                !   S_attn_out(lp_i) = S_init(lp_i)*exp(-alpha(lp_i)*L)
                   do lp_j=1,nth
                    S_attn_out(lp_i+nw*(lp_j-1)) = S_init(lp_i+nw*(lp_j-1))* &
                     exp(-alpha(lp_i)*cos(th(lp_j))*L)
                   end do      
                  end do
                  tmt = 1
                ```
                
            - Routine: `fn_apxAttn`, lines 777-820.
                
                Approximate the attenuation coefficient.
                
                ```fortran
                Rsq = 0.002       ! Reflected energy (10s period, 1m thickness)
                 A = -2*log(1-Rsq) ! long-floe limit
                
                 B = 5/(lambda/2)**gamma ! tanh(5) = 1
                
                 fn_apxAttn = A*tanh(B*(fl_len**gamma))
                ```
                
            - Routine: `fn_splitPDF`, lines 826-868.
                
                The split PDF of Toyota et al. (2011).
                
                ```fortran
                if (D.ge.floe_sz_crit) then
                    fn_splitPDF = (1-P0)*beta1*gamma1/D**(gamma1+1d0)
                 elseif (D.lt.floe_sz_crit.and.D.ge.floe_sz_min) then    
                    fn_splitPDF = P0*beta0*gamma0/D**(gamma0+1d0)
                 else
                    fn_splitPDF=0
                 end if
                ```
                
            - Routine: `fn_AvAttn`, lines 874-926.
                
                The averaged attenuation coefficient with respect to floe diameter.
                
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
                
            - Routine: `fn_IntAttn`, lines 932-1008.
                
                The integrated attenuation coefficient with respect to floe diameter.
                
                ```fortran
                ! mid-points
                 dD = (floe_sz_max-floe_sz_min)/(res+1)
                 do loop_w=1,res+1
                  dum_D_vec(loop_w) = floe_sz_min + (loop_w-1)*dD
                 end do
                 do loop_w=1,res
                  D_vec(loop_w) = (dum_D_vec(loop_w+1)+dum_D_vec(loop_w))/2d0
                 end do
                 
                 ! FSD params
                 
                 beta0 = 1/(floe_sz_min**(-gamma0) - floe_sz_max**(-gamma1))
                 beta1 = floe_sz_crit**gamma1
                 
                 if (floe_sz_max.ge.floe_sz_crit) then
                  P0 = 1d0-0.05*((floe_sz_max/floe_sz_crit)**gamma1)
                 else
                  P0=0d0
                 end  if
                  
                 if (P0.ge.1.0) then
                  P0=1d0
                 elseif (P0.le.0.0) then
                  P0=0d0
                 end if
                 
                 do loop_w=1,res
                  if (ATTEN_MODEL.eq.0) then
                   alpha_vec(loop_w) = fn_apxAttn(D_vec(loop_w),lam,1d0)
                  elseif (ATTEN_MODEL.eq.1) then
                   alpha_vec(loop_w) = fn_Attn_WIM_v1(dum_om, hice, dum_idl)
                  end if 
                  p_vec(loop_w) = fn_splitPDF(D_vec(loop_w),beta0,beta1,P0)
                 end do 
                 
                 fn_IntAttn = dD*conc*dot_product(alpha_vec,p_vec/D_vec)
                ```
                
            - Routine: `fn_proplng`, lines 1014-1066.
                
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
                
            - Routine: `fn_proplng_uncoupled`, lines 1072-1118.
                
                `fn_proplng` for problem when attn/floe length are uncoupled.
                
                ```fortran
                L = dum_L
                 
                 do loop_w=1,nw
                  alpha(loop_w)  = conc*fn_Attn_MBK(om(loop_w))/0.7d0 
                  S_attn(loop_w) = S_init(loop_w)*exp(-alpha(loop_w)*L)
                 end do 
                 
                 call sub_StrainSpec(S_attn, dum_Es)
                 
                 fn_proplng_uncoupled = dum_Es - epsc*sqrt(-2/log(Pc))
                ```
                
            - Routine: `fn_SpecMoment`, lines 1124-1261.
                
                ```fortran
                dom_local = om_in(2)-om_in(1)
                 
                 if (nth_in.eq.1) then
                  dth_local = 3d0         ! for Simpson's rule
                 else 
                  dth_local = th_in(2)-th_in(1)
                 end if 
                 
                 dum_simp(1) = 1d0
                 dum_simp(nw_in) = 1d0
                 
                 do loop_w=2,nw_in-1,2
                  dum_simp(loop_w) = 4d0
                 end do
                
                 do loop_w=3,nw_in-1,2
                  dum_simp(loop_w) = 2d0
                 end do
                 
                 dum_simp_th(1) = 1d0
                 dum_simp_th(nth_in) = 1d0
                 
                 do loop_th=2,nth_in-1,2
                  dum_simp_th(loop_th) = 4d0
                 end do
                
                 do loop_th=3,nth_in-1,2
                  dum_simp_th(loop_th) = 2d0
                 end do
                 
                 do loop_w=1,nw_in
                  do loop_th=1,nth_in
                   wt_simp(loop_w+nw_in*(loop_th-1)) = dum_simp(loop_w)*dum_simp_th(loop_th)
                  end do
                 end do
                 
                 do loop_w=1,nw_in
                  do loop_th=1,nth_in
                   wt_int(loop_w+nw_in*(loop_th-1)) = (dom_local/3d0)*(dth_local/3d0)* &
                                              wt_simp(loop_w+nw_in*(loop_th-1))
                   dum_v(loop_w+nw_in*(loop_th-1))  = dum_S(loop_w+nw_in*(loop_th-1))* &
                                              (om_in(loop_w)**mom)
                  end do
                 end do
                 
                
                 fn_SpecMoment   = dot_product(wt_int,dum_v)
                 
                ```
                
            - Routine: `sub_WavelenSpec`, lines 1267-1381.
                
                `kappa` := frequency parameter
                
                `k_ice`
                
                From Williams et al. (2013a)
                
                The mean square displacement of the ice is approximately $\langle\eta_{ice}^2 \rang = m_0[\eta_{ice}]$, where
                
                 $m_0 = \int^\infty_0 S(\omega)W^2(\omega)d\omega$
                
                 $m_n = \int^\infty_0 S(\omega)\omega^nW^2(\omega)d\omega$
                
                so, 
                
                ```fortran
                mom0     = dot_product(wt_int,dum_u) 
                mom2     = dot_product(wt_int,dum_v)
                ```
                
                is just a discretised version of this integration.
                
                From Williams et al. (2013a)
                
                The wave number is
                
                $\kappa (\omega) = \frac{\omega ^2}{g}$
                
                ```fortran
                dum_simp(1) = 1d0
                 dum_simp(nw) = 1d0
                 
                 do loop_w=2,nw-1,2
                  dum_simp(loop_w) = 4d0
                 end do
                
                 do loop_w=3,nw-1,2
                  dum_simp(loop_w) = 2d0
                 end do
                 
                 dum_simp_th(1) = 1d0
                 dum_simp_th(nth) = 1d0
                 
                 do loop_th=2,nth-1,2
                  dum_simp_th(loop_th) = 4d0
                 end do
                
                 do loop_th=3,nth-1,2
                  dum_simp_th(loop_th) = 2d0
                 end do
                 
                 do loop_w=1,nw
                  do loop_th=1,nth
                   wt_simp(loop_w+nw*(loop_th-1)) = dum_simp(loop_w)*dum_simp_th(loop_th)
                  end do
                 end do
                 
                 do loop_w=1,nw
                  do loop_th=1,nth
                   wt_int(loop_w+nw*(loop_th-1)) = (dom/3d0)*(dth/3d0)*wt_simp(loop_w+nw*(loop_th-1))
                   !om_lng(loop_w+nw*(loop_th-1)) = om(loop_w)
                   if (hice.eq.0d0) then
                    !F(loop_w+nw*(loop_th-1)) = 1d0
                    dum_u(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))
                    dum_v(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))*(om(loop_w)**2)
                   else
                    if (Wsq_METH.eq.0) then
                     !F(loop_w+nw*(loop_th-1)) = 1d0
                     dum_u(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))
                     dum_v(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))*(om(loop_w)**2)
                    elseif (Wsq_METH.eq.1) then 
                     !F(loop_w+nw*(loop_th-1)) = k_ice(loop_w)/k_wtr(loop_w)
                     dum_u(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))* &
                      ((k_ice(loop_w)/k_wtr(loop_w))**2)
                     dum_v(loop_w+nw*(loop_th-1)) = dum_S(loop_w+nw*(loop_th-1))*(om(loop_w)**2)* &
                      ((k_ice(loop_w)/k_wtr(loop_w))**2)
                    end if
                   end if 
                  end do
                 end do
                
                 mom0     = dot_product(wt_int,dum_u) 
                 mom2     = dot_product(wt_int,dum_v)
                
                 om_crit = sqrt(mom2/mom0)
                 
                 !if (out_ct.eq.1) then
                 ! print*, 'wt_int=', wt_int
                 ! print*, 'S     =', dum_S
                 ! print*, 'om    =', om
                 ! print*, 'mom0,mom2,om_crit=', mom0,mom2,om_crit
                 ! out_ct = out_ct+1
                 !end if 
                 
                 kappa = om_crit**2d0/gravity
                  
                 wlng_crest = 2d0*pi/zero(0d0,max(kappa/tanh(kappa),sqrt(sqrt(kappa*mass/fr))), &
                                          toli,toli,fn_DispRel_ice_inf)
                ```
                
            - Routine: `sub_StrainSpec`, lines 1387-1490.
                
                Calculates the strain imposed on the ice by the wave spectrum.
                
                ```fortran
                dum_simp(1) = 1d0
                 dum_simp(nw) = 1d0
                 
                 do loop_w=2,nw-1,2
                  dum_simp(loop_w) = 4d0
                 end do
                
                 do loop_w=3,nw-1,2
                  dum_simp(loop_w) = 2d0
                 end do
                 
                 dum_simp_th(1) = 1d0
                 dum_simp_th(nth) = 1d0
                 
                 do loop_th=2,nth-1,2
                  dum_simp_th(loop_th) = 4d0
                 end do
                
                 do loop_th=3,nth-1,2
                  dum_simp_th(loop_th) = 2d0
                 end do
                 
                 do loop_w=1,nw
                  do loop_th=1,nth
                   wt_simp(loop_w+nw*(loop_th-1)) = dum_simp(loop_w)*dum_simp_th(loop_th)
                  end do
                 end do
                 
                 do loop_w=1,nw
                  do loop_th=1,nth
                   wt_int(loop_w+nw*(loop_th-1))   = (dom/3d0)*(dth/3d0)*wt_simp(loop_w+nw*(loop_th-1))
                   if (hice.eq.0d0) then
                    !F(loop_w+nw*(loop_th-1)) = 1d0
                    dum_vec(loop_w+nw*(loop_th-1)) = 0.25d0*(hice**2d0)*dum_S(loop_w+nw*(loop_th-1)) &
                         *(k_ice(loop_w)**4d0)
                   else
                    if (Wsq_METH.eq.0) then
                     !F(loop_w+nw*(loop_th-1)) = 1d0
                     dum_vec(loop_w+nw*(loop_th-1)) = 0.25d0*(hice**2d0)*dum_S(loop_w+nw*(loop_th-1)) &
                         *(k_ice(loop_w)**4d0)
                    elseif (Wsq_METH.eq.1) then 
                     !F(loop_w+nw*(loop_th-1)) = k_ice(loop_w)/k_wtr(loop_w)
                     dum_vec(loop_w+nw*(loop_th-1)) = 0.25d0*(hice**2d0)*dum_S(loop_w+nw*(loop_th-1)) &
                         *(k_ice(loop_w)**4d0)*((k_ice(loop_w)/k_wtr(loop_w))**2d0)
                    end if
                   end if 
                  end do
                 end do
                
                 mom0_eps = dot_product(wt_int,dum_vec)
                
                 Es = 2d0*sqrt(mom0_eps)
                ```
                
            - Routine: `fn_DispRel_ice_inf`, lines 1496-1532.
                
                Infinite depth ice-coupled dispersion relation.
                
                ```fortran
                fn_DispRel_ice_inf = (1d0-(mass*kappa)+(fr*(dum_k**4d0)))*dum_k - kappa
                ```
                
            - Routine: `fn_wlng`, lines 1538-1592.
                
                Calculate wavelength.
                
                $\text{fn\_wlng} = \frac{\lambda}{2} - D_1$
                
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
                
            - Routine: `SDF_Bretschneider`, lines 1598-1648.
                
                Bretschneider wave spectrum.
                
                ```fortran
                ! Hs - significant wave height
                ! Tm - peak period
                om_m = 2d0*pi/Tm
                 tau  = 2d0*pi/omega
                 
                 f1 = (5d0/16d0)*(Hs**2)*(om_m**4)
                 f2 = omega**(moment_no-5)
                 f3 = exp(-1.25d0*((tau/Tm)**4))
                 
                 SDF_Bretschneider = f1*f2*f3
                ```
                
- `forcing`
- `held_files`
- `rundir`
- `scripts`
