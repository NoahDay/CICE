!=======================================================================
!
!BOP
!
! !MODULE: ice_floe - Floe size distribution - from the cicewithwaves
!                                                  project
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors Noah Day, using the work from
!         Luke Bennetts and Siobhan O'Farrell
! date: 25-08-21
!
! !INTERFACE:
!
      module ice_floe
!
!
! !USES:
!
        use icepack_kinds
        use icepack_parameters, only: p01, p5, c0, c1, c2, c3, c4, c10
        use icepack_parameters, only: bignum, puny, gravit, pi
        use icepack_tracers, only: nt_fsd
        use icepack_warnings, only: warnstr, icepack_warnings_add,  icepack_warnings_aborted
        use ice_domain_size, only: ncat, nfsd, nfreq
!use ice_state ! noah day wim 004
        !use ice_constants
        use ice_fileunits
        use ice_read_write
        use ice_restart_shared, only: lenstr, restart_dir, restart_file, &
                               pointer_file, runtype
        use ice_communicate, only: my_task, master_task
        use ice_exit, only: abort_ice
        use m_prams_waveice, only: waveicedatadir, fname_ww3, WAVE_METH, ww3_lat, ww3_lon, &
       ww3_dir, ww3_tm, ww3_swh, ww3_fp, ATTEN_METH, ATTEN_MODEL, attn_fac, do_coupled, &
       OVERWRITE_DIRS, ww3_dir_full, ww3_swh_full, ww3_fp_full, nww3_dt
       !use netcdf
       use ice_constants, only: eps1, eps3
!
!EOP
!
  implicit none

  !logical (kind=log_kind) :: &
     !tr_fsd,       & ! if .true., use floe tracer
     !restart_fsd,   & ! if .true., read floe tracer restart file
     !tr_wspec,       & ! if .true., use wave spec tracer
     !restart_spec      ! if .true., read wave spec tracer restart file

!=======================================================================

    contains

!=======================================================================
!BOP
!
! !ROUTINE: init_floe
!
! !DESCRIPTION:
!
!  Initialize ice floe tracer (call prior to reading restart data)
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
subroutine init_floe
!
! !USES:
 use ice_arrays_column, only: floe_rad_c
 use ice_state, only: trcrn
 use icepack_tracers, only: nt_fsd, tr_fsd
!
!EOP
!
!  integer (kind=int_kind)            :: i, j


!real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
!   trcrn           ! tracer array
! noah day wim 005 -------

    call init_floe_0

! --------------------

!  if (trim(runtype) == 'continue') restart_fsd = .true.



!     call init_floe_0
!         trcrn(i,j,nt_fsd,:,:) = c0
!         do j = 1, wavemask
!          do i = 1, nx_block
!           trcrn(i,j,nt_fsd,:,:) = c300
!          enddo
!         enddo
!         floe_rad_c(i,j,:) = c0
!         do j = 1, wavemask
!          do i = 1, nx_block
!           floe_rad_c(i,j,:) = c300
!          enddo
!         enddo


  end subroutine init_floe

!=======================================================================
!BOP
!
! !ROUTINE: init_floe_0
!
! !DESCRIPTION:
!
!  Initialize ice floe tracer
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
  subroutine init_floe_0
!
! !USES:
!
    use ice_state, only: trcrn!, nt_ifloe
    use icepack_tracers, only: nt_fsd, nfsd
    use ice_flux,  only: ifd
    use ice_blocks, only: nx_block, ny_block
    use ice_constants, only: c0, c1, c2, c3, c4, c10, c30, c110, c300, c1000, p25, p01
!
!EOP
!
  integer (kind=int_kind)            :: i, j
! This is broken - Noah day 20/3/22
!trcrn(:,:,nt_fsd,:,:) = c0
ifd(:,:,:) = c0
do j = 1, ny_block
  do i = 1, nx_block
    !trcrn(:,:,nt_fsd,1,:) = c300 ! Noah Day WIM, instead of :,: at the end
    ifd(i,j,:)            = c300
  enddo
enddo


end subroutine init_floe_0

!=======================================================================

!BOP
!
! !ROUTINE: increment_floe
!
! !DESCRIPTION:
!
!  Increase ice floe tracer by scaled timestep length.
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
        subroutine increment_floe (nx_block, ny_block,                   &
                                  dt, tmask, Lcell,                      &
                                  loc_swh, loc_ppd, loc_mwd,             &
                                  ifloe,  afice,  vfice, dum_wavemask, wave_spec_blk ) ! Noah Day 4/11/21

! !USES:
!
      use m_prams_waveice
      use m_waveice
      use m_waveattn
      use ice_fileunits, only: nu_diag

! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block       ! block dimensions

!    ! time step
      real (kind=dbl_kind), intent(in) :: &
         dt

!     ! ice concentration & thickness (1st category)
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         afice, vfice

!     ! ice floe size param, significant wave height, wave peak period
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
          loc_swh, loc_ppd, loc_mwd, ifloe

!     ! land/boundary mask, thickness (T-cell)
      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          Lcell

	   integer (kind=int_kind), intent(in)   :: dum_wavemask
!
!  local variables
!
!     ! Counters
      integer (kind=int_kind)            :: i, j, ij

!	  ! For the wave-ice code:
!     ! counter
      integer :: lp_i, lp_j
!     ! points in freq & angular domains
      integer, parameter          :: nw_in = 16, nth_in = 1 !ND: nw_in=31, nth_in=1
      !integer :: nw_in
!     ! max floe size param
      real(kind=dbl_kind)                :: D1
!     ! cell length
!      real(kind=dbl_kind), parameter     :: Lcell = 25*1000d0
!     ! open ocean wave params
      real(kind=dbl_kind), parameter                  :: fmin = 1d0/1000d0, fmax = 1d0 ! Noah Day changing from [1/16, 1/6], 1/s
                                                  ! freq min/max
      real(kind=dbl_kind), parameter                  :: om1=2*pi*fmin, om2=2*pi*fmax ! Angular frequencies, rad/s
                                                  ! ang freqs
      real(kind=dbl_kind), parameter                  :: om_0 = (om2 - om1)/(nw_in-1) ! rad/s
      real(kind=dbl_kind), dimension(nw_in)           :: om_in, T
                                                  ! freqency
      real(kind=dbl_kind), dimension(nth_in)          :: th_in
                                                  ! direction (not used yet!!!!)
 	  real(kind=dbl_kind), dimension(nw_in)           :: S_init_in
 	                                              ! initial energy spec

 	  real(kind=dbl_kind), dimension(nw_in)           :: lam_wtr_in, k_wtr_in
 	                                              ! wavelen/num

 	  real(kind=dbl_kind)                             :: dum_sm0, dum_sm2
 	  											  ! dummy spectral moments

 	  real(kind=dbl_kind), parameter                  :: ws_tol = 1.0e-1_dbl_kind

!     ! spectrum passed through ice-covered ocean
      real (kind=dbl_kind), dimension(nx_block,nw_in) :: &
         wspec_row

!     ! dummy version of wspec_row to help with directionality
      real (kind=dbl_kind), dimension(nx_block,nw_in) :: &
         wspec_row_hld

!     ! sine mean wave directions for row
      real (kind=dbl_kind), dimension(nx_block) :: &
         mwd_row

!     ! array to hold values for calculation of mean wave directions
      real (kind=dbl_kind), dimension(2,nx_block) :: &
         mwd_hld


 ! Noah Day 4/11/21
      real (kind=dbl_kind), dimension(nx_block,ny_block,nw_in), intent(inout) :: &
        wave_spec_blk ! Noah Day 4/11/21, wave spectrum for the block

      integer, dimension(nx_block)                    :: tmt, tmt_hld ! WIM termination flag

!     ! attenuation
	  integer              :: idd_alp
      real(kind=dbl_kind)  :: dum_alp

       ! Noah Day 2/2/24 Amplitude drop
      logical, dimension(nx_block, nw_in) :: flg_amp_drop_global

!EOP
!

!!! Define wave mask edge and propagate waves into ice:
     if (cmt.ne.0) then
      write(nu_diag,*) '    in increment_floe -> wavemask=', dum_wavemask
     endif

!!! Define attenuation coefficient coefficients

	 if (ATTEN_MODEL.eq.1) then
	  write(nu_diag,*) 'opening: ', waveicedatadir , fname_alp, '...'

	  open(newunit=idd_alp,file=waveicedatadir // fname_alp)
	  do lp_j=1,Ncheb_h+1
	   do lp_i=1,Ncheb_f+1
		read(idd_alp,*) dum_alp
		alp_coeffs(lp_i,lp_j)=dum_alp
	   end do
	  end do
	  close(idd_alp)
	  write(nu_diag,*) '... done'
	 endif


! Set frequency bin number equal to CICE's from ice_in
   ! nw_in = nfreq
   !!! Calculate wavenumbers & wavelengths

   	  do lp_i=1,nw_in
          om_in(lp_i)        = om1 + (lp_i-1)*om_0
          T(lp_i)            = 2d0*pi/om_in(lp_i)
          lam_wtr_in(lp_i)   = gravity*(T(lp_i)**2d0)/2d0/pi
          k_wtr_in(lp_i)     = 2d0*pi/lam_wtr_in(lp_i)
         end do

         if (nth_in.ne.1) then
          do lp_i=1,nth_in
           th_in(lp_i)        = -pi/2 + (lp_i-1)*pi/(nth_in-1)
          end do
         else
          th_in = 0d0
         end if

    	 tmt(:)             = 0      ! flag to terminate routine
    	 tmt_hld(:)         = 1      ! flag to terminate routine
       wspec_row(:,:)     = c0     ! a dummy vector
       wspec_row_hld(:,:) = c0     ! a dummy dummy vector
       mwd_hld(:,:)       = c0     ! another dummy dummy vector, noah day uncommented
       !loc_mwd(:,:)       = c0     ! another dummy dummy vector
      S_init_in(:)       = c0

       !!! Begin at wavemask (only difference is initialisation)
       if (cmt.ne.0) then
        write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
        write(nu_diag,*) '                      -> starting wavemask loop'
        write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
       endif

       j = dum_wavemask

        !     loc_mwd(:,j)    = pi/6d0   !    initialise, noah day uncommented

       ! A0. initialise with: Bretschneider
       do i=1,nx_block
        if (tmask(i,j)) then
        ! initialise direction
        mwd_row(i) = loc_mwd(i,j)
        if (loc_swh(i,j).lt.ws_tol) then
         if (cmt.ne.0) then
          write(nu_diag,*) '>>>--------------------------------------------->>>'
          write(nu_diag,*) '                         i, j = ', i, j
          write(nu_diag,*) '                         NOT running wave-ice routine'
          write(nu_diag,*) '                      -> sig wave ht <', ws_tol
          write(nu_diag,*) '<<<---------------------------------------------<<<'
         endif
         tmt(i) = 1
         if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
          ifloe(i,j) = floe_sz_pancake
          if (cmt.ne.0) then
           if (vfice(i,j).lt.tolh) then
            !write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
           endif
           if (afice(i,j).lt.tolice) then
            write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
           endif
           !write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
          endif ! END IF COMMENT
         endif  ! END IF h<tolh OR c<tolc
        else ! IF there are waves:
         do lp_i=1,nw_in
          S_init_in(lp_i) = SDF_Bretschneider(om_in(lp_i),0,loc_swh(i,j),loc_ppd(i,j))
          wave_spec_blk(i,j,lp_i) = S_init_in(lp_i) ! Noah Day 4/11/21
         end do
         ! want consistency between SWH definitions:
         !if (i.eq.1.and.cmt.ne.0) then

         if (cmt.ne.0) then
         if (loc_swh(i,j).gt.c0) then
            !write(nu_diag,*) '                      -> check: swh ', loc_swh(i,j)
            !write(nu_diag,*) '                      ->        ppd ', loc_ppd(i,j)
            dum_sm0        = fn_SpecMoment(S_init_in,nw_in,nth_in,om_in,th_in,0,nu_diag)
   	        dum_sm2        = fn_SpecMoment(S_init_in,nw_in,nth_in,om_in,th_in,2,nu_diag)
            !write(nu_diag,*) '                      ->        swh ', 4d0*(dum_sm0**0.5d0)
            !write(nu_diag,*) '                      ->        ppd ', &
            !     											2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
            !write(nu_diag,*) '     -> om_in ', om_in
            !write(nu_diag,*) '     -> S_init_in ', S_init_in
            !write(nu_diag,*) '     -> wave_spec_blk', wave_spec_blk(i,j,:)
        endif ! cmt
        endif ! END COMMENT

          ! A1. Use WIM to update wave spectrum and floe sizes

        if (cmt.ne.0) then
           write(nu_diag,*) '>>>--------------------------------------------->>>'
           write(nu_diag,*) '        INPUT         -> i,j   =', i, j
 		       write(nu_diag,*) '                      -> ifloe =', ifloe(i,j)
           write(nu_diag,*) '                      -> swh   =', loc_swh(i,j)
           write(nu_diag,*) '                      -> ppd   =', loc_ppd(i,j)
           write(nu_diag,*) '                      -> mwd   =', 180d0*mwd_row(i)/pi
    	   endif
    	if (do_coupled.ne.0) then
    	      call sub_Balance(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j), &
 		  nw_in,nth_in,om_in,th_in,k_wtr_in,S_init_in,wspec_row_hld(i,:),tmt(i),nu_diag)
 		 else
 		   call sub_Uncoupled(ifloe(i,j),Lcell(i,j),vfice(i,j),afice(i,j), &
 		   nw_in,nth_in,om_in,th_in,k_wtr_in,S_init_in,wspec_row_hld(i,:),tmt(i),nu_diag, flg_amp_drop_global(i,:))
 		 endif
 		! Noah Day 20/3/21 ifloe(i,j) = D1
 		endif ! ENDIF ws_tol
        else
         if (cmt.ne.0) then
          write(nu_diag,*) '>>>--------------------------------------------->>>'
          write(nu_diag,*) '                         i, j = ', i, j
          write(nu_diag,*) '                      -> LAND'
          write(nu_diag,*) '<<<---------------------------------------------<<<'
         endif
         tmt(i) = 1
        endif ! ENDIF tmask(i,j)
       enddo ! ENDDO i=1,nx_block

       if (cmt.ne.0) then
        write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
        write(nu_diag,*) '                      -> do wave directions', dum_wavemask
        write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
       endif

       !print*, 'doing wave directions on wavemask loop'

       ! A2. Redistribute energy according to mean directions
        ! LH Boundary:
         i=1
         ! if the cell isn't land ...
       if (tmask(i,j).and.tmt(i).ne.1) then
        ! Calculate a dummy swh to weight directions:
        dum_sm0        = &
           fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
 	   dum_sm0        = 4d0*(dum_sm0**0.5d0)
        if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
         !write(nu_diag,*) '                      -> going south', 1
         !print*, '                      -> going south', 1
         wspec_row(1,:)  = wspec_row(1,:) + &
        						(1d0-sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
         mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(1)*dum_sm0
         mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
         tmt_hld(1)      = 0
        endif ! ENDIF -pi/4<mwd<pi/4
        ! if wave energy needs to be advected to the east ...
        if (tmask(2,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
         !write(nu_diag,*) '                      -> going east', 1
         wspec_row(2,:) = wspec_row(2,:) + &
         				(sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
         mwd_hld(1,2)   = mwd_hld(1,2) + mwd_row(1)*dum_sm0
         mwd_hld(2,2)   = mwd_hld(2,2) + dum_sm0
         tmt_hld(2)     = 0
         ! if wave energy needs to be advected to the west (wrap around vector) ...
        elseif (tmask(nx_block,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
         !write(nu_diag,*) '                      -> going west', 1
         wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                                (sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
         mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + mwd_row(1)*dum_sm0
         mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
         tmt_hld(nx_block)     = 0
        endif
 	  endif ! END IF TMASK

  ! RH Boundary:
    i=nx_block
      ! if the cell isn't land ...
    if (tmask(i,j).and.tmt(i).ne.1) then
     dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
   dum_sm0        = 4d0*(dum_sm0**0.5d0)
     if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
      !write(nu_diag,*) '                      -> going south', nx_block
      wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                            (1d0-sin(2d0*mwd_row(nx_block))**2d0)*wspec_row_hld(nx_block,:)
      mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + &
                             mwd_row(nx_block)*dum_sm0
      mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
      tmt_hld(nx_block)     = 0
     endif   ! ENDIF -pi/4<mwd<pi/4
      ! if wave energy needs to be advected to the WEST ...
     if (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
      !write(nu_diag,*) '                         going west', nx_block
      wspec_row(i-1,:) = wspec_row(i-1,:) + &
     	                         (sin(2d0*mwd_row(i))**2)*wspec_row_hld(i,:)
      mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
      mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
      tmt_hld(i-1)     = 0
      ! if wave energy needs to be advected to the EAST ...
     elseif (tmask(1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
      !write(nu_diag,*) '                         going east', nx_block
      wspec_row(1,:)  = wspec_row(1,:) + &
                             (sin(2d0*mwd_row(nx_block))**2)*wspec_row(nx_block,:)
      mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(nx_block)*dum_sm0
      mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
      tmt_hld(1)      = 0
     endif ! ENDIF (tmask(nx_block-1,j).and.sinmwd_row(nx_block).gt.c0)
  endif ! END IF TMASK

    ! loop the inner cells...
    do i=2,nx_block-1
     if (tmask(i,j).and.tmt(i).ne.1) then
      dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
    dum_sm0        = 4d0*(dum_sm0**0.5d0)
      if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
       ! South
       !write(nu_diag,*) '                      -> going south', i
       wspec_row(i,:) = wspec_row(i,:) + &
       			(1d0-sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
       mwd_hld(1,i)   = mwd_hld(1,i) + mwd_row(i)*dum_sm0
       mwd_hld(2,i)   = mwd_hld(2,i) + dum_sm0
       tmt_hld(i)     = 0
      endif ! ENDIF -pi/4<mwd<pi/4
      ! East
      if (tmask(i+1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
       !write(nu_diag,*) '                      -> going east', i
       wspec_row(i+1,:) = wspec_row(i+1,:) + &
       					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
       mwd_hld(1,i+1)   = mwd_hld(1,i+1) + mwd_row(i)*dum_sm0
       mwd_hld(2,i+1)   = mwd_hld(2,i+1) + dum_sm0
       tmt_hld(i+1)     = 0
      elseif (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
       !write(nu_diag,*) '                      -> going west', i
       wspec_row(i-1,:) = wspec_row(i-1,:) + &
       				(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
       mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
       mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
       tmt_hld(i-1)     = 0
      end if
     endif ! END IF TMASK
  enddo

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> update mean values', dum_wavemask
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
    endif

    ! A3. Calculate mean parameters
    do i=1,nx_block

     if (tmask(i,j)) then
      if (mwd_hld(2,i).gt.c0) then
        mwd_row(i) = mwd_hld(1,i)/mwd_hld(2,i)
      else
        mwd_row(i) = c0
      endif
      tmt(i)         = tmt_hld(i)
      loc_mwd(i,j)   = mwd_row(i)
      if (tmt(i).eq.1) then
       loc_swh(i,j)   = 0d0
       loc_ppd(i,j)   = 0d0
      else
	 dum_sm0        = fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
	 dum_sm2        = fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,2,nu_diag)
	 loc_swh(i,j)   = 4d0*(dum_sm0**0.5d0)
	 loc_ppd(i,j)   = 2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
   wave_spec_blk(i,j,:) = wspec_row(i,:) ! Noah Day
	endif
	if (cmt.ne.0) then
       !write(nu_diag,*) '           OUTPUT     -> i,j   =', i, j
  	   !write(nu_diag,*) '                      -> tmt      =', tmt(i)
       !write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
       !write(nu_diag,*) '                      -> swh      =', loc_swh(i,j)
       !write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j)
       !write(nu_diag,*) '                      -> mwd      =', 180d0*loc_mwd(i,j)/pi
       !write(nu_diag,*) '                      -> wave_spec_blk(i,j,:)      =', wave_spec_blk(i,j,:)
       !write(nu_diag,*) '<<<---------------------------------------------<<<'
      endif
     endif ! ENDIF tmask(i,j)
    enddo ! ENDDO i=1,nblock

    wspec_row_hld(:,:) = c0     ! reset the dummy dummy vector
    mwd_hld(:,:)       = c0     ! ditto
    tmt_hld(:)         = 1

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> starting wave propagation'
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
    endif

    do j=dum_wavemask-1,1,-1
     if (cmt.ne.0) then
      write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
      write(nu_diag,*) '                      -> j=', j
      write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
     endif

     ! B1. Use WIM to update wave spectrum and floe sizes
     do i=1,nx_block
      !print*, 'Hi: i,j=', i, j
      if (tmt(i).eq.0) then
      !print*, '... tmt=0'
      if (tmask(i,j)) then
       if (loc_swh(i,j+1).lt.ws_tol) then
      !  if (cmt.ne.0) then
         write(nu_diag,*) '>>>--------------------------------------------->>>'
         write(nu_diag,*) '                         i, j = ', i, j
         write(nu_diag,*) '                         NOT running wave-ice routine'
         write(nu_diag,*) '                      -> sig wave ht <', ws_tol
      !  endif
        tmt(i) = 1
        if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
         ifloe(i,j) = floe_sz_pancake
         if (cmt.ne.0) then
          if (vfice(i,j).lt.tolh) then
           write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
          endif
          if (afice(i,j).lt.tolice) then
           write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
          endif
          !write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
         endif ! END IF COMMENT
        endif  ! END IF h<tolh OR c<tolc

       else ! IF there are waves:
       !print*, '... tmask=true'
      ! if (cmt.ne.0) then
        write(nu_diag,*) '>>>--------------------------------------------->>>'
        write(nu_diag,*) '       INPUT          -> i,j      =', i, j
        !print*, '       INPUT          -> i,j      =', i, j
	      write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
	      !print*, '                      -> ifloe    =', ifloe(i,j)
        write(nu_diag,*) '                      -> swh      =', loc_swh(i,j+1)
        !print*, '                      -> swh      =', loc_swh(i,j+1)
        write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j+1)
        !print*, '                      -> ppd      =', loc_ppd(i,j+1)
        write(nu_diag,*) '                      -> mwd      =', 180d0*mwd_row(i)/pi
        !print*, '                      -> mwd      =', 180d0*mwd_row(i)/pi
        write(nu_diag,*) '                      -> afice      =', afice(i,j)
      ! endif  ! ENDIF if (cmt.ne.0) then
       if (do_coupled.ne.0) then
         call sub_Balance(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j),nw_in,&
          nth_in,om_in,th_in,k_wtr_in,wspec_row(i,:),wspec_row_hld(i,:),tmt(i),nu_diag)
       else
         call sub_Uncoupled(ifloe(i,j),Lcell(i,j),vfice(i,j),afice(i,j),nw_in,&
           nth_in,om_in,th_in,k_wtr_in,wspec_row(i,:),wspec_row_hld(i,:),tmt(i),nu_diag, flg_amp_drop_global(i,:))
       endif ! ENDIF (do_coupled.ne.0)
       !print*, '... done wave-ice routine'
       !!ifloe(i,j) = D1
        !print*, '... set ifloe(i,j)=', ifloe(i,j), 'D1=', D1
       endif ! ENDIF ws_tol
       else
        if (cmt.ne.0) then
         write(nu_diag,*) '>>>--------------------------------------------->>>'
         write(nu_diag,*) '                         i, j = ', i, j
         write(nu_diag,*) '                      -> LAND'
         write(nu_diag,*) '<<<---------------------------------------------<<<'
        endif
       tmt(i)         = 1
      endif ! ENDIF tmask(i,j)
      else
       if (cmt.ne.0) then
        write(nu_diag,*) '>>>--------------------------------------------->>>'
        write(nu_diag,*) '                         i, j = ', i, j
        write(nu_diag,*) '                      -> tmt(i)=1'
       endif
       if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
        ifloe(i,j) = floe_sz_pancake
        if (cmt.ne.0) then
         if (vfice(i,j).lt.tolh) then
          write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
         endif
         if (afice(i,j).lt.tolice) then
          write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
         endif
         !write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
        endif ! END IF COMMENT
       endif  ! END IF h<tolh OR c<tolc
      endif ! ENDIF tmt(i)
     enddo ! ENDDO i=1,nblock

     wspec_row(:,:) = c0     ! reset dummy vector

     if (cmt.ne.0) then
      write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
      write(nu_diag,*) '                      -> do wave directions', j
      write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
     endif

     ! B2. Redistribute energy according to mean directions
      ! LH Boundary:
      i=1
       ! if the cell isn't land ...
     if (tmask(1,j).and.tmt(1).ne.1) then
      dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
    dum_sm0        = 4d0*(dum_sm0**0.5d0)
      if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
       wspec_row(1,:)  = wspec_row(1,:) + &
      				(1d0-sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
       mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(1)*dum_sm0
       mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
       tmt_hld(1)      = 0
      endif ! ENDIF -pi/4<mwd<pi/4
      ! if wave energy needs to be advected to the east ...
      if (tmask(2,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
       wspec_row(2,:) = wspec_row(2,:) + &
       				(sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
       mwd_hld(1,2)   = mwd_hld(1,2) + mwd_row(1)*dum_sm0
       mwd_hld(2,2)   = mwd_hld(2,2) + dum_sm0
       tmt_hld(2)     = 0
       ! if wave energy needs to be advected to the west (wrap around vector) ...
      elseif (tmask(nx_block,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
       wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                             (sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(1,:)
       mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + mwd_row(1)*dum_sm0
       mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
       tmt_hld(nx_block)     = 0
      endif ! ENDIF (tmask(2,j).and.sinmwd_row(1).lt.c0)
     endif ! END IF TMASK
     ! RH Boundary:
     i=nx_block
      ! if the cell isn't land ...
     if (tmask(nx_block,j).and.tmt(nx_block).ne.1) then
      dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
    dum_sm0        = 4d0*(dum_sm0**0.5d0)
      if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
       wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                            (1d0-sin(2d0*mwd_row(nx_block))**2d0)*wspec_row_hld(nx_block,:)
       mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + &
                             mwd_row(nx_block)*dum_sm0
       mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
       tmt_hld(nx_block)     = 0
      endif ! ENDIF -pi/4<mwd<pi/4
      ! if wave energy needs to be advected to the west ...
      if (tmask(nx_block-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
       wspec_row(nx_block-1,:) = wspec_row(nx_block-1,:) + &
     	                         (sin(2d0*mwd_row(nx_block))**2)*wspec_row_hld(nx_block,:)
       mwd_hld(1,nx_block-1)   = mwd_hld(1,nx_block-1) + &
     							  mwd_row(nx_block)*dum_sm0
       mwd_hld(2,nx_block-1)   = mwd_hld(1,nx_block-1) + dum_sm0
       tmt_hld(nx_block-1)     = 0
      ! if wave energy needs to be advected to the east ...
      elseif (tmask(1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
       wspec_row(1,:) = wspec_row(1,:) + &
       				(sin(2d0*mwd_row(nx_block))**2)*wspec_row_hld(nx_block,:)
       mwd_hld(1,1)   = mwd_hld(1,1) + mwd_row(nx_block)*dum_sm0
       mwd_hld(2,1)   = mwd_hld(2,1) + dum_sm0
       tmt_hld(1)     = 0
      endif ! ENDIF (tmask(nx_block-1,j).and.sinmwd_hld(nx_block).gt.c0)
     endif ! ENDIF TMASK
     ! loop the inner cells...
     do i=2,nx_block-1
      if (tmask(i,j).and.tmt(i).ne.1) then
       dum_sm0        = &
        fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
     dum_sm0        = 4d0*(dum_sm0**0.5d0)
       if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
       ! South
        wspec_row(i,:) = wspec_row(i,:) + &
       				(1d0-sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
        mwd_hld(1,i)   = mwd_hld(1,i) + mwd_row(i)*dum_sm0
        mwd_hld(2,i)   = mwd_hld(2,i) + dum_sm0
        tmt_hld(i)     = 0
       endif ! ENDIF -pi/4<mwd<pi/4
       ! East
       if (tmask(i+1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
        wspec_row(i+1,:) = wspec_row(i+1,:) + &
        					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
        mwd_hld(1,i+1)   = mwd_hld(1,i+1) + mwd_row(i)*dum_sm0
        mwd_hld(2,i+1)   = mwd_hld(2,i+1) + dum_sm0
        tmt_hld(i+1)     = 0
       elseif (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
        wspec_row(i-1,:) = wspec_row(i-1,:) + &
        					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
        mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
        mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
        tmt_hld(i-1)     = 0
       endif ! ENDIF (tmask(i+1,j).and.sinmwd_row(i).lt.c0)
      endif ! ENDIF TMASK
     end do ! ENDDO i=2,nx_block-1

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> update mean values', j
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
    endif

    ! B3. Calculate mean parameters
     do i=1,nx_block
      if (tmask(i,j)) then
       if (mwd_hld(2,i).gt.c0) then
        mwd_row(i) = mwd_hld(1,i)/mwd_hld(2,i)
       else
        mwd_row(i) = c0
       endif
       loc_mwd(i,j) = mwd_row(i)
       tmt(i)         = tmt_hld(i)
       if (tmt(i).eq.1) then
        loc_swh(i,j)   = 0d0
        loc_ppd(i,j)   = 0d0
       else
	  dum_sm0        = &
         fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
	  dum_sm2        = &
	   fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,2,nu_diag)
	  loc_swh(i,j)   = 4d0*(dum_sm0**0.5d0)
	  loc_ppd(i,j)   = 2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
    wave_spec_blk(i,j,:) = wspec_row(i,:) ! Noah Day
       endif

      ! if (cmt.ne.0) then
        write(nu_diag,*) '         OUTPUT       -> i,j      =', i, j
        write(nu_diag,*) '                      -> tmt      =', tmt(i)
        write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
        write(nu_diag,*) '                      -> swh      =', loc_swh(i,j)
        write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j)
        write(nu_diag,*) '                      -> mwd      =', 180d0*loc_mwd(i,j)/pi
        write(nu_diag,*) '                      -> wspec_row=', wspec_row(i,:)
        write(nu_diag,*) '                      -> wave_spec_blk(i,j,:)      =', wave_spec_blk(i,j,:)
        write(nu_diag,*) '<<<---------------------------------------------<<<'
     !  endif
       endif ! END IF tmask(i,j)
      enddo ! END do i=1,nx_block

     wspec_row_hld(:,:) = c0     ! reset the dummy dummy vector
     mwd_hld(:,:)       = c0     ! ditto
     tmt_hld(:)         = 1

    enddo ! END do j=dum_wavemask-1,2,-1

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> finished increment_floe'
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
    endif





      end subroutine increment_floe


  !=======================================================================

  !BOP
  !
  ! !ROUTINE: increment_floe_long
  !
  ! !DESCRIPTION:
  !
  !  Increase ice floe tracer by scaled timestep length along meridonal lines.
  !
  ! !REVISION HISTORY: same as module
  !
  ! !INTERFACE:
  !
          subroutine increment_floe_long (nx_block, ny_block,                   &
                                    dt, tmask, Lcell,                      &
                                    loc_swh, loc_ppd, loc_mwd,             &
                                    ifloe,  afice,  vfice, dum_wavemask, dum_wavemask_vec, &
                                    wave_spec_blk ) ! Noah Day 4/11/21

  ! !USES:
  !
        use m_prams_waveice
        use m_waveice
        use m_waveattn
        use ice_fileunits, only: nu_diag

  ! !INPUT/OUTPUT PARAMETERS:
  !
        integer (kind=int_kind), intent(in) :: &
           nx_block, ny_block       ! block dimensions

  !    ! time step
        real (kind=dbl_kind), intent(in) :: &
           dt

  !     ! ice concentration & thickness (1st category)
        real (kind=dbl_kind), dimension(nx_block,ny_block), &
           intent(in) :: &
           afice, vfice

  !     ! ice floe size param, significant wave height, wave peak period
        real (kind=dbl_kind), dimension(nx_block,ny_block), &
           intent(inout) :: &
            loc_swh, loc_ppd, loc_mwd

! Noah Day moving ifloe from inout to just in
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
          ifloe ! ice floe size param
  !     ! land/boundary mask, thickness (T-cell)
        logical (kind=log_kind), dimension (nx_block,ny_block), &
            intent(in) :: &
            tmask

        real (kind=dbl_kind), dimension (nx_block,ny_block), &
            intent(in) :: &
            Lcell
       integer (kind=int_kind), intent(in)   :: dum_wavemask
  	    integer (kind=int_kind), dimension (nx_block), intent(in)   :: dum_wavemask_vec

  !
  !  local variables
  !
  !     ! Counters
        integer (kind=int_kind)            :: i, j, ij, jj

        ! Wavemask
        integer (kind=int_kind)            :: max_wavemask ! maximum value of the wavemask vector

  !	  ! For the wave-ice code:
  !     ! counter
        integer :: lp_i, lp_j
  !     ! points in freq & angular domains
        integer, parameter          :: nw_in=31, nth_in=1
  !     ! max floe size param
        real(kind=dbl_kind)                :: D1
  !     ! cell length
  !      real(kind=dbl_kind), parameter     :: Lcell = 25*1000d0
  !     ! open ocean wave params
        real(kind=dbl_kind), parameter                  :: fmin = 1d0/1000d0, fmax = 1d0 ! ND: changing from 1/16, 1/6
                                                    ! freq min/max
        real(kind=dbl_kind), parameter                  :: om1=2*pi*fmin, om2=2*pi*fmax
                                                    ! ang freqs
        real(kind=dbl_kind), parameter                  :: om_0 = (om2 - om1)/(nw_in-1)
        real(kind=dbl_kind), dimension(nw_in)           :: om_in, T
                                                    ! freqency
        real(kind=dbl_kind), dimension(nth_in)          :: th_in
                                                    ! direction (not used yet!!!!)
   	  real(kind=dbl_kind), dimension(nw_in)           :: S_init_in
   	                                              ! initial energy spec

   	  real(kind=dbl_kind), dimension(nw_in)           :: lam_wtr_in, k_wtr_in
   	                                              ! wavelen/num

   	  real(kind=dbl_kind)                             :: dum_sm0, dum_sm2
   	  											  ! dummy spectral moments

   	  real(kind=dbl_kind), parameter                  :: ws_tol = 1.0e-11_dbl_kind ! ND: changing tolerance from 1.0e-1

  !     ! spectrum passed through ice-covered ocean
        real (kind=dbl_kind), dimension(nx_block,nw_in) :: &
           wspec_row

  !     ! dummy version of wspec_row to help with directionality
        real (kind=dbl_kind), dimension(nx_block,nw_in) :: &
           wspec_row_hld

  !     ! sine mean wave directions for row
        real (kind=dbl_kind), dimension(nx_block) :: &
           mwd_row

  !     ! array to hold values for calculation of mean wave directions
        real (kind=dbl_kind), dimension(2,nx_block) :: &
           mwd_hld


   ! Noah Day 4/11/21
        real (kind=dbl_kind), dimension(nx_block,ny_block,nw_in), intent(inout) :: &
          wave_spec_blk ! Noah Day 4/11/21, wave spectrum for the block m^2s/rad
        real(kind=dbl_kind)                :: int_d ! Integral of the directional spectrum


        integer, dimension(nx_block)                    :: tmt, tmt_hld ! WIM termination flag

   ! Noah Day 2/2/24 Amplitude drop
      logical, dimension(nx_block, nw_in) :: flg_amp_drop_global

  !     ! attenuation
  	  integer              :: idd_alp
      real(kind=dbl_kind)  :: dum_alp

  !EOP
  !   
      ! Initial flags for amplitude drop
      flg_amp_drop_global(:,:) = .false.

  !!! Define wave mask edge and propagate waves into ice:
       if (cmt.ne.0) then
        write(nu_diag,*) '    in increment_floe_long -> wavemask=', dum_wavemask_vec
       endif

  !!! Define attenuation coefficient coefficients

  	 if (ATTEN_MODEL.eq.1) then
  	  write(nu_diag,*) 'opening: ', waveicedatadir , fname_alp, '...'

  	  open(newunit=idd_alp,file=waveicedatadir // fname_alp)
  	  do lp_j=1,Ncheb_h+1
  	   do lp_i=1,Ncheb_f+1
  		read(idd_alp,*) dum_alp
  		alp_coeffs(lp_i,lp_j)=dum_alp
  	   end do
  	  end do
  	  close(idd_alp)
  	  write(nu_diag,*) '... done'
  	 endif



  !!! Calculate wavenumbers & wavelengths

  do lp_i=1,nw_in
      om_in(lp_i)        = om1 + (lp_i-1)*om_0
      T(lp_i)            = 2d0*pi/om_in(lp_i)
      lam_wtr_in(lp_i)   = gravity*(T(lp_i)**2d0)/2d0/pi
      k_wtr_in(lp_i)     = 2d0*pi/lam_wtr_in(lp_i)
   end do

   if (nth_in.ne.1) then
      do lp_i=1,nth_in
         th_in(lp_i)        = -pi/2 + (lp_i-1)*pi/(nth_in-1)
      end do
   else
      th_in = 0d0
   end if

	tmt(:)             = 0      ! flag to terminate routine
	tmt_hld(:)         = 1      ! flag to terminate routine
   wspec_row(:,:)     = c0     ! a dummy vector
   wspec_row_hld(:,:) = c0     ! a dummy dummy vector
   mwd_hld(:,:)       = c0     ! another dummy dummy vector, noah day uncommented
   !loc_mwd(:,:)       = c0     ! another dummy dummy vector
	 S_init_in(:)       = c0   ! Wave energy spectrum, m^2s/rad

   !!! Begin at wavemask (only difference is initialisation)
   if (cmt.ne.0) then
    write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
    write(nu_diag,*) '                      -> starting wavemask loop'
    write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
   endif


    !     loc_mwd(:,j)    = pi/6d0   !    initialise, noah day uncommented

   

   ! A0. initialise with: Bretschneider
   do i=1,nx_block
     j = dum_wavemask_vec(i)
     if (j.gt.0) then
        !write(nu_diag,*) 'i is:', i
         !write(nu_diag,*) 'j is:', j
         if (tmask(i,j)) then
           ! initialise direction
           mwd_row(i) = loc_mwd(i,j)
           if (loc_swh(i,j).lt.ws_tol) then
             if (cmt.ne.0) then
                write(nu_diag,*) '>>>--------------------------------------------->>>'
                write(nu_diag,*) '                         i, j = ', i, j
                write(nu_diag,*) '                         NOT running wave-ice routine'
                write(nu_diag,*) '                      -> sig wave ht <', ws_tol
                write(nu_diag,*) '<<<---------------------------------------------<<<'
             endif
             tmt(i) = 1
             if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
               !ifloe(i,j) = floe_sz_pancake
               if (cmt.ne.0) then
                 if (vfice(i,j).lt.tolh) then
                  !write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
                 endif
                 if (afice(i,j).lt.tolice) then
                  !write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
                 endif
                 write(nu_diag,*) ' For cell (i,j)', i,j
                 write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
               endif ! END IF COMMENT
             endif  ! END IF h<tolh OR c<tolc
           else ! IF there are waves:
             do lp_i=1,nw_in
               ! ND: Initiailise the spectrum, S(omega)
               S_init_in(lp_i) = SDF_Bretschneider(om_in(lp_i),0,loc_swh(i,j),loc_ppd(i,j)) ! m^2s/rad
             end do
             ! ND: Initialise and integrate the directional spectrum, D(theta)
             int_d = dir_spec_integral(S_init_in, loc_mwd(i,j), nu_diag) ! if WIM_DIR = 0 then int_d=1 i.e. all energy is taken

             do lp_i=1,nw_in
               ! ND: E(omega, theta) = S(omega)D(theta)
               wave_spec_blk(i,j,lp_i) = S_init_in(lp_i)*int_d
             enddo
             ! want consistency between SWH definitions:
             !if (i.eq.1.and.cmt.ne.0) then

             if (cmt.ne.0) then
             if (loc_swh(i,j).gt.c0) then
                write(nu_diag,*) '                      -> check: swh ', loc_swh(i,j)
                write(nu_diag,*) '                      ->        ppd ', loc_ppd(i,j)
                  dum_sm0        = fn_SpecMoment(S_init_in,nw_in,nth_in,om_in,th_in,0,nu_diag)
       	          dum_sm2        = fn_SpecMoment(S_init_in,nw_in,nth_in,om_in,th_in,2,nu_diag)
                write(nu_diag,*) '                      ->        swh ', 4d0*(dum_sm0**0.5d0)
                write(nu_diag,*) '                      ->        ppd ', &
                       											2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
                write(nu_diag,*) '     -> om_in ', om_in
                write(nu_diag,*) '     -> S_init_in ', S_init_in
                write(nu_diag,*) '     -> wave_spec_blk', wave_spec_blk(i,j,:)
             endif ! cmt
             endif ! END COMMENT

            ! A1. Use WIM to update wave spectrum and floe sizes

             if (cmt.ne.0) then
              write(nu_diag,*) '>>>--------------------------------------------->>>'
              write(nu_diag,*) '        INPUT A1        -> i,j   =', i, j
      	       write(nu_diag,*) '                      -> ifloe =', ifloe(i,j)
              write(nu_diag,*) '                      -> swh   =', loc_swh(i,j)
              write(nu_diag,*) '                      -> ppd   =', loc_ppd(i,j)
              write(nu_diag,*) '                      -> mwd   =', 180d0*mwd_row(i)/pi
      	     endif

             if (do_coupled.ne.0) then
          	      call sub_Balance(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j), &
       		  nw_in,nth_in,om_in,th_in,k_wtr_in,S_init_in,wspec_row_hld(i,:),tmt(i),nu_diag)
            !write(nu_diag,*) '>>>--------------------------------------------->>>'
            !write(nu_diag,*) 'finished sub_balance'
       		   else
       		        call sub_Uncoupled(ifloe(i,j),Lcell(i,j),vfice(i,j),afice(i,j), &
       		   nw_in,nth_in,om_in,th_in,k_wtr_in,S_init_in,wspec_row_hld(i,:),tmt(i),nu_diag, flg_amp_drop_global(i,:))
       		   endif
       		   !ifloe(i,j) = D1
       		endif ! ENDIF ws_tol
         else ! Else TMASK
               if (cmt.ne.0) then
                write(nu_diag,*) '>>>--------------------------------------------->>>'
                write(nu_diag,*) '                         i, j = ', i, j
                write(nu_diag,*) '                      -> LAND'
                write(nu_diag,*) '<<<---------------------------------------------<<<'
               endif
               tmt(i) = 1
        endif ! ENDIF tmask(i,j)
      endif ! j > 0
  enddo ! ENDDO i=1,nx_block

 if (cmt.ne.0) then
  write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
  write(nu_diag,*) '                      -> do wave directions'
  write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
 endif

 !print*, 'doing wave directions on wavemask loop'

! ND: end of loop ------------------------------------------------------------------------------------------




 ! A2. Redistribute energy according to mean directions
  ! LH Boundary:
   i=1
   j = dum_wavemask_vec(i)
   !write(nu_diag,*) 'j is: ', j
     ! if the cell isn't land ...
   if (tmask(i,j).and.tmt(i).ne.1) then
        ! Calculate a dummy swh to weight directions:
          dum_sm0        = &
           fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
          dum_sm0        = 4d0*(dum_sm0**0.5d0)
          !if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
           !write(nu_diag,*) '                      -> going south', 1
           !print*, '                      -> going south', 1
           wspec_row(1,:)  = wspec_row(1,:) + &
          						(1d0-sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
           mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(1)*dum_sm0
           mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
           tmt_hld(1)      = 0
          !endif ! ENDIF -pi/4<mwd<pi/4
 !         ! if wave energy needs to be advected to the east ...
 !         if (tmask(2,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
 !          !write(nu_diag,*) '                      -> going east', 1
 !          wspec_row(2,:) = wspec_row(2,:) + &
 !          				(sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
 !          mwd_hld(1,2)   = mwd_hld(1,2) + mwd_row(1)*dum_sm0
 !          mwd_hld(2,2)   = mwd_hld(2,2) + dum_sm0
 !          tmt_hld(2)     = 0
 !          ! if wave energy needs to be advected to the west (wrap around vector) ...
 !         elseif (tmask(nx_block,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
 !          !write(nu_diag,*) '                      -> going west', 1
 !          wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
 !                                 (sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
 !          mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + mwd_row(1)*dum_sm0
 !          mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
 !          tmt_hld(nx_block)     = 0
 !         endif
   	endif ! END IF TMASK

    ! RH Boundary:
      i=nx_block
      j = dum_wavemask_vec(i)
        ! if the cell isn't land ...
      if (tmask(i,j).and.tmt(i).ne.1) then
       dum_sm0        = &
          fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
     dum_sm0        = 4d0*(dum_sm0**0.5d0)
 !          if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
            !write(nu_diag,*) '                      -> going south', nx_block
            wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                                  (1d0-sin(2d0*mwd_row(nx_block))**2d0)*wspec_row_hld(nx_block,:)
            mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + &
                                   mwd_row(nx_block)*dum_sm0
            mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
            tmt_hld(nx_block)     = 0
 !          endif   ! ENDIF -pi/4<mwd<pi/4
 !           ! if wave energy needs to be advected to the WEST ...
 !          if (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
 !           !write(nu_diag,*) '                         going west', nx_block
 !           wspec_row(i-1,:) = wspec_row(i-1,:) + &
 !          	                         (sin(2d0*mwd_row(i))**2)*wspec_row_hld(i,:)
 !           mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
 !           mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
 !           tmt_hld(i-1)     = 0
 !           ! if wave energy needs to be advected to the EAST ...
 !          elseif (tmask(1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
 !           !write(nu_diag,*) '                         going east', nx_block
 !           wspec_row(1,:)  = wspec_row(1,:) + &
 !                                  (sin(2d0*mwd_row(nx_block))**2)*wspec_row(nx_block,:)
 !           mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(nx_block)*dum_sm0
 !           mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
 !           tmt_hld(1)      = 0
 !          endif ! ENDIF (tmask(nx_block-1,j).and.sinmwd_row(nx_block).gt.c0)
    endif ! END IF TMASK

    ! loop the inner cells...
    do i=2,nx_block-1
      j = dum_wavemask_vec(i)
       if (tmask(i,j).and.tmt(i).ne.1) then
        dum_sm0        = &
          fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
      dum_sm0        = 4d0*(dum_sm0**0.5d0)
 !           if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
             ! South
             !write(nu_diag,*) '                      -> going south', i
             wspec_row(i,:) = wspec_row(i,:) + &
             			(1d0-sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
             mwd_hld(1,i)   = mwd_hld(1,i) + mwd_row(i)*dum_sm0
             mwd_hld(2,i)   = mwd_hld(2,i) + dum_sm0
             tmt_hld(i)     = 0
 !           endif ! ENDIF -pi/4<mwd<pi/4
 !           ! East
 !           if (tmask(i+1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
 !            !write(nu_diag,*) '                      -> going east', i
 !            wspec_row(i+1,:) = wspec_row(i+1,:) + &
 !            					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
 !            mwd_hld(1,i+1)   = mwd_hld(1,i+1) + mwd_row(i)*dum_sm0
 !            mwd_hld(2,i+1)   = mwd_hld(2,i+1) + dum_sm0
 !            tmt_hld(i+1)     = 0
 !           elseif (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
 !            !write(nu_diag,*) '                      -> going west', i
 !            wspec_row(i-1,:) = wspec_row(i-1,:) + &
 !            				(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
 !            mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
 !            mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
 !            tmt_hld(i-1)     = 0
 !           end if
       endif ! END IF TMASK
    enddo ! i = nx_block

    if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> update mean values', dum_wavemask_vec
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
    endif
! ND: end of loop ------------------------------------------------------------------------------------------




    ! A3. Calculate mean parameters
  do i=1,nx_block
      j = dum_wavemask_vec(i)
        if (tmask(i,j)) then
            if (mwd_hld(2,i).gt.c0) then
              mwd_row(i) = mwd_hld(1,i)/mwd_hld(2,i)
            else
              mwd_row(i) = c0
            endif
            tmt(i)         = tmt_hld(i)
            loc_mwd(i,j)   = mwd_row(i)
            if (tmt(i).eq.1) then
               loc_swh(i,j)   = 0d0
               loc_ppd(i,j)   = 0d0
            else
            	 dum_sm0        = fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
            	 dum_sm2        = fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,2,nu_diag)
            	 loc_swh(i,j)   = 4d0*(dum_sm0**0.5d0)
            	 loc_ppd(i,j)   = 2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
               wave_spec_blk(i,j,:) = wspec_row(i,:) ! Noah Day
           endif
        	if (cmt.ne.0) then
                 write(nu_diag,*) '           OUTPUT A3    -> i,j   =', i, j
            	   write(nu_diag,*) '                      -> tmt      =', tmt(i)
                 write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
                 write(nu_diag,*) '                      -> swh      =', loc_swh(i,j)
                 write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j)
                 write(nu_diag,*) '                      -> mwd      =', 180d0*loc_mwd(i,j)/pi
                 write(nu_diag,*) '                      -> wave_spec_blk(i,j,:)      =', wave_spec_blk(i,j,:)
                 write(nu_diag,*) '<<<---------------------------------------------<<<'
            endif
       endif ! ENDIF tmask(i,j)
  enddo ! ENDDO i=1,nblock

! ND: end of loop ------------------------------------------------------------------------------------------




  wspec_row_hld(:,:) = c0     ! reset the dummy dummy vector
  mwd_hld(:,:)       = c0     ! ditto
  tmt_hld(:)         = 1

  if (cmt.ne.0) then
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) '                      -> starting wave propagation'
     write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
     write(nu_diag,*) 'MINIMUM: ', minval(dum_wavemask_vec)
     write(nu_diag,*) 'max dum_wavemask is:', maxval(dum_wavemask_vec)
  endif

!  if (maxval(dum_wavemask_vec).gt.0) then
!    ! There is ice
!    max_wavemask = maxval(dum_wavemask_vec)
!  else
!    ! No ice
!    max_wavemask = 50
!  end if ! maxval

max_wavemask = dum_wavemask
!write(nu_diag,*) 'max dum_wavemask is:', max_wavemask
  do jj=1,max_wavemask!maxval(dum_wavemask_vec)!dum_wavemask-1,2,minval(dum_wavemask_vec)+1!-1
       if (cmt.ne.0) then
          write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
          write(nu_diag,*) '                      -> jj=', jj
          write(nu_diag,*) 'dum_wavemask is:', dum_wavemask

          write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooooo'
       endif

       ! B1. Use WIM to update wave spectrum and floe sizes
       do i=1,nx_block
         j = dum_wavemask_vec(i) - jj!(dum_wavemask_vec(i) - jj) ! wavemask - (wavemask-j) = j
         !write(nu_diag,*) '>>>--------------------------------------------->>>'
         !write(nu_diag,*) 'B1. j is:', j
         !write(nu_diag,*) ' i is:', i
         !write(nu_diag,*) 'dum_wavemask_vec is:', dum_wavemask_vec(i)
         !write(nu_diag,*) 'jj is:', jj
         !write(nu_diag,*) '<<<---------------------------------------------<<<'
         if (j.gt.0) then
            !print*, 'Hi: i,j=', i, j
            if (tmt(i).eq.0) then
                  !print*, '... tmt=0'
                  if (tmask(i,j)) then
                       if (loc_swh(i,j+1).lt.ws_tol) then
                          if (cmt.ne.0) then
                             write(nu_diag,*) '>>>--------------------------------------------->>>'
                             write(nu_diag,*) '                         i, j = ', i, j
                             write(nu_diag,*) '                         NOT running wave-ice routine'
                             write(nu_diag,*) '                      -> sig wave ht <', ws_tol
                          endif
                          tmt(i) = 1
                          if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
                              !ifloe(i,j) = floe_sz_pancake
                                  if (cmt.ne.0) then
                                      if (vfice(i,j).lt.tolh) then
                                        write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
                                        write(nu_diag,*) '                       -> coords', i,j
                                      endif
                                      if (afice(i,j).lt.tolice) then
                                        write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
                                        write(nu_diag,*) '                       -> coords', i,j
                                      endif
                                      !write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
                                  endif ! END IF COMMENT
                          endif  ! END IF h<tolh OR c<tolc
                      else ! IF there are waves:
                         !print*, '... tmask=true'
                         if (cmt.ne.0) then
                              write(nu_diag,*) '>>>--------------------------------------------->>>'
                              write(nu_diag,*) '       INPUT  B1        -> i,j      =', i, j
                              !print*, '       INPUT          -> i,j      =', i, j
                      	      write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
                      	      !print*, '                      -> ifloe    =', ifloe(i,j)
                              write(nu_diag,*) '                      -> swh      =', loc_swh(i,j+1)
                              !print*, '                      -> swh      =', loc_swh(i,j+1)
                              write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j+1)
                              !print*, '                      -> ppd      =', loc_ppd(i,j+1)
                              write(nu_diag,*) '                      -> mwd      =', 180d0*mwd_row(i)/pi
                              !print*, '                      -> mwd      =', 180d0*mwd_row(i)/pi
                          endif  ! ENDIF if (cmt.ne.0) then
                       if (do_coupled.ne.0) then
                         call sub_Balance(ifloe(i,j),D1,Lcell(i,j),vfice(i,j),afice(i,j),nw_in,&
                          nth_in,om_in,th_in,k_wtr_in,wspec_row(i,:),wspec_row_hld(i,:),tmt(i),nu_diag)
                       else
                         call sub_Uncoupled(ifloe(i,j),Lcell(i,j),vfice(i,j),afice(i,j),nw_in,&
                           nth_in,om_in,th_in,k_wtr_in,wspec_row(i,:),wspec_row_hld(i,:),tmt(i),nu_diag, flg_amp_drop_global(i,:))
                       endif ! ENDIF (do_coupled.ne.0)
                       !print*, '... done wave-ice routine'
                       ! Noah Day 20/3/22 ifloe(i,j) = D1
                        !print*, '... set ifloe(i,j)=', ifloe(i,j), 'D1=', D1
                     endif ! ENDIF ws_tol
                   else ! tmask
                    if (cmt.ne.0) then
                       write(nu_diag,*) '>>>--------------------------------------------->>>'
                       write(nu_diag,*) '                         i, j = ', i, j
                       write(nu_diag,*) '                      -> LAND'
                       write(nu_diag,*) '<<<---------------------------------------------<<<'
                    endif
                    tmt(i)         = 1
                  endif ! ENDIF tmask(i,j)
            else ! tmt
                   if (cmt.ne.0) then
                      write(nu_diag,*) '>>>--------------------------------------------->>>'
                      write(nu_diag,*) '                         i, j = ', i, j
                      write(nu_diag,*) '                      -> tmt(i)=1'
                   endif
                   if (afice(i,j).lt.tolice.or.vfice(i,j).lt.tolh) then
                      !ifloe(i,j) = floe_sz_pancake
                      if (cmt.ne.0) then
                         if (vfice(i,j).lt.tolh) then
                            write(nu_diag,*) '                       -> no ice in this cell: h<', tolh
                         endif
                         if (afice(i,j).lt.tolice) then
                            write(nu_diag,*) '                       -> no ice in this cell: c<', tolice
                         endif
                         write(nu_diag,*) '                       -> ifloe(i,j)=', ifloe(i,j)
                      endif ! END IF COMMENT
                   endif  ! END IF h<tolh OR c<tolc
            endif ! ENDIF tmt(i)
         endif ! j>0
       enddo ! ENDDO i=1,nblock

       wspec_row(:,:) = c0     ! reset dummy vector

       if (cmt.ne.0) then
          write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
          write(nu_diag,*) '                      -> do wave directions 2'!, j
          write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
       endif

       ! B2. Redistribute energy according to mean directions
       ! LH Boundary:
       i=1
       j = dum_wavemask_vec(i) - jj!dum_wavemask_vec(i) - (dum_wavemask_vec(i) - jj) ! wavemask - (wavemask-j) = j
       !write(nu_diag,*) ' Redistribute energy according to mean directions '
       !write(nu_diag,*) ' LH Boundary: j: ', j
       ! if the cell isn't land ...
       if (j.gt.0) then
           if (tmask(1,j).and.tmt(1).ne.1) then
               dum_sm0        = &
               fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
               dum_sm0        = 4d0*(dum_sm0**0.5d0)
!              if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
               ! South
               wspec_row(1,:)  = wspec_row(1,:) + &
              				(1d0-sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
               mwd_hld(1,1)    = mwd_hld(1,1) + mwd_row(1)*dum_sm0
               mwd_hld(2,1)    = mwd_hld(2,1) + dum_sm0
               tmt_hld(1)      = 0
!              endif ! ENDIF -pi/4<mwd<pi/4
!              ! if wave energy needs to be advected to the east ...
!              if (tmask(2,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
!               wspec_row(2,:) = wspec_row(2,:) + &
!               				(sin(2d0*mwd_row(1))**2d0)*wspec_row_hld(1,:)
!               mwd_hld(1,2)   = mwd_hld(1,2) + mwd_row(1)*dum_sm0
!               mwd_hld(2,2)   = mwd_hld(2,2) + dum_sm0
!               tmt_hld(2)     = 0
!               ! if wave energy needs to be advected to the west (wrap around vector) ...
!              elseif (tmask(nx_block,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
!               wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
!                                     (sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(1,:)
!               mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + mwd_row(1)*dum_sm0
!               mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
!               tmt_hld(nx_block)     = 0
              !endif ! ENDIF (tmask(2,j).and.sinmwd_row(1).lt.c0)
           endif ! END IF TMASK
        endif ! j > 0
       ! RH Boundary:
       i=nx_block
       j = dum_wavemask_vec(i) - jj!dum_wavemask_vec(i) - (dum_wavemask_vec(i) - jj) ! wavemask - (wavemask-j) = j
        ! if the cell isn't land ...
        if (j.gt.0) then
           if (tmask(nx_block,j).and.tmt(nx_block).ne.1) then
                dum_sm0        = &
                  fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
              dum_sm0        = 4d0*(dum_sm0**0.5d0)
!                if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
                     ! South
                     wspec_row(nx_block,:) = wspec_row(nx_block,:) + &
                                          (1d0-sin(2d0*mwd_row(nx_block))**2d0)*wspec_row_hld(nx_block,:)
                     mwd_hld(1,nx_block)   = mwd_hld(1,nx_block) + &
                                           mwd_row(nx_block)*dum_sm0
                     mwd_hld(2,nx_block)   = mwd_hld(2,nx_block) + dum_sm0
                     tmt_hld(nx_block)     = 0
!                endif ! ENDIF -pi/4<mwd<pi/4
!                ! if wave energy needs to be advected to the west ...
!                if (tmask(nx_block-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
!                     wspec_row(nx_block-1,:) = wspec_row(nx_block-1,:) + &
!                   	                         (sin(2d0*mwd_row(nx_block))**2)*wspec_row_hld(nx_block,:)
!                     mwd_hld(1,nx_block-1)   = mwd_hld(1,nx_block-1) + &
!                   							  mwd_row(nx_block)*dum_sm0
!                     mwd_hld(2,nx_block-1)   = mwd_hld(1,nx_block-1) + dum_sm0
!                     tmt_hld(nx_block-1)     = 0
!                    ! if wave energy needs to be advected to the east ...
!                elseif (tmask(1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
!                     wspec_row(1,:) = wspec_row(1,:) + &
!                     				(sin(2d0*mwd_row(nx_block))**2)*wspec_row_hld(nx_block,:)
!                     mwd_hld(1,1)   = mwd_hld(1,1) + mwd_row(nx_block)*dum_sm0
!                     mwd_hld(2,1)   = mwd_hld(2,1) + dum_sm0
!                     tmt_hld(1)     = 0
               !endif ! ENDIF (tmask(nx_block-1,j).and.sinmwd_hld(nx_block).gt.c0)
           endif ! ENDIF TMASK
        endif ! j > 0

       ! loop the inner cells...
       do i=2,nx_block-1
         j = dum_wavemask_vec(i) - jj!dum_wavemask_vec(i) - (dum_wavemask_vec(i) - jj) ! wavemask - (wavemask-j) = j
         if (j.gt.0) then
              if (tmask(i,j).and.tmt(i).ne.1) then
                     dum_sm0        = &
                      fn_SpecMoment(wspec_row_hld(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
                    dum_sm0        = 4d0*(dum_sm0**0.5d0)
!                     if (mwd_row(i).gt.3d0*pi/4d0.and.mwd_row(i).lt.5d0*pi/4d0) then
                         ! South
                          wspec_row(i,:) = wspec_row(i,:) + &
                         				(1d0-sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
                          mwd_hld(1,i)   = mwd_hld(1,i) + mwd_row(i)*dum_sm0
                          mwd_hld(2,i)   = mwd_hld(2,i) + dum_sm0
                          tmt_hld(i)     = 0
!                     endif ! ENDIF -pi/4<mwd<pi/4
!                     ! East
!                     if (tmask(i+1,j).and.mwd_row(i).gt.pi/2d0.and.mwd_row(i).lt.pi) then
!                          wspec_row(i+1,:) = wspec_row(i+1,:) + &
!                          					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
!                          mwd_hld(1,i+1)   = mwd_hld(1,i+1) + mwd_row(i)*dum_sm0
!                          mwd_hld(2,i+1)   = mwd_hld(2,i+1) + dum_sm0
!                          tmt_hld(i+1)     = 0
!                     elseif (tmask(i-1,j).and.mwd_row(i).gt.pi.and.mwd_row(i).lt.3*pi/2d0) then
!                          wspec_row(i-1,:) = wspec_row(i-1,:) + &
!                          					(sin(2d0*mwd_row(i))**2d0)*wspec_row_hld(i,:)
!                          mwd_hld(1,i-1)   = mwd_hld(1,i-1) + mwd_row(i)*dum_sm0
!                          mwd_hld(2,i-1)   = mwd_hld(2,i-1) + dum_sm0
!                          tmt_hld(i-1)     = 0
!                     endif ! ENDIF (tmask(i+1,j).and.sinmwd_row(i).lt.c0)
              endif ! ENDIF TMASK
          endif !j > 0
       end do ! ENDDO i=2,nx_block-1

      if (cmt.ne.0) then
       write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
       write(nu_diag,*) '                      -> update mean values'!, j
       write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooo'
      endif

      ! B3. Calculate mean parameters
       do i=1,nx_block
          j = dum_wavemask_vec(i) - jj!dum_wavemask_vec(i) - (dum_wavemask_vec(i) - jj) ! wavemask - (wavemask-j) = j
          if (j.gt.0) then
              if (tmask(i,j)) then
                   if (mwd_hld(2,i).gt.c0) then
                      mwd_row(i) = mwd_hld(1,i)/mwd_hld(2,i)
                  else
                    mwd_row(i) = c0
                  endif
                 loc_mwd(i,j) = mwd_row(i)
                 tmt(i)         = tmt_hld(i)
                 if (tmt(i).eq.1) then
                    loc_swh(i,j)   = 0d0
                    loc_ppd(i,j)   = 0d0
                 else
                	  dum_sm0        = &
                         fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,0,nu_diag)
                	  dum_sm2        = &
                	   fn_SpecMoment(wspec_row(i,:),nw_in,nth_in,om_in,th_in,2,nu_diag)
                	  loc_swh(i,j)   = 4d0*(dum_sm0**0.5d0)
                	  loc_ppd(i,j)   = 2d0*pi*((dum_sm0/dum_sm2)**0.5d0)
                    wave_spec_blk(i,j,:) = wspec_row(i,:) ! Noah Day
                endif

               if (cmt.ne.0) then
                    write(nu_diag,*) '         OUTPUT  B3     -> i,j      =', i, j
                    write(nu_diag,*) '                      -> tmt      =', tmt(i)
                    write(nu_diag,*) '                      -> ifloe    =', ifloe(i,j)
                    write(nu_diag,*) '                      -> swh      =', loc_swh(i,j)
                    write(nu_diag,*) '                      -> ppd      =', loc_ppd(i,j)
                    write(nu_diag,*) '                      -> mwd      =', 180d0*loc_mwd(i,j)/pi
                    write(nu_diag,*) '                      -> wspec_row=', wspec_row(i,:)
                    write(nu_diag,*) '                      -> wave_spec_blk(i,j,:)      =', wave_spec_blk(i,j,:)
                    write(nu_diag,*) '<<<---------------------------------------------<<<'
               endif
              endif ! END IF tmask(i,j)
          endif ! j > 0
       enddo ! END do i=1,nx_block

       wspec_row_hld(:,:) = c0     ! reset the dummy dummy vector
       mwd_hld(:,:)       = c0     ! ditto
       tmt_hld(:)         = 1

enddo ! END do j = dum_wavemask-1,2,-1
! ND: end of loop ------------------------------------------------------------------------------------------





if (cmt.ne.0) then
   write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
   write(nu_diag,*) '                      -> finished increment_floe_long'
   write(nu_diag,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'
endif





        end subroutine increment_floe_long

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_floe - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_floe(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for restarting
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: msec, mmonth, mday, myear, istep1, &
                              timesecs, idate, year_init ! Noah Day WIM, removed , time_forc
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = myear + year_init - 1
         imonth = mmonth
         iday = mday

         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.floe.', &
              iyear,'-',mmonth,'-',mday,'-',msec
      end if

      ! begin writing restart data
      !call ice_open(nu_dump_floe,filename,0)

      if (my_task == master_task) then
        !write(nu_dump_floe) istep1,timesecs!,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         !call ice_write(nu_dump_floe,0,trcrn(:,:,nt_fsd,n,:),'ruf8',diag)
      enddo

      !if (my_task == master_task) close(nu_dump_floe)

      end subroutine write_restart_floe

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_floe - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_floe(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for an ice floe restart
!
! !REVISION HISTORY:
!
! author S O'Farrell
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: msec, mmonth, mday, myear, istep1, &
                              timesecs, idate, year_init ! Noah Day WIM, removed , time_forc
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag

      if (my_task == master_task) then
         open(nu_rst_pointer,file=pointer_file)
         read(nu_rst_pointer,'(a)') filename0
         filename = trim(filename0)
         close(nu_rst_pointer)

         ! reconstruct path/file
         n = index(filename0,trim(restart_file))
         if (n == 0) call abort_ice('ifloe restart: filename discrepancy')
         string1 = trim(filename0(1:n-1))
         string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
         write(filename,'(a,a,a,a)') &
            string1(1:lenstr(string1)), &
            restart_file(1:lenstr(restart_file)),'.floe.', &
            string2(1:lenstr(string2))
      endif ! master_task

      !call ice_open(nu_restart_floe,filename,0)

      if (my_task == master_task) then
        !read(nu_restart_floe) istep1,timesecs!,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do n = 1, ncat
         !call ice_read(nu_restart_floe,0,trcrn(:,:,nt_fsd,n,:),'ruf8',diag)
      enddo

      !if (my_task == master_task) close(nu_restart_floe)

      end subroutine read_restart_floe


end module ice_floe

!=======================================================================
