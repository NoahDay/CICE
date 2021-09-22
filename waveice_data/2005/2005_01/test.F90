program read_ww3
  use netcdf     
  implicit none
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

  ! !INPUT/OUTPUT PARAMETERS:
  !
        integer :: &
           mth ! month 1-12

        integer :: &
           N_tm, N_lat, N_lon

        ! WW3 variables
  	  integer                               :: ncid, varid, dimid, numDims
  	  !integer, dimension(nf90_max_var_dims) :: rhDimIds

  ! This is a comment line; it is ignored by the compiler
  print *, waveicedatadir

        mth = 1



       print *, '    sub_WW3_dataread WAVE_METH=1', mth
       !print*, '    sub_WW3_dataread WAVE_METH=1', mth
        if (mth.eq.1) then
         call check( nf90_open(waveicedatadir // fname_ww3 // 'dir_20050101.txt', NF90_NOWRITE, ncid) )
        elseif (mth.eq.2) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '02_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.3) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '03_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.4) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '04_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.5) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '05_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.6) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '06_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.7) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '07_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.8) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '08_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.9) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '09_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.10) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '10_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.11) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '11_full.nc', NF90_NOWRITE, ncid) )
        elseif (mth.eq.12) then
         call check( nf90_open(waveicedatadir // fname_ww3 // '12_full.nc', NF90_NOWRITE, ncid) )
        endif
        print*, '1st check done'



end program read_ww3
