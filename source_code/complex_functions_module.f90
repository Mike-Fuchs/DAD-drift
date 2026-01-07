! Copyright (C) 2024 BASF SE
! 
! This file is part of DAD-drift.
! 
! DAD drift is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; version 2.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; If not, see <https://www.gnu.org/licenses/>.

module complex_functions_module
  !! ~~~~ description ~~~~
  !! This module containes the raster prediction function
  use constants_module
  
  implicit none
  private
  public :: predict_raster
 
  contains
    function predict_raster (i_ds) result(raster_tmp)
	!! ~~~~ description ~~~~
    !! This function will precit the drift raster based on the droplet model run
	
    use data_module
    use constants_module
    use ascii_grid_module
	use functions_module
    
    implicit none
    integer, intent(in) :: i_ds 								![m]					|counter over droplet spectrum
    real(i_kind), dimension(:,:), allocatable :: xyz			![-]					|xyz table
    real(i_kind), dimension(:), allocatable :: lat_coord		![m]					|lateral coordinates
    real(i_kind), dimension(:), allocatable :: long_coord		![m]					|longitudinal coordinates
    real(i_kind), dimension(:), allocatable :: time				![s]					|time vector
    real(i_kind), dimension(:), allocatable :: time1			![s]					|time vector
	real(i_kind), dimension(:), allocatable :: time2			![s]					|time vector
	real(i_kind), dimension(:), allocatable :: xyz_tmp			![-]					|xyz table
	real(i_kind), dimension(5) :: mm							![m]					|min-max vector
    type(ascii_grid) :: raster_tmp								![-]					|temporal raster
    real(i_kind), dimension(3) :: res							![m]					|resolution vector (coordinate, time)
    real(i_kind) :: edv											![m/s]					|effective deposition velocity
    real(i_kind) :: z_off										![m]					|offset in z-direction as function of x as calculated by droplet model
    real(i_kind) :: x_mean										![m]					|mean puff position at deposition height
    real(i_kind) :: dt											![s]					|time delta
	real(i_kind) :: t_0											![s]					|minimal time
	real(i_kind) :: t_0_min										![s]					|minimal time
	real(i_kind) :: mag											![kg]					|mass above ground
	real(i_kind) :: mag_t										![kg]					|mass above ground
	real(i_kind) :: corr										![-]					|correction factor
	real(i_kind) :: sigma_v										![m]					|vertical dispersion parameter
	real(i_kind) :: sigma_h										![m]					|horizontal dispersion parameter
	real(i_kind) :: sigma_v_tmp									![m]					|temporal horizontal dispersion parameter
	real(i_kind) :: fract_AI									![-]					|fraction AI
	integer :: i,j												![-]					|counter
    integer :: t_min											![s]					|minimal time
	integer :: t_max											![s]					|maximal time
	logical :: flag												![-]					|flag
	character(len=20) :: tmp_head								![-]					|temporal header
	
	!allocate
	allocate(raster_tmp%data(0,0))
	
	!reset values
	mm = 0._i_kind
	res = 0._i_kind
	t_0 = 0._i_kind

	!determine min and max for x and y coord
    mm = xyt_profiler(i_ds)
	
	!time vector construction
	if(any(drop_mod(i_ds)%z < control_dat%z_0))then
	  t_0 = interpolate(drop_mod(i_ds)%z,drop_mod(i_ds)%time,control_dat%z_0)
	  if(t_0 < mm(5))then
	    !calculate resolution for first period
		res(2) = time_res2((/mm(4),t_0/))
		!calculate resolution for first period
		res(3) = time_res2((/t_0,mm(5)/))
		
		!construct first time vector
		t_min = floor(mm(4)/res(2))
        dt = (real(t_min,i_kind)*res(2))-mm(4)
	    t_max = floor((t_0+dt)/res(2))
	    time1 = (real((/(i,i=t_min,t_max)/),i_kind)*res(2))-dt 
		
		!construct second time vector
		t_0_min = maxval(time1)+res(2)/2._i_kind+res(3)/2._i_kind
		t_min = floor(t_0_min/res(3))
        t_max = ceiling(mm(5)/res(3))
	    dt = (real(t_min,i_kind)*res(3))-t_0_min
	    time2 = (real((/(i,i=t_min,t_max)/),i_kind)*res(3))-dt 
		
		!combine time vectors
		time = (/time1,time2/)
		flag = .true.
	  else
	    ! 1 step time vector
		!calculate resolution of time
        res(2) = time_res1(mm(4:5))
		
	    !create time vector
        t_min = floor(mm(4)/res(2))
        t_max = ceiling(mm(5)/res(2))
	    dt = (real(t_min,i_kind)*res(2))-mm(4)
	    time = (real((/(i,i=t_min,t_max)/),i_kind)*res(2))-dt
		flag = .false.
	  end if
	else
	  ! 1 step time vector
	  !calculate resolution of time
      res(2) = time_res1(mm(4:5))
	  
	  !create time vector
      t_min = floor(mm(4)/res(2))
      t_max = ceiling(mm(5)/res(2))
	  dt = (real(t_min,i_kind)*res(2))-mm(4)
	  time = (real((/(i,i=t_min,t_max)/),i_kind)*res(2))-dt
	  flag = .false.
	end if
	
	!calculate resolution of raster
    res(1) = xy_res(mm(1:3))
    
	!create xyz table based on wind direction
    call xyz_create(mm(1:3),res(1),xyz,lat_coord,long_coord)
	allocate(xyz_tmp(size(xyz(:,5))))
	
	!set mass above ground
	mag = app_dat%app_rate/10000._i_kind
	mag_t = mag
  	
	!start time loop
    do i=1,size(time)
	  !reset xyz_tmp
	  xyz_tmp = 0._i_kind
	  !interpolatte z_off
  	  z_off = interpolate(drop_mod(i_ds)%time,drop_mod(i_ds)%z,time(i))
  	  !interpolate effective deposition velocity
  	  edv = interpolate(drop_mod(i_ds)%time,drop_mod(i_ds)%edv,time(i))
	  !interpolate fract_AI
	  fract_AI = interpolate(drop_mod(i_ds)%time,drop_mod(i_ds)%fract_AI,time(i))
	  !interpolate x_mean
	  x_mean = interpolate(drop_mod(i_ds)%time,drop_mod(i_ds)%x_mean,time(i))
	  
	  !compute sigma
	  if(z_off > control_dat%z_0)then
	    sigma_h = souton_sigma_h(env_dat%sig_h,time(i),(time(i)*env_dat%U_wind))
	    sigma_v = souton_sigma_v(env_dat%sig_v,time(i))
	    sigma_v_tmp = sigma_v
	  else
	    sigma_h = souton_sigma_h(env_dat%sig_h,time(i),(time(i)*env_dat%U_wind))
	    sigma_v = sigma_v_tmp
	  end if
	  
	  !correction factor
	  corr = max(min(mag/mag_t,1._i_kind),0._i_kind)
	  
	  !start coordinate loop
	  if(flag)then
	    if(time(i) < t_0)then
		  !time smaler t_0
  	      do j=1,size(xyz(:,1))
  	        if (xyz(j,3)>mm(1) .and. xyz(j,3)<=mm(2) .and. xyz(j,4)>=-mm(3) .and. xyz(j,4)<=mm(3)) then
  	          xyz_tmp(j) = (diffusion_model_xyz (xyz(j,3),xyz(j,4),control_dat%z_0,app_dat%app_rate/10000._i_kind,z_off,x_mean,edv,sigma_h,sigma_v)) * res(1)**2._i_kind * res(2)
  	        end if
	      end do
		else
		  !time larger t_0
  	      do j=1,size(xyz(:,1))
  	        if (xyz(j,3)>mm(1) .and. xyz(j,3)<=mm(2) .and. xyz(j,4)>=-mm(3) .and. xyz(j,4)<=mm(3)) then
  	          xyz_tmp(j) = (diffusion_model_xyz (xyz(j,3),xyz(j,4),control_dat%z_0,app_dat%app_rate/10000._i_kind,z_off,x_mean,edv,sigma_h,sigma_v)) * res(1)**2._i_kind * res(3)
  	        end if
	      end do
		end if
	  else
	    !no two time steps
  	    do j=1,size(xyz(:,1))
  	      if (xyz(j,3)>mm(1) .and. xyz(j,3)<=mm(2) .and. xyz(j,4)>=-mm(3) .and. xyz(j,4)<=mm(3)) then
  	        xyz_tmp(j) = (diffusion_model_xyz (xyz(j,3),xyz(j,4),control_dat%z_0,app_dat%app_rate/10000._i_kind,z_off,x_mean,edv,sigma_h,sigma_v)) * res(1)**2._i_kind * res(2)
  	      end if
	    end do
	  end if
	  
	  !update xyz
	  xyz(:,5) = xyz(:,5) + xyz_tmp * corr
	  
	  !update mag
	  mag = mag - sum(xyz_tmp * corr)
	  
	  !compute mass above ground
  	  mag_t = (1._i_kind-cdf_normal_sigma(sigma_v,z_off,control_dat%z_0))*(app_dat%app_rate/10000._i_kind)*fract_AI
	  
    end do
	
    !transform xyz to ascii_grid
    raster_tmp = llz_to_ascii_grid_with_ll(xyz(:,(/1,2,5/)), res(1), lat_coord, long_coord)
	
    !debugging when devmode is active
	if(control_dat%dev_mode)then
	  print *, i_ds, mm, res, size(xyz(:,1)), size(time), minval(time), maxval(time)
	end if
	
	!change resolution of raster to 1 meter  
    if(raster_tmp%header%cellsize /= 1._i_kind) then
    !change resolution
    raster_tmp = rescale_ascii_grid(raster_tmp, 1._i_kind, 'sum')
    end if
	
	!deallocate arrays
	if(allocated(xyz)) deallocate(xyz)
	if(allocated(lat_coord)) deallocate(lat_coord)
	if(allocated(long_coord)) deallocate(long_coord)
	if(allocated(time)) deallocate(time)
	if(allocated(time1)) deallocate(time1)
	if(allocated(time2)) deallocate(time2)
	if(allocated(xyz_tmp)) deallocate(xyz_tmp)
	
	!debugging when devmode is active
    if(control_dat%dev_mode)then
	  write (tmp_head,128) i_ds
	  call write_ascii_grid(raster_tmp, 'output/ds_'//trim(tmp_head)//'.asc')
    end if

	128 format (1I3.3)	
	
    end function predict_raster

end module complex_functions_module