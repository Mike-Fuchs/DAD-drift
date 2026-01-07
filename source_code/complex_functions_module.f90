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
    real(i_kind), dimension(:), allocatable :: xyz_tmp			![-]					|xyz table
    type(ascii_grid) :: raster_tmp								![-]					|temporal raster
    real(i_kind), dimension(2) :: res							![m]					|resolution vector (coordinate, time)
    real(i_kind) :: V_z											![m/s]					|effective deposition velocity
    real(i_kind) :: z_mean										![m]					|offset in z-direction as function of x as calculated by droplet model
    real(i_kind) :: x_mean										![m]					|mean puff position at deposition height
    real(i_kind) :: dt											![s]					|time delta
	real(i_kind) :: sigma_v										![m]					|vertical dispersion parameter
	real(i_kind) :: sigma_h										![m]					|horizontal dispersion parameter
	real(i_kind) :: fract_AI									![-]					|fraction AI
	integer :: i,j												![-]					|counter
    integer :: t_min											![s]					|minimal time
	integer :: t_max											![s]					|maximal time
	integer :: t_idx_1											![]						|
	integer :: t_idx_2											![]						|
	character(len=20) :: tmp_head								![-]					|temporal header
	
	!allocate
	allocate(raster_tmp%data(0,0))
	
	!reset values
	res = 0._i_kind
	
	!determine time resolution
	res(2) = time_res(drop_mod(i_ds)%puff_dim%t1,drop_mod(i_ds)%puff_dim%t2)
	
	!construct time vector
	t_min = floor(drop_mod(i_ds)%puff_dim%t1 / res(2))
    t_max = floor(drop_mod(i_ds)%puff_dim%t2 / res(2))
    dt = (real(t_min, i_kind) * res(2)) - drop_mod(i_ds)%puff_dim%t1
    time = (real((/(i, i = t_min, t_max)/), i_kind) * res(2)) - dt
	
	!calculate resolution of raster
    res(1) = min(xy_res(drop_mod(i_ds)%puff_dim%x_min, drop_mod(i_ds)%puff_dim%x_max, drop_mod(i_ds)%puff_dim%y_max),8._i_kind)
    
    !create xyz table based on wind direction
    call xyz_create(drop_mod(i_ds)%puff_dim%x_min, drop_mod(i_ds)%puff_dim%x_max, drop_mod(i_ds)%puff_dim%y_max, res(1), xyz, lat_coord, long_coord)
    
	!allocate xyz
    allocate(xyz_tmp(size(xyz(:,5))))
	
	!start time loop
    do i=1,size(time)
	  !reset xyz_tmp
	  xyz_tmp = 0._i_kind
	  
	  !define time index t_idx_1 and t_idx_2
	  t_idx_1 = maxloc(drop_mod(i_ds)%time, dim=1, mask=drop_mod(i_ds)%time <= time(i))
      t_idx_2 = t_idx_1 + 1
	  
	  !interpolatte z_mean
  	  z_mean = interpolate_idx(drop_mod(i_ds)%time,drop_mod(i_ds)%z,t_idx_1,t_idx_2,time(i))
  	  !interpolate effective deposition velocity
  	  V_z = interpolate_idx(drop_mod(i_ds)%time,drop_mod(i_ds)%V_z,t_idx_1,t_idx_2,time(i))
	  !interpolate fract_AI
	  fract_AI = interpolate_idx(drop_mod(i_ds)%time,drop_mod(i_ds)%fract_AI,t_idx_1,t_idx_2,time(i))
	  
	  !interpolate x_mean
	  x_mean = interpolate_2d_geo(drop_mod(i_ds)%time, drop_mod(i_ds)%puff_slc%z, drop_mod(i_ds)%puff_slc%x_trav, t_idx_1, t_idx_2, drop_mod(i_ds)%z_idx(t_idx_1,1), drop_mod(i_ds)%z_idx(t_idx_1,2), drop_mod(i_ds)%z_idx(t_idx_2,1), drop_mod(i_ds)%z_idx(t_idx_2,2), time(i), control_dat%z_0)
	  !interpolate sigma_h
	  sigma_h = interpolate_2d_geo(drop_mod(i_ds)%time, drop_mod(i_ds)%puff_slc%z, drop_mod(i_ds)%puff_slc%sigma_h_out, t_idx_1, t_idx_2, drop_mod(i_ds)%z_idx(t_idx_1,1), drop_mod(i_ds)%z_idx(t_idx_1,2), drop_mod(i_ds)%z_idx(t_idx_2,1), drop_mod(i_ds)%z_idx(t_idx_2,2), time(i), control_dat%z_0)
	  
	  !start coordinate loop
	  do j=1,size(xyz(:,1))
  	    if (xyz(j,3)>=drop_mod(i_ds)%puff_dim%x_min .and. xyz(j,3)<=drop_mod(i_ds)%puff_dim%x_max .and. xyz(j,4)>=-drop_mod(i_ds)%puff_dim%y_max .and. xyz(j,4)<=drop_mod(i_ds)%puff_dim%y_max) then
  	      xyz_tmp(j) = (diffusion_model_xyz_integrated(xyz(j,3),xyz(j,4),control_dat%z_0,app_dat%app_rate/10000._i_kind*(1._i_kind/0.997300203936740_i_kind)*fract_AI,z_mean,x_mean,V_z,sigma_h,time(i),res(1))) * res(2)
  	    end if
	  end do
	    
	  !update xyz
	  xyz(:,5) = xyz(:,5) + xyz_tmp
    end do
	
	!transform xyz to ascii_grid
    raster_tmp = llz_to_ascii_grid_with_ll(xyz(:,(/1,2,5/)), res(1), lat_coord, long_coord)
	
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
	if(allocated(xyz_tmp)) deallocate(xyz_tmp)
	
	128 format (1I3.3)	
	
    end function predict_raster

end module complex_functions_module