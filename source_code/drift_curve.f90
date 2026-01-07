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

subroutine drift_curve
  !! ~~~~ description ~~~~
  !! This function will predict the drift curve based on defined field dimensions
  
  use data_module
  use constants_module
  use functions_module
  use ascii_grid_module
  
  implicit none
  type(ascii_grid) :: drift_curve_raster										![-]					|drift curve raster
  real(i_kind), dimension(:,:), allocatable :: field_matrix						![-]					|field matrix
  real(i_kind), dimension(:,:), allocatable :: drift_curve_raster_matrix		![-]					|drift curve raster matrix
  real(i_kind), dimension(:), allocatable :: lat								![m]					|vector of lateral coordinates
  real(i_kind), dimension(:), allocatable :: long								![m]					|vector of longitudinal coordinates
  integer :: i																	![-]					|counter
  integer :: width																![m]					|field width
  integer :: length																![m]					|field length
  
  print *, 'calculating drift curve'
  
  !create lat and long arrays
  width = int(app_dat%b_width * real(app_dat%n_pass,i_kind))
  length = int(app_dat%f_length)
  lat = real([(i,i=1,length)],i_kind)
  long = real([(i,i=1,width)],i_kind)
  
  !allocate and create field matrix
  allocate(field_matrix((size(lat)*size(long)),2))
  field_matrix(:,1:2) = expand_grid(lat,long)
  
  !add pattern to field
  drift_curve_raster = add_ascii_grid_at(drift_pattern,field_matrix)
  
  !create xyz matrix of drift raster
  drift_curve_raster_matrix = ascii_grid_to_llz(drift_curve_raster)
  
  !get the middel of the field at lat = 50
  drift_curve_matrix = drift_curve_raster_matrix(pack([(i,i=1,size(drift_curve_raster_matrix(:,1)))], mask = (abs(drift_curve_raster_matrix(:,1) - real(nint((app_dat%f_length/2._i_kind)),i_kind)) < 1e-6_i_kind)),2:3)
  
  !calculate drift as fraction
  drift_curve_matrix(:,2) = drift_curve_matrix(:,2) / (app_dat%app_rate/10000._i_kind)
  
  !flip distance to be positive
  drift_curve_matrix(:,1) = drift_curve_matrix(:,1) * (-1._i_kind)
  
  !shift distance by 0.5m, so 0m is at the border between the last field cell and the first non-taget cell
  drift_curve_matrix(:,1) = drift_curve_matrix(:,1) + 0.5_i_kind
  
  !flip order of element
  drift_curve_matrix = drift_curve_matrix([(i,i=size(drift_curve_matrix(:,1)),1,-1)],:)
  
  !deallocate arrays
  if(allocated(field_matrix)) deallocate(field_matrix)
  if(allocated(drift_curve_raster_matrix)) deallocate(drift_curve_raster_matrix)
  if(allocated(lat)) deallocate(lat)
  if(allocated(long)) deallocate(long)
  
  end subroutine drift_curve