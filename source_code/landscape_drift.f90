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

subroutine landscape_drift
  !! ~~~~ description ~~~~
  !! This function will predict the downwind drift deposition based on one or multiple field input raster
  
  use data_module
  use constants_module
  use functions_module
  use ascii_grid_module
  use file_module
  
  implicit none
  real(i_kind), dimension(:,:), allocatable :: landscape_input_matrix				![-]					|landscape input matrix
  real(i_kind), dimension(:,:), allocatable :: landscape_input_matrix_tmp			![-]					|temporal landscape input matrix
  real(i_kind), dimension(:,:), allocatable :: tmp_1_mat							![-]					|temporal matrix
  type(ascii_grid) :: tmp_raster													![-]					|temporal raster
  type(ascii_grid) :: tmp_1															![-]					|temporal raster
  type(ascii_grid) :: tmp_2															![-]					|temporal raster
  real(i_kind) :: res_r																![m]					|cellsize
  integer :: res_i																	![-]					|cellsize rounded
  integer :: i																		![-]					|counter
  
  print *, 'calculating landsape drift pattern'
  
  !check if landscape raster cellsize is unequal to 1
  if(control_dat%cellsize /= 1)then
    !creat ascii grid with suiting cellsize
	res_i = int(control_dat%cellsize)
	res_r = real(control_dat%cellsize,i_kind)
	tmp_1 = create_ascii_grid(res_i, res_i, (-res_r/2._i_kind), (-res_r/2._i_kind), 1._i_kind)
	
	!transforme to matrix
	tmp_1_mat = ascii_grid_to_llz(tmp_1)
	
	!call add_ascii_grid_at
	tmp_2 = add_ascii_grid_at(drift_pattern,tmp_1_mat)
	
	!scale to target resolution
	drift_pattern = rescale_ascii_grid_centered(tmp_2, res_r, 'sum')	
  end if
  
  !loop over all field input raster
  do i=1,control_dat%field_count
    !read landscape_input_raster
	landscape_input_raster = read_ascii_grid(in_files%landscape_file_name(i))
	
	!convert landscape data to matrix
    landscape_input_matrix_tmp = ascii_grid_to_llz(landscape_input_raster)
    
    !select all cells with a value of 1
    landscape_input_matrix = landscape_input_matrix_tmp((pack([(i,i=1,size(landscape_input_matrix_tmp(:,3)))], mask = (abs(landscape_input_matrix_tmp(:,3)-1._i_kind) < 1.e-5_i_kind))),:)
    
    !generate temporal landscape_drift_raster
    tmp_raster = add_ascii_grid_at(drift_pattern,landscape_input_matrix)
	
	!add to output raster
	if(i == 1)then
	  landscape_drift_raster = tmp_raster
	else
	  landscape_drift_raster = add_ascii_grid(landscape_drift_raster,tmp_raster)
	end if
	
	print *, 'finished field', i
  end do
    
  !deallocate arrays
  if(allocated(landscape_input_matrix)) deallocate(landscape_input_matrix)
  if(allocated(landscape_input_matrix_tmp)) deallocate(landscape_input_matrix_tmp)
  if(allocated(tmp_1_mat)) deallocate(tmp_1_mat)
  
  end subroutine landscape_drift