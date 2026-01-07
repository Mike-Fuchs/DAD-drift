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
  
  !flip y coordinates
  drift_pattern%data = drift_pattern%data(drift_pattern%header%nrows:1:-1, :)

  !check if landscape raster cellsize is unequal to 1
  if (control_dat%cellsize /= 1._i_kind) then
    !drift_pattern = rescale_drift_pattern_fast(drift_pattern, control_dat%cellsize)
  end if
  
  !loop over all field input raster
  do i=1,control_dat%field_count
    !read landscape_input_raster
	landscape_input_raster = read_ascii_grid(in_files%landscape_file_name(i))
	
	!generate temporal landscape_drift_raster
	tmp_raster = add_ascii_grid_landscape(landscape_input_raster, drift_pattern)
	
	!add to output raster
	if(i == 1)then
	  landscape_drift_raster = tmp_raster
	else
	  landscape_drift_raster = add_ascii_grid_fast(landscape_drift_raster,tmp_raster)
	end if
	
	print *, 'finished field', i
  end do
    
  !deallocate arrays
  if(allocated(landscape_input_matrix)) deallocate(landscape_input_matrix)
  if(allocated(landscape_input_matrix_tmp)) deallocate(landscape_input_matrix_tmp)
  if(allocated(tmp_1_mat)) deallocate(tmp_1_mat)
  
  end subroutine landscape_drift