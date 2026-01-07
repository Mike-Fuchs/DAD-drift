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

subroutine predict_drift_pattern
  !! ~~~~ description ~~~~
  !! This function will predict the downwind drift deposition pattern based on a single sqaure meter as emitter.
  
  use data_module
  use ascii_grid_module
  use complex_functions_module
  
  implicit none
  type(ascii_grid) :: raster_tmp1										![-]					|temporal raster
  integer :: i_ds														![-]					|counter
  integer :: n_dc														![-]					|counter
  integer :: i															![-]					|counter
  
  !loop over droplet spectrum
  n_dc = size(drop_spec%ds)
  do i_ds = 1, n_dc
	grid_list(i_ds) = predict_raster(i_ds)
	print '(1A20,1F5.1,1A11)', ' pattern prediction ', (real(i_ds)/real(n_dc))*100. , '% completed'
  end do
  
  !define raster as first raster in grid_list with fraction larger zero
  drift_pattern = grid_list(1)
  drift_pattern%data = drift_pattern%data * drop_spec%f(1)
  
  !loop over the remaining droplet sizes
  do i = 2, size(drop_spec%f)
    !define tmp2
    raster_tmp1 = grid_list(i)
	raster_tmp1%data = raster_tmp1%data * drop_spec%f(i)
	
	!add raster
	drift_pattern = add_ascii_grid(drift_pattern, raster_tmp1)
  end do
  
  print *, 'drift pattern calculated'
  
  end subroutine predict_drift_pattern