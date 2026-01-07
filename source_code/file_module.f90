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

module file_module
  !! ~~~~ description ~~~~
  !! In this module the input and output file names are defined.
    
  implicit none
  
  type input_files
    character(len=40) :: env_input = "input/environment_input.txt"				!unit = 100
	character(len=40) :: spec_input												!unit = 101
	character(len=40) :: app_input = "input/application_input.txt"				!unit = 102
	character(len=40) :: gaus_input = "input/gaussian_input.txt"				!unit = 103
	character(len=40) :: controle_input = "input/control_input.txt"				!unit = 104
	character(len=80), dimension(:), allocatable :: landscape_file_name			!unit = raster
  end type input_files
  type(input_files) :: in_files
  
  
  type output_files
    character(len=40) :: env_output = "output/environmental_data.txt"			!unit = 110
	character(len=40) :: drop_output = "output/droplet_output.txt"				!unit = 111
	character(len=40) :: predict_output = "output/prediction_output.txt"		!unit = 112
	character(len=40) :: drift_curve_output = "output/drift_curve_output.txt"	!unit = 113
	character(len=40) :: drift_pattern = "output/drift_pattern.asc"				!unit = raster
	character(len=40) :: landscape_drift = "output/landscape_drift.asc"			!unit = raster
  end type output_files
  type(output_files) :: out_files
  
  end module file_module