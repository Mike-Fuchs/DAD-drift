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

module header_module
  !! ~~~~ description ~~~~
  !! In this module various header for result printing are defined.
    
  implicit none
  
  type :: header_droplet_matrix
    character(len=4) :: a = 'time'
	character(len=8) :: b = 'diameter'
	character(len=6) :: c = 'volume'
	character(len=5) :: d = 'm_h2o'
	character(len=4) :: e = 'm_AI'
	character(len=13):: f = 'm_AI/m_AI(t0)'
	character(len=7) :: g = 'x_water'
	character(len=1) :: h = 'x'
	character(len=1) :: i = 'z'
	character(len=3) :: j = 'V_x'
	character(len=3) :: k = 'V_z'
	character(len=15):: l = 'reynolds_number'
	character(len=1) :: m = 'f'
	character(len=5) :: n = 'tau_m'
	character(len=6) :: o = 'Sh_h2o'
	character(len=5) :: p = 'Sh_AI'
	character(len=5) :: q = 'F_h2o'
	character(len=4) :: r = 'F_AI'
	character(len=7) :: s = 'dV_x/dt'
	character(len=7) :: t = 'dV_z/dt'
	character(len=9) :: u = 'dm_h2o/dt'
	character(len=8) :: v = 'dm_AI/dt'
	character(len=2) :: w = 'dt'
  end type header_droplet_matrix
  type(header_droplet_matrix) :: head_dm

  type :: units_droplet_matrix
    character(len=3) :: a = '[s]'
	character(len=3) :: b = '[m]'
	character(len=4) :: c = '[m³]'
	character(len=4) :: d = '[kg]'
	character(len=4) :: e = '[kg]'
	character(len=3) :: f = '[-]'
	character(len=3) :: g = '[-]'
	character(len=3) :: h = '[m]'
	character(len=3) :: i = '[m]'
	character(len=5) :: j = '[m/s]'
	character(len=5) :: k = '[m/s]'
	character(len=3) :: l = '[-]'
	character(len=3) :: m = '[-]'
	character(len=3) :: n = '[s]'
	character(len=3) :: o = '[-]'
	character(len=3) :: p = '[-]'
	character(len=8) :: q = '[kg/m²s]'
	character(len=8) :: r = '[kg/m²s]'
	character(len=6) :: s = '[m/s²]'
	character(len=6) :: t = '[m/s²]'
	character(len=6) :: u = '[kg/s]'
	character(len=6) :: v = '[kg/s]'
	character(len=3) :: w = '[s]'
  end type units_droplet_matrix
  type(units_droplet_matrix) :: units_dm
  
  type :: header_prediction_matrix
    character(len=3) :: a = 'lat'
	character(len=4) :: b = 'long'
	character(len=1) :: f = 'x'
	character(len=1) :: g = 'y'	
	character(len=13) :: c = 'concentration'
	character(len=9) :: d = 'mass_flux'
	character(len=5) :: e = 'drift'	
  end type header_prediction_matrix
  type(header_prediction_matrix) :: head_pm
  
  type :: units_prediction_matrix
    character(len=3) :: a = '[m]'
	character(len=3) :: b = '[m]'
	character(len=3) :: f = '[m]'
	character(len=3) :: g = '[m]'
	character(len=7) :: c = '[kg/m³]'
	character(len=8) :: d = '[kg/m²s]'
	character(len=3) :: e = '[-]'	
  end type units_prediction_matrix
  type(units_prediction_matrix) :: units_pm
  
  type :: header_drift_curve
    character(len=8) :: a = 'distance'
	character(len=5) :: b = 'drift'
  end type header_drift_curve
  type(header_drift_curve) :: head_dc
  
  type :: units_drift_curve
    character(len=3) :: a = '[m]'
	character(len=3) :: b = '[-]'
  end type units_drift_curve
  type(units_drift_curve) :: units_dc
  
  end module header_module