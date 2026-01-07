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

subroutine application_initialisation
  !! ~~~~ description ~~~~
  !! In this subroutine the application parameters are initialised.
    
  use data_module
  use constants_module
    
  implicit none
  real(i_kind) :: tpa									![h/ha]					|time per area
  real(i_kind) :: m_AI									![kg]					|mass of active ingredient
  real(i_kind) :: m_H2O									![kg]					|mass of water
  real(i_kind) :: V_AI									![m³]					|volume of active ingredient
  real(i_kind) :: V_H2O									![m³]					|volume of water
  
  !calculate density of spray solution
  m_AI = app_dat%c_solution * 1._i_kind
  V_AI = m_AI / app_dat%rho_AI
  V_H2O = 1._i_kind - V_AI
  m_H2O = V_H2O * env_dat%rho_h2o
  app_dat%rho_solution = (m_AI+m_H2O)/(V_AI+V_H2O)
  
  !calculate initial droplet sheet velocity
  app_dat%V_initial = (2._i_kind*app_dat%app_pressure/app_dat%rho_solution)**0.5_i_kind
  
  !calculate horizontal droplet velocity
  app_dat%V_horizontal = cosd(app_dat%nozzle_angle)*app_dat%V_initial
  
  !calculate application rate [kg/ha] if not defined in input
  if (abs(app_dat%app_rate +99._i_kind) < 1e-6_i_kind) then
    !if volume per ha is not defined
	if (abs(app_dat%app_rate_mha +99._i_kind) < 1e-6_i_kind) then
	  !calculate time per area
	  tpa = (1._i_kind / (app_dat%v_trac * app_dat%b_width)) / 3600._i_kind * 10000._i_kind
	  
	  !calculate volume per area
	  app_dat%app_rate_mha = app_dat%app_rate_mh * tpa
	end if
	
	!application rate
	app_dat%app_rate = app_dat%app_rate_mha * app_dat%c_solution
  end if

  !calculate application rate [m³/ha] if not defined in input or already calculated
  if (abs(app_dat%app_rate_mha +99._i_kind) < 1e-6_i_kind) then
    app_dat%app_rate_mha = app_dat%app_rate / app_dat%c_solution
  end if
   
  print *, 'application initialisation successful'
  
  end subroutine application_initialisation