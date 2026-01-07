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

subroutine environment_initialisation
  !! ~~~~ description ~~~~
  !! In this subroutine the environmental parameters are initialised.
    
  use data_module
  use functions_module
  use constants_module
    
  implicit none
  
  ! compute temperature in K
  env_dat%T_K = env_dat%T_C + K_0
  
  ! computerelative hunidity as percent
  env_dat%Rh_P = env_dat%Rh * 100._i_kind
  
  ! compute dynamic viscosity of air
  env_dat%mu_air = sutherland_formula(env_dat%T_K)
  
  ! compute saturation vapor pressure at ambient temperature
  env_dat%Psat = saturated_vapor_pressure(env_dat%T_C)
  
  ! compute air density
  env_dat%rho_air = air_density(env_dat%Psat,env_dat%RH,env_dat%P,env_dat%T_K)
  
  ! compute water density at ambient temperature
  env_dat%rho_h2o = water_density(env_dat%T_C)
  
  !compute absolute humidity at ambient temperature
  env_dat%ah = absolute_humidity(env_dat%rh,env_dat%Psat,env_dat%T_K)
  
  !compute friction velocity
  env_dat%U_fric = (Kc*env_dat%U_wind)/(log(env_dat%z_wind/env_dat%z0))
  
  !! compute wind profile
  !compute d
  env_dat%d = env_dat%Hv * 0.63_i_kind
  
  !calculate wind speed at the top of the canopy
  if(env_dat%Hv > (env_dat%d + 2 * env_dat%z0)) then
    env_dat%Uh = (env_dat%U_fric/Kc)*log((env_dat%Hv-env_dat%d)/env_dat%z0)
  else
    env_dat%Hv = env_dat%d + 2 * env_dat%z0
	env_dat%Uh = (env_dat%U_fric/Kc)*log((env_dat%Hv-env_dat%d)/env_dat%z0)
  end if

  print *, 'environment initialisation successful' 
  
  end subroutine environment_initialisation