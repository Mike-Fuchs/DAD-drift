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

subroutine droplet_model_allo (n, n_dc, n_x)
  !! ~~~~ description ~~~~
  !! This function will allocate the droplet and environment data frames.
  
  use data_module
  
  implicit none
  integer, intent(in) :: n								![-]					|number of time steps
  integer, intent(in) :: n_dc							![-]					|number of droplet classes
  integer, intent(in) :: n_x							![-]					|length of sd vector
  integer :: i											![-]					|counter
  
  !allocate local environment
  allocate(loc_env%time(n))
  allocate(loc_env%volume(n))
  allocate(loc_env%d_volume(n))
  allocate(loc_env%m_air(n))
  allocate(loc_env%d_m_air(n))
  allocate(loc_env%m_h2o(n))
  allocate(loc_env%d_m_h2o(n))
  allocate(loc_env%d_E_air(n))
  allocate(loc_env%d_E_h2o(n))
  allocate(loc_env%E_tot_t1(n))
  allocate(loc_env%m_h2o_drop(n))
  allocate(loc_env%d_m_h2o_evap(n))
  allocate(loc_env%d_E_h2o_evap(n))
  allocate(loc_env%T_sys_K(n))
  allocate(loc_env%T_sys_C(n))
  allocate(loc_env%Psat(n))
  allocate(loc_env%Rh(n))
  allocate(loc_env%mu_air(n))
  allocate(loc_env%rho_air(n))
  allocate(loc_env%T_wb_C(n))
  allocate(loc_env%T_f_C(n))
  allocate(loc_env%Psat_f(n))
  allocate(loc_env%mu_air_f(n))
  allocate(loc_env%rho_air_f(n))
  allocate(loc_env%rho_h2o_f(n))
  allocate(loc_env%D_h2o_f(n))
  allocate(loc_env%Sc_h2o_f(n))
  allocate(loc_env%dt(n))
  
  !allocate number of droplet classes to drop_mod
  allocate(drop_mod(n_dc))
  do i=1, n_dc
    allocate(drop_mod(i)%time(n))
	allocate(drop_mod(i)%diameter(n))
	allocate(drop_mod(i)%volume(n))
	allocate(drop_mod(i)%m_h2o(n))
	allocate(drop_mod(i)%m_AI(n))
	allocate(drop_mod(i)%fract_AI(n))
	allocate(drop_mod(i)%x_h2o(n))
	allocate(drop_mod(i)%x(n))
	allocate(drop_mod(i)%z(n))
	allocate(drop_mod(i)%V_x(n))
	allocate(drop_mod(i)%V_z(n))
	allocate(drop_mod(i)%Re_s(n))
	allocate(drop_mod(i)%f(n))
	allocate(drop_mod(i)%tau_m(n))
	allocate(drop_mod(i)%Re_f(n))
	allocate(drop_mod(i)%Sh_h2o_f(n))
	allocate(drop_mod(i)%F_h2o(n))
	allocate(drop_mod(i)%F_AI(n))
	allocate(drop_mod(i)%dV_x_dt(n))
	allocate(drop_mod(i)%dV_z_dt(n))
	allocate(drop_mod(i)%dm_h2o_dt(n))
	allocate(drop_mod(i)%dm_AI_dt(n))
	allocate(drop_mod(i)%dm_h2o(n))
	allocate(drop_mod(i)%dm_AI(n))
	allocate(drop_mod(i)%AI_a(n))
	allocate(drop_mod(i)%dm_h2o_a(n))
	allocate(drop_mod(i)%V_s(n))
	allocate(drop_mod(i)%edv(n))
	allocate(drop_mod(i)%dag(n))
	allocate(drop_mod(i)%dag_res(n))
	allocate(drop_mod(i)%dag_res_chg(n))
	allocate(drop_mod(i)%dt(n))
	allocate(drop_mod(i)%Uz(n))
	allocate(drop_mod(i)%Ud(n))
	allocate(drop_mod(i)%x_mean(n))
	allocate(drop_mod(i)%sigma_h(n))
	allocate(drop_mod(i)%sigma_v(n))
	allocate(drop_mod(i)%x_s(n,n_x))
  end do
  
  end subroutine droplet_model_allo