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

subroutine droplet_model
  !! ~~~~ description ~~~~
  !! This function will run the droplet and the climate model in parallel
  
  use data_module
  use constants_module
  use functions_module
  
  implicit none
  
  real(i_kind), dimension(:), allocatable :: dt_min, dt_vel, dt_mass, dt_dist	![s]					|
  real(i_kind) :: dt_tmp														![s]					|temporal time step width
  real(i_kind), dimension(:), allocatable :: vol_fact							![-]					|vector of volume factors for each droplet class
  real(i_kind), dimension(:), allocatable :: min_x
  real(i_kind) :: xmol															![-]					|molar fraction
  real(i_kind) :: ratio															![-]					|
  real(i_kind) :: mh2o															![-]					|mass of water
  real(i_kind) :: evap_fact														![-]					|evaporation fraction
  real(i_kind) :: fract_limit
  integer, dimension(:), allocatable :: sel_vec									![-]					|vector for data selection
  integer, dimension(:), allocatable :: iter_vec								![-]					|vector of maximal iteration steps
  integer :: s_max																![-]					|maximal step count
  integer :: n_dc																![-]					|number of droplet classes
  integer :: t																	![-]					|counter looping over time
  integer :: d																	![-]					|counter looping over droplet classes
  integer :: i																	![-]					|counter
  integer :: s																	![-]					|counter
  integer, dimension(:), allocatable :: exit_code								![-]					|counter
  logical, dimension(:), allocatable :: ac										![-]					|logical vector if class is activly simulated, droplets z-position is above a given threshold
  logical, dimension(:), allocatable :: vv										![-]					|logical vector if velocity is deviating from ideal speed
  
  !fraction change limiter
  fract_limit = 100._i_kind
  
  !set s_max
  s_max = 15000
  
  !set evap_fact
  evap_fact = 1.0_i_kind
  
  !number of droplet classes
  n_dc = size(drop_spec%ds)
  
  !!total volume correction factor
  allocate(vol_fact(n_dc))
  do d=1, n_dc
	!calculate application rate in m³/m² and calculate how many droplets of a given size are needed to fill that volume
	vol_fact(d) = (drop_spec%f(d)*app_dat%app_rate_mha/10000._i_kind)/((1._i_kind/6._i_kind)*pi*((drop_spec%ds(d)**3._i_kind)))	
  end do
  
  !allocate active classes vector
  allocate(ac(n_dc))
  ac = .true.
  
  !allocate velocity vector
  allocate(vv(n_dc))
  vv = .true.
  
  !allocate exit code
  allocate(exit_code(n_dc))
  exit_code = 0
  
  ! allocate min distance vector
  allocate(min_x(n_dc))
  
  !allocate delta time vectors
  allocate(dt_min(n_dc))
  allocate(dt_vel(n_dc))
  allocate(dt_mass(n_dc))
  allocate(dt_dist(n_dc))
  
  !allocate interation count vector
  allocate(iter_vec(n_dc))
  iter_vec = s_max
  
  !allocate droplet model and loval environment
  call droplet_model_allo(s_max, n_dc)
  
  !!initialise environment
  !time
  loc_env%time(1) = 0._i_kind
  !volume
  loc_env%volume(1) = 1._i_kind * 1._i_kind * (app_dat%b_height + souton_sigma_v(env_dat%sig_v,loc_env%time(1)))
  !mass water
  loc_env%m_h2o(1) = env_dat%ah * loc_env%volume(1)
  !mass air
  loc_env%m_air(1) = loc_env%volume(1)*env_dat%rho_air - loc_env%m_h2o(1)
  !total internal energy
  loc_env%E_tot_t1(1) = (loc_env%m_h2o(1)*env_dat%T_K*Cp_h2o) + (loc_env%m_air(1)*env_dat%T_K*Cp_air)
  !systemn temperature in Kelvin
  loc_env%T_sys_K(1) = env_dat%T_K
  !system temperature in Celcius
  loc_env%T_sys_C(1) = env_dat%T_C
  !saturation vapour pressure
  loc_env%Psat(1) = saturated_vapor_pressure(loc_env%T_sys_C(1))
  !relative humidty
  loc_env%Rh(1) = relative_humidity(loc_env%Psat(1),loc_env%m_h2o(1),loc_env%T_sys_K(1),loc_env%volume(1))
  !dynamic viscosity of air
  loc_env%mu_air(1) = sutherland_formula(loc_env%T_sys_K(1))
  !air density
  loc_env%rho_air(1) = air_density(loc_env%Psat(1),loc_env%Rh(1),env_dat%P,loc_env%T_sys_K(1))
  !wet-bulb temperature
  loc_env%T_wb_C(1) = wet_bulb_temperature(loc_env%T_sys_C(1), loc_env%Rh(1)*100._i_kind)
  !film temperature
  loc_env%T_f_C(1) = (loc_env%T_wb_C(1)+loc_env%T_sys_C(1))/2._i_kind
  !saturation vapour pressure at film temperature
  loc_env%Psat_f(1) = saturated_vapor_pressure(loc_env%T_f_C(1))
  !dynamic viscosity of air at film temperature
  loc_env%mu_air_f(1) = sutherland_formula(loc_env%T_f_C(1)+K_0)
  !air density at film temperature
  loc_env%rho_air_f(1) = air_density(loc_env%Psat_f(1),loc_env%Rh(1),env_dat%P,loc_env%T_f_C(1)+K_0)
  !water density at film temperature
  loc_env%rho_h2o_f(1) = water_density(loc_env%T_f_C(1))
  !diffusion coefficient of water at film temperature
  loc_env%D_h2o_f(1) = diffusion_coefficient_water(loc_env%T_f_C(1))
  !Schmidt number at film temperature
  loc_env%Sc_h2o_f(1) = loc_env%mu_air_f(1)/(loc_env%rho_air_f(1)*loc_env%D_h2o_f(1))
   
  
  !!setting unused variables to zero
  loc_env%d_volume(1) = 0._i_kind
  loc_env%d_m_air(1) = 0._i_kind
  loc_env%d_m_h2o(1) = 0._i_kind
  loc_env%d_E_air(1) = 0._i_kind
  loc_env%d_E_h2o(1) = 0._i_kind
  loc_env%m_h2o_drop = 0._i_kind
  loc_env%d_m_h2o_evap = 0._i_kind
  
  !!initialising droplet classes
  do d=1, n_dc
	!time
	drop_mod(d)%time(1) = 0._i_kind
	!diameter
	drop_mod(d)%diameter(1) = drop_spec%ds(d)
	!volume
	drop_mod(d)%volume(1) = (1._i_kind/6._i_kind)*pi*((drop_mod(d)%diameter(1)**3._i_kind))
	!mass AI
	drop_mod(d)%m_AI(1) = app_dat%c_solution*drop_mod(d)%volume(1)
	!mass h2o
	drop_mod(d)%m_h2o(1) = env_dat%rho_h2o*(drop_mod(d)%volume(1)-(drop_mod(d)%m_AI(1)/app_dat%rho_AI))
	!fraction of AI
	drop_mod(d)%fract_AI(1) = 1._i_kind
	!molar fraction of water
	drop_mod(d)%x_h2o(1) = (drop_mod(d)%m_h2o(1)/mm_h2o)/((drop_mod(d)%m_h2o(1)/mm_h2o)+(drop_mod(d)%m_AI(1)/app_dat%mm_AI))
	!x-position
	drop_mod(d)%x(1) = 0._i_kind
	!z-position
	drop_mod(d)%z(1) = app_dat%b_height
	!x-velocity
	drop_mod(d)%V_x(1) = app_dat%V_vertical
	!z-velocity
	drop_mod(d)%V_z(1) = app_dat%V_horizontal
	!Reynolds number at system temperature
	drop_mod(d)%Re_s(1) = (drop_mod(d)%diameter(1)*loc_env%rho_air(1)*(((env_dat%U_wind-drop_mod(d)%V_x(1))**2._i_kind+(drop_mod(d)%V_z(1))**2._i_kind)**0.5_i_kind))/(loc_env%mu_air(1))
	!Boothroyd's f
	drop_mod(d)%f(1) = 1.0_i_kind+0.15_i_kind*drop_mod(d)%Re_s(1)**0.687_i_kind
	!tau_m
	drop_mod(d)%tau_m(1) = (((drop_mod(d)%m_AI(1)+drop_mod(d)%m_h2o(1))/drop_mod(d)%volume(1))*drop_mod(d)%diameter(1)**2._i_kind)/(18._i_kind*loc_env%mu_air(1)*drop_mod(d)%f(1))
	!Reynolds number at film temperature
	drop_mod(d)%Re_f(1) = (drop_mod(d)%diameter(1)*loc_env%rho_air_f(1)*(((env_dat%U_wind-drop_mod(d)%V_x(1))**2._i_kind+(drop_mod(d)%V_z(1))**2._i_kind)**0.5_i_kind))/(loc_env%mu_air_f(1))
	!Sherwood number h2o
	drop_mod(d)%Sh_h2o_f(1) = 1.0_i_kind+0.276_i_kind*drop_mod(d)%Re_f(1)**0.5_i_kind*loc_env%Sc_h2o_f(1)**(1._i_kind/3._i_kind)
	!flux of water
	drop_mod(d)%F_h2o(1) = ((drop_mod(d)%Sh_h2o_f(1)*loc_env%D_h2o_f(1)*mm_h2o)/(drop_mod(d)%diameter(1)*R*(loc_env%T_f_C(1)+K_0)))*((loc_env%Psat(1)*loc_env%Rh(1))-(loc_env%Psat_f(1)*drop_mod(d)%x_h2o(1))) * evap_fact
	!flux of AI
	drop_mod(d)%F_AI(1) = -4.07E-7_i_kind*app_dat%P_AI*app_dat%mm_AI*(1._i_kind-drop_mod(d)%x_h2o(1))
	!rate of change of x-velocity
	drop_mod(d)%dV_x_dt(1) = (env_dat%U_wind-drop_mod(d)%V_x(1))/drop_mod(d)%tau_m(1)
	!rate of change of z-velocity
	drop_mod(d)%dV_z_dt(1) = g -(drop_mod(d)%V_z(1)/drop_mod(d)%tau_m(1))
	!rate of change of water mass
	drop_mod(d)%dm_h2o_dt(1) = pi*drop_mod(d)%diameter(1)**2._i_kind*drop_mod(d)%F_h2o(1)
	!rate of change of AI mass
	drop_mod(d)%dm_AI_dt(1) = pi*drop_mod(d)%diameter(1)**2._i_kind*drop_mod(d)%F_AI(1)
	!settling velocity
	!drop_mod(d)%V_s(1) = g*drop_mod(d)%tau_m(1)
	drop_mod(d)%V_s(1) = settling_velocity (drop_mod(d)%diameter(1),drop_mod(d)%V_z(1),loc_env%rho_air(1),loc_env%mu_air(1),((drop_mod(d)%m_AI(1)+drop_mod(d)%m_h2o(1))/drop_mod(d)%volume(1)))
	!effective deposition velocity
	drop_mod(d)%edv(1) = (drop_mod(d)%V_z(1))*drop_mod(d)%fract_AI(1)
	! sigma calculation
	drop_mod(d)%sigma_v(1) = 0._i_kind
	!calculate distribution above ground
    drop_mod(d)%dag(1) = 1._i_kind
    !calculate distribution above ground for resolution
    drop_mod(d)%dag_res(1) = 1._i_kind
    !
    drop_mod(d)%dag_res_chg(1) = 0._i_kind
	
	! puff slices
	do s = 1, N_slices
      ! vertical position of slice
      drop_mod(d)%puff_slc%z(1,s) = drop_mod(d)%z(1) + (slice_multi(s) * drop_mod(d)%sigma_v(1))
  
      ! wind speed at slice height
      drop_mod(d)%puff_slc%Uz(1,s) = wind_profile(drop_mod(d)%puff_slc%z(1,s))
  
      ! initialize travelled distance (none yet)
      drop_mod(d)%puff_slc%x_trav(1,s) = 0._i_kind
  
      ! horizontal dispersion at t=0 (zero or small epsilon)
      drop_mod(d)%puff_slc%sigma_h_out(1,s) = 0._i_kind
    end do
	
	! initialize puff dimensions
    drop_mod(d)%puff_dim%t1     = 0._i_kind
    drop_mod(d)%puff_dim%t2     = 0._i_kind
    drop_mod(d)%puff_dim%x_min  = 0._i_kind
    drop_mod(d)%puff_dim%x_max  = 0._i_kind
    drop_mod(d)%puff_dim%y_max  = 0._i_kind
	
	!determine time step width
	if(vv(d)) then
	  !calculate dt's
	  dt_vel(d) = abs((drop_mod(d)%V_z(1)/abs(drop_mod(d)%dV_z_dt(1))))/fract_limit
	  dt_mass(d) = abs((drop_mod(d)%m_h2o(1) / drop_mod(d)%dm_h2o_dt(1)))/fract_limit
	  if(drop_mod(d)%V_x(1) < 1e-3)then
		    dt_dist(d) = 1._i_kind
		  else
		    dt_dist(d) = control_dat%max_dist_step / drop_mod(d)%V_x(1)
		  end if
	  !check if velocity is still different from direct solution
	  if(abs((drop_mod(d)%V_z(1)-g*drop_mod(d)%tau_m(1))/drop_mod(d)%V_z(1)) < 0.001_i_kind)then
	    vv(d) = .false.
	  end if
	  
	  !set minimal dt
	  dt_min(d) = max(minval((/dt_vel(d),dt_mass(d),dt_dist(d)/)),1e-6)
	else
	  !calculate dt_mass
	  dt_mass(d) = abs(drop_mod(d)%m_h2o(1) / drop_mod(d)%dm_h2o_dt(1))/fract_limit
	  
	  !set minimal dt
	  dt_min(d) = max(minval((/dt_mass(d),dt_dist(d)/)),1e-6)
	end if
  end do
  
  !set loc time step width
  loc_env%dt(1) = minval(dt_min)
	
  !!calculate evaporating masses
  do d=1, n_dc
    !change in mass of h2o
    drop_mod(d)%dm_h2o(1) = drop_mod(d)%dm_h2o_dt(1) * loc_env%dt(1)
    !change in mass of AI
    drop_mod(d)%dm_AI(1) = drop_mod(d)%dm_AI_dt(1) * loc_env%dt(1)
    !calculate AI loss to atmosphere
	drop_mod(d)%AI_a(1) = abs(drop_mod(d)%dm_AI(1))
	!calculate h2o mass added to ambient system
    drop_mod(d)%dm_h2o_a(1) = -drop_mod(d)%dm_h2o(1)*vol_fact(d)*drop_mod(d)%dag(1)
    !summ mass
    loc_env%m_h2o_drop(1) = loc_env%m_h2o_drop(1) + drop_mod(d)%m_h2o(1)*vol_fact(d)*drop_mod(d)%dag(1)
    !summ evaporated mass
    loc_env%d_m_h2o_evap(1) = loc_env%d_m_h2o_evap(1) + drop_mod(d)%dm_h2o_a(1)
  end do
  
  !calculate energy losst due to evaporation
  loc_env%d_E_h2o_evap(1) = -(loc_env%d_m_h2o_evap(1)*latent_heat_of_vaporization(loc_env%T_sys_K(1))) + ((loc_env%T_wb_C(1)+K_0)*loc_env%d_m_h2o_evap(1)*Cp_h2o)
    
  !!time loop
  t = 2
  do while (any(ac))
    !!update environment
    !time
    loc_env%time(t) = loc_env%time(t-1) + loc_env%dt(t-1)
    !volume
    loc_env%volume(t) = 1._i_kind * 1._i_kind * (app_dat%b_height + souton_sigma_v(env_dat%sig_v,loc_env%time(t)))
    !change in volume
    loc_env%d_volume(t) = loc_env%volume(t) - loc_env%volume(t-1)
    !change of water mass
    loc_env%d_m_h2o(t) = env_dat%ah * loc_env%d_volume(t)
    !change of air mass
    loc_env%d_m_air(t) = loc_env%d_volume(t)*env_dat%rho_air - loc_env%d_m_h2o(t)
    !energy contained in added air
    loc_env%d_E_h2o(t) = loc_env%d_m_h2o(t)*env_dat%T_K*Cp_h2o
    !energy contained in added water
    loc_env%d_E_air(t) = loc_env%d_m_air(t)*env_dat%T_K*Cp_air
    !update total internal energy
    loc_env%E_tot_t1(t) = loc_env%E_tot_t1(t-1)+loc_env%d_E_h2o(t)+loc_env%d_E_air(t)+loc_env%d_E_h2o_evap(t-1)
    if(t>2)then
      loc_env%E_tot_t1(t) = loc_env%E_tot_t1(t) + ((loc_env%T_wb_C(t-2)-loc_env%T_wb_C(t-1))*loc_env%m_h2o_drop(t-1)*Cp_h2o)
    end if
    !update mass of water
    loc_env%m_h2o(t) = loc_env%m_h2o(t-1)+loc_env%d_m_h2o(t)+loc_env%d_m_h2o_evap(t-1)
    !update mass of air
    loc_env%m_air(t) = loc_env%m_air(t-1)+loc_env%d_m_air(t)
    !systemn temperature in Kelvin
    loc_env%T_sys_K(t) = loc_env%E_tot_t1(t)/((Cp_air*loc_env%m_air(t))+(Cp_h2o*loc_env%m_h2o(t)))
    !system temperature in Celcius
    loc_env%T_sys_C(t) = loc_env%T_sys_K(t)-K_0
    !saturation vapour pressure
    loc_env%Psat(t) = saturated_vapor_pressure(loc_env%T_sys_C(t))
    !relative humidty
    loc_env%Rh(t) = relative_humidity(loc_env%Psat(t),loc_env%m_h2o(t),loc_env%T_sys_K(t),loc_env%volume(t))
    !dynamic viscosity of air
    loc_env%mu_air(t) = sutherland_formula(loc_env%T_sys_K(t))
    !air density
    loc_env%rho_air(t) = air_density(loc_env%Psat(t),loc_env%Rh(t),env_dat%P,loc_env%T_sys_K(t))
    !wet-bulb temperature
    loc_env%T_wb_C(t) = wet_bulb_temperature(loc_env%T_sys_C(t), loc_env%Rh(t)*100._i_kind)
    !film temperature
    loc_env%T_f_C(t) = (loc_env%T_wb_C(t)+loc_env%T_sys_C(t))/2._i_kind
    !saturation vapour pressure at film temperature
    loc_env%Psat_f(t) = saturated_vapor_pressure(loc_env%T_f_C(t))
    !dynamic viscosity of air at film temperature
    loc_env%mu_air_f(t) = sutherland_formula(loc_env%T_f_C(t)+K_0)
    !air density at film temperature
    loc_env%rho_air_f(t) = air_density(loc_env%Psat_f(t),loc_env%Rh(t),env_dat%P,loc_env%T_f_C(t)+K_0)
    !water density at film temperature
    loc_env%rho_h2o_f(t) = water_density(loc_env%T_f_C(t))
    !diffusion coefficient of water at film temperature
    loc_env%D_h2o_f(t) = diffusion_coefficient_water(loc_env%T_f_C(t))
    !Schmidt number at film temperature
    loc_env%Sc_h2o_f(t) = loc_env%mu_air_f(t)/(loc_env%rho_air_f(t)*loc_env%D_h2o_f(t))
    !set time step width
	loc_env%dt(t) = minval(dt_min, mask=ac)

    !!update droplet classes
    do d=1, n_dc
  	  if(ac(d))then
  	    !time
  	    drop_mod(d)%time(t) = loc_env%time(t)
  	    !update mass h2o 
  	    drop_mod(d)%m_h2o(t) = drop_mod(d)%m_h2o(t-1)+drop_mod(d)%dm_h2o(t-1)
  	    if(drop_mod(d)%m_h2o(t)<0._i_kind)then
  	      drop_mod(d)%m_h2o(t) = 0._i_kind
  	    end if
  	    !update mass AI
  	    drop_mod(d)%m_AI(t) = drop_mod(d)%m_AI(t-1)+drop_mod(d)%dm_AI(t-1)
  	    if(drop_mod(d)%m_AI(t)<0._i_kind)then
  	      drop_mod(d)%m_AI(t) = 0._i_kind
  	    end if
  	    !update volume
  	    drop_mod(d)%volume(t) = (drop_mod(d)%m_h2o(t)/loc_env%rho_h2o_f(t)) + (drop_mod(d)%m_AI(t)/app_dat%rho_AI)
  	    !update diameter
  	    drop_mod(d)%diameter(t) = (6._i_kind*drop_mod(d)%volume(t)/pi)**(1._i_kind/3._i_kind)
  	    !update fraction of AI
  	    drop_mod(d)%fract_AI(t) = drop_mod(d)%m_AI(t)/drop_mod(d)%m_AI(1)
  	    !update molar fraction of water
  	    drop_mod(d)%x_h2o(t) = (drop_mod(d)%m_h2o(t)/mm_h2o)/((drop_mod(d)%m_h2o(t)/mm_h2o)+(drop_mod(d)%m_AI(t)/app_dat%mm_AI))
  	    !update x-position
  	    drop_mod(d)%x(t) = drop_mod(d)%x(t-1) + (drop_mod(d)%V_x(t-1)*loc_env%dt(t-1))
  	    !update z-position
  	    drop_mod(d)%z(t) = drop_mod(d)%z(t-1) - (drop_mod(d)%V_z(t-1)*loc_env%dt(t-1))
  	    !update x-velocity
  	    if(abs(env_dat%U_wind-drop_mod(d)%V_x(t-1)) < 0.0001_i_kind)then
		  drop_mod(d)%V_x(t) = env_dat%U_wind
		else
		  drop_mod(d)%V_x(t) = drop_mod(d)%V_x(t-1) + (drop_mod(d)%dV_x_dt(t-1)*loc_env%dt(t-1))
		end if
		!update z-velocity
  	    if(vv(d))then
		  drop_mod(d)%V_z(t) = drop_mod(d)%V_z(t-1) + (drop_mod(d)%dV_z_dt(t-1)*loc_env%dt(t-1))
		else
		  drop_mod(d)%V_z(t) = g*drop_mod(d)%tau_m(t-1)
		end if
		!Reynolds number at system temperature
  	    drop_mod(d)%Re_s(t) = (drop_mod(d)%diameter(t)*loc_env%rho_air(t)*(((env_dat%U_wind-drop_mod(d)%V_x(t))**2._i_kind+(drop_mod(d)%V_z(t))**2._i_kind)**0.5_i_kind))/(loc_env%mu_air(t))
  	    !Boothroyd's f
  	    drop_mod(d)%f(t) = 1.0_i_kind+0.15_i_kind*drop_mod(d)%Re_s(t)**0.687_i_kind
		!tau_m
  	    drop_mod(d)%tau_m(t) = (((drop_mod(d)%m_AI(t)+drop_mod(d)%m_h2o(t))/drop_mod(d)%volume(t))*drop_mod(d)%diameter(t)**2._i_kind)/(18._i_kind*loc_env%mu_air(t)*drop_mod(d)%f(t))
  	    !Reynolds number at film temperature
  	    drop_mod(d)%Re_f(t) = (drop_mod(d)%diameter(t)*loc_env%rho_air_f(t)*(((env_dat%U_wind-drop_mod(d)%V_x(t))**2._i_kind+(drop_mod(d)%V_z(t))**2._i_kind)**0.5_i_kind))/(loc_env%mu_air_f(t))
  	    !Sherwood number h2o
  	    drop_mod(d)%Sh_h2o_f(t) = 1.0_i_kind+0.276_i_kind*drop_mod(d)%Re_f(t)**0.5_i_kind*loc_env%Sc_h2o_f(t)**(1._i_kind/3._i_kind)
  	    !flux of water
  	    drop_mod(d)%F_h2o(t) = ((drop_mod(d)%Sh_h2o_f(t)*loc_env%D_h2o_f(t)*mm_h2o)/(drop_mod(d)%diameter(t)*R*(loc_env%T_f_C(t)+K_0)))*((loc_env%Psat(t)*loc_env%Rh(t))-(loc_env%Psat_f(t)*drop_mod(d)%x_h2o(t))) * evap_fact
  	    !flux of AI
  	    drop_mod(d)%F_AI(t) = -4.07E-7_i_kind*app_dat%P_AI*app_dat%mm_AI*(1._i_kind-drop_mod(d)%x_h2o(t))
  	    !rate of change of x-velocity
		if(abs(env_dat%U_wind-drop_mod(d)%V_x(t)) < 0.0001_i_kind)then
		  drop_mod(d)%dV_x_dt(t) = 0._i_kind
		else
		  drop_mod(d)%dV_x_dt(t) = (env_dat%U_wind-drop_mod(d)%V_x(t))/drop_mod(d)%tau_m(t)
		end if
  	    !rate of change of z-velocity
  	    drop_mod(d)%dV_z_dt(t) = g -(drop_mod(d)%V_z(t)/drop_mod(d)%tau_m(t))
  	    !rate of change of water mass
  	    drop_mod(d)%dm_h2o_dt(t) = pi*drop_mod(d)%diameter(t)**2._i_kind*drop_mod(d)%F_h2o(t)
  	    !rate of change of AI mass
  	    drop_mod(d)%dm_AI_dt(t) = pi*drop_mod(d)%diameter(t)**2._i_kind*drop_mod(d)%F_AI(t)
  	    !settling velocity
		!drop_mod(d)%V_s(t) = g*drop_mod(d)%tau_m(t)
		drop_mod(d)%V_s(t) = settling_velocity (drop_mod(d)%diameter(t),drop_mod(d)%V_s(t-1),loc_env%rho_air(t),loc_env%mu_air(t),((drop_mod(d)%m_AI(t)+drop_mod(d)%m_h2o(t))/drop_mod(d)%volume(t)))
		!effective deposition velocity
  	    drop_mod(d)%edv(t) = (drop_mod(d)%V_z(t))*drop_mod(d)%fract_AI(t)
		! sigma_v calculation
	    if(drop_mod(d)%z(t) > 0._i_kind)then
		  !Souton 2001
		  drop_mod(d)%sigma_v(t) = souton_sigma_v(env_dat%sig_v,drop_mod(d)%time(t))
	    else
		  !Souton 2001
		  drop_mod(d)%sigma_v(t) = drop_mod(d)%sigma_v(t-1)
	    end if
		!calculate distribution above ground
        drop_mod(d)%dag(t) = 1._i_kind-cdf_normal_sigma(drop_mod(d)%sigma_v(t),drop_mod(d)%z(t),0._i_kind)
		!calculate distribution above ground for resolution
        drop_mod(d)%dag_res(t) = 1._i_kind-cdf_normal_sigma(drop_mod(d)%sigma_v(t),drop_mod(d)%z(t),control_dat%z_0)
	    !calculate distribution above ground for resolution
        drop_mod(d)%dag_res_chg(t) = drop_mod(d)%dag_res(t)-drop_mod(d)%dag_res(t-1)
		
		! puff slices
    	do s = 1, N_slices
          ! vertical position of slice
          drop_mod(d)%puff_slc%z(t,s) = drop_mod(d)%z(t) + (slice_multi(s) * drop_mod(d)%sigma_v(t))
      
          ! wind speed at slice height
          drop_mod(d)%puff_slc%Uz(t,s) = wind_profile(drop_mod(d)%puff_slc%z(t,s))
      
          ! initialize travelled distance (none yet)
          drop_mod(d)%puff_slc%x_trav(t,s) = drop_mod(d)%puff_slc%x_trav(t-1,s) + drop_mod(d)%puff_slc%Uz(t,s) * loc_env%dt(t-1)
      
          ! horizontal dispersion at t=0 (zero or small epsilon)
          drop_mod(d)%puff_slc%sigma_h_out(t,s) = drop_mod(d)%puff_slc%sigma_h_out(t-1,s) + (souton_sigma_h(env_dat%sig_h,drop_mod(d)%time(t),drop_mod(d)%puff_slc%x_trav(t,s))-souton_sigma_h(env_dat%sig_h,drop_mod(d)%time(t-1),drop_mod(d)%puff_slc%x_trav(t-1,s)))
        end do
		
		! indexing for interpolation step
    	drop_mod(d)%z_idx(t,1) = max(1, min(maxloc(drop_mod(d)%puff_slc%z(t,:), dim=1, mask=drop_mod(d)%puff_slc%z(t,:) <= control_dat%z_0), N_slices-1))
        drop_mod(d)%z_idx(t,2) = drop_mod(d)%z_idx(t,1) + 1
		
		! update puff dimensions
        ! bottom (first contact)
        if (abs(drop_mod(d)%puff_dim%t1) < 1e-6_i_kind) then
            if (drop_mod(d)%puff_slc%z(t,1) < control_dat%z_0) then
                drop_mod(d)%puff_dim%t1 = drop_mod(d)%time(t)
            end if
        end if
        
        ! assigne x and y limits
		drop_mod(d)%puff_dim%x_min = min( drop_mod(d)%puff_dim%x_min, minval(drop_mod(d)%puff_slc%x_trav(t,(/drop_mod(d)%z_idx(t,1),drop_mod(d)%z_idx(t,2)/)) - drop_mod(d)%puff_slc%sigma_h_out(t,(/drop_mod(d)%z_idx(t,1),drop_mod(d)%z_idx(t,2)/)) * 3._i_kind))
        drop_mod(d)%puff_dim%x_max = max( drop_mod(d)%puff_dim%x_max, maxval(drop_mod(d)%puff_slc%x_trav(t,(/drop_mod(d)%z_idx(t,1),drop_mod(d)%z_idx(t,2)/)) + drop_mod(d)%puff_slc%sigma_h_out(t,(/drop_mod(d)%z_idx(t,1),drop_mod(d)%z_idx(t,2)/)) * 3._i_kind))
        drop_mod(d)%puff_dim%y_max = max( drop_mod(d)%puff_dim%y_max, maxval(drop_mod(d)%puff_slc%sigma_h_out(t,(/drop_mod(d)%z_idx(t,1),drop_mod(d)%z_idx(t,2)/)) * 3._i_kind))

        !determine time step width
		if(vv(d)) then
		  !calculate dt's
		  dt_vel(d) = abs(drop_mod(d)%V_z(t)/drop_mod(d)%dV_z_dt(t))/fract_limit
		  dt_mass(d) = abs(drop_mod(d)%m_h2o(t) / drop_mod(d)%dm_h2o_dt(t))/fract_limit
		  if(drop_mod(d)%V_x(t) < 1e-3)then
		    dt_dist(d) = 1._i_kind
		  else
		    dt_dist(d) = control_dat%max_dist_step / drop_mod(d)%V_x(t)
		  end if
		  !check if velocity is still different from direct solution
		  if(abs((drop_mod(d)%V_z(t)-g*drop_mod(d)%tau_m(t))/drop_mod(d)%V_z(t)) < 0.001_i_kind)then
		    vv(d) = .false.
		  end if
		  
		  !set minimal dt
		  dt_min(d) = max(minval((/dt_vel(d),dt_mass(d),dt_dist(d)/)),1e-6)
		else
		  !calculate dt_mass
		  dt_mass(d) = abs(drop_mod(d)%m_h2o(t) / drop_mod(d)%dm_h2o_dt(t))/fract_limit
		  
		  !set minimal dt
		  dt_min(d) = max(minval((/dt_mass(d),dt_dist(d)/)),1e-6)
		end if
		
		!check for exit conditions
		if (t == s_max) then
		  exit_code(d) = 1
		else if (drop_mod(d)%time(t) > control_dat%max_time) then
		  exit_code(d) = 2
		else if (drop_mod(d)%puff_slc%z(t,N_slices) <= control_dat%z_0) then
		  exit_code(d) = 3
		else if (minval(drop_mod(d)%puff_slc%x_trav(t,:) - drop_mod(d)%puff_slc%sigma_h_out(t,:) * 3._i_kind) > control_dat%max_dist) then
		  exit_code(d) = 4
		end if
		
		!execute exit logic
		if (exit_code(d) /= 0) then
		  ac(d) = .false.
  	  	  iter_vec(d) = t
		  drop_mod(d)%puff_dim%t2 = drop_mod(d)%time(t-1)
		end if
		
	  end if
    end do
    
    
    !!calculate evaporating masses
    do d=1, n_dc
      if(ac(d))then
        !change in mass of h2o
        drop_mod(d)%dm_h2o(t) = drop_mod(d)%dm_h2o_dt(t) * loc_env%dt(t)
		!calculate x_mol based on psat,Rh and psatwb
		xmol = loc_env%Psat(t)*loc_env%Rh(t)/loc_env%Psat_f(t)
		!calculate mh2o based on psat,Rh and psatwb
		mh2o = (-(xmol*(drop_mod(d)%m_AI(t)/app_dat%mm_AI))/(xmol-1._i_kind))*mm_h2o
		if(drop_mod(d)%m_h2o(t)+drop_mod(d)%dm_h2o(t) < mh2o)then
		  drop_mod(d)%dm_h2o(t) = mh2o - drop_mod(d)%m_h2o(t)
		end if
		!change in mass of AI
        drop_mod(d)%dm_AI(t) = drop_mod(d)%dm_AI_dt(t) * loc_env%dt(t)
        !calculate AI loss to atmosphere
		drop_mod(d)%AI_a(t) = drop_mod(d)%AI_a(t) + abs(drop_mod(d)%dm_AI(t))
		if(drop_mod(d)%dag(t) > 0._i_kind)then
  	      !calculate h2o mass added to ambient system
          drop_mod(d)%dm_h2o_a(t) = -drop_mod(d)%dm_h2o(t)*vol_fact(d)*drop_mod(d)%dag(t)
  	      !summ mass
          loc_env%m_h2o_drop(t) = loc_env%m_h2o_drop(t) + drop_mod(d)%m_h2o(t)*vol_fact(d)*drop_mod(d)%dag(t)
          !summ evaporated mass
          loc_env%d_m_h2o_evap(t) = loc_env%d_m_h2o_evap(t) + drop_mod(d)%dm_h2o_a(t)
  	    else
  	      !calculate h2o mass added to ambient system
          drop_mod(d)%dm_h2o_a(t) = 0._i_kind
  	    end if
  	  end if
    end do
	
	!calculate energy losst due to evaporation
	loc_env%d_E_h2o_evap(t) = -(loc_env%d_m_h2o_evap(t)*latent_heat_of_vaporization(loc_env%T_sys_K(t))) + ((loc_env%T_wb_C(t)+K_0)*loc_env%d_m_h2o_evap(t)*Cp_h2o)
	
	!update t
	t = t+1
  end do
  
  !limit puff_dim to user inputs
  do d=1, n_dc
    ratio = abs(drop_mod(d)%puff_dim%x_min)/abs(drop_mod(d)%puff_dim%x_max)
	if(d == 1)then
	  drop_mod(d)%puff_dim%x_min = max(drop_mod(d)%puff_dim%x_min, (-control_dat%max_dist*ratio))
	else
	  drop_mod(d)%puff_dim%x_min = max(drop_mod(d)%puff_dim%x_min, (-control_dat%max_dist*ratio), drop_mod(d-1)%puff_dim%x_min)
	end if
	drop_mod(d)%puff_dim%x_max = min(drop_mod(d)%puff_dim%x_max, control_dat%max_dist)
	drop_mod(d)%puff_dim%y_max = min(drop_mod(d)%puff_dim%y_max, (control_dat%max_dist/2._i_kind))	
  end do
  
  !limit data point count in droplet models to last valid execution
  !drop model
  do d=1, n_dc
  sel_vec = (/(i,i=1,(min(iter_vec(d),t-1)))/)
  drop_mod(d)%time = drop_mod(d)%time(sel_vec)
  drop_mod(d)%diameter = drop_mod(d)%diameter(sel_vec)
  drop_mod(d)%volume = drop_mod(d)%volume(sel_vec)
  drop_mod(d)%m_h2o = drop_mod(d)%m_h2o(sel_vec)
  drop_mod(d)%m_AI = drop_mod(d)%m_AI(sel_vec)
  drop_mod(d)%fract_AI = drop_mod(d)%fract_AI(sel_vec)
  drop_mod(d)%x_h2o = drop_mod(d)%x_h2o(sel_vec)
  drop_mod(d)%x = drop_mod(d)%x(sel_vec)
  drop_mod(d)%z = drop_mod(d)%z(sel_vec)
  drop_mod(d)%V_x = drop_mod(d)%V_x(sel_vec)
  drop_mod(d)%V_z = drop_mod(d)%V_z(sel_vec)
  drop_mod(d)%Re_s = drop_mod(d)%Re_s(sel_vec)
  drop_mod(d)%f = drop_mod(d)%f(sel_vec)
  drop_mod(d)%tau_m = drop_mod(d)%tau_m(sel_vec)
  drop_mod(d)%Re_f = drop_mod(d)%Re_f(sel_vec)
  drop_mod(d)%Sh_h2o_f = drop_mod(d)%Sh_h2o_f(sel_vec)
  drop_mod(d)%F_h2o = drop_mod(d)%F_h2o(sel_vec)
  drop_mod(d)%F_AI = drop_mod(d)%F_AI(sel_vec)
  drop_mod(d)%dV_x_dt = drop_mod(d)%dV_x_dt(sel_vec)
  drop_mod(d)%dV_z_dt = drop_mod(d)%dV_z_dt(sel_vec)
  drop_mod(d)%dm_h2o_dt = drop_mod(d)%dm_h2o_dt(sel_vec)
  drop_mod(d)%dm_AI_dt = drop_mod(d)%dm_AI_dt(sel_vec)
  drop_mod(d)%dm_h2o = drop_mod(d)%dm_h2o(sel_vec)
  drop_mod(d)%dm_AI = drop_mod(d)%dm_AI(sel_vec)
  drop_mod(d)%AI_a = drop_mod(d)%AI_a(sel_vec)
  drop_mod(d)%dm_h2o_a = drop_mod(d)%dm_h2o_a(sel_vec)
  drop_mod(d)%V_s = drop_mod(d)%V_s(sel_vec)
  drop_mod(d)%edv = drop_mod(d)%edv(sel_vec)
  drop_mod(d)%dag = drop_mod(d)%dag(sel_vec)
  drop_mod(d)%dag_res = drop_mod(d)%dag_res(sel_vec)
  drop_mod(d)%dag_res_chg = drop_mod(d)%dag_res_chg(sel_vec)
  drop_mod(d)%dt = drop_mod(d)%dt(sel_vec)
  drop_mod(d)%sigma_v = drop_mod(d)%sigma_v(sel_vec)
  drop_mod(d)%puff_slc%z = drop_mod(d)%puff_slc%z(sel_vec,:)
  drop_mod(d)%puff_slc%Uz = drop_mod(d)%puff_slc%Uz(sel_vec,:)
  drop_mod(d)%puff_slc%x_trav = drop_mod(d)%puff_slc%x_trav(sel_vec,:)
  drop_mod(d)%puff_slc%sigma_h_out = drop_mod(d)%puff_slc%sigma_h_out(sel_vec,:)
  end do
  
  !local environment
  sel_vec = (/(i,i=1,(t-1))/)
  loc_env%time = loc_env%time(sel_vec)
  loc_env%volume = loc_env%volume(sel_vec)
  loc_env%d_volume = loc_env%d_volume(sel_vec)
  loc_env%m_air = loc_env%m_air(sel_vec)
  loc_env%d_m_air = loc_env%d_m_air(sel_vec)
  loc_env%m_h2o = loc_env%m_h2o(sel_vec)
  loc_env%d_m_h2o = loc_env%d_m_h2o(sel_vec)
  loc_env%d_E_air = loc_env%d_E_air(sel_vec)
  loc_env%d_E_h2o = loc_env%d_E_h2o(sel_vec)
  loc_env%E_tot_t1 = loc_env%E_tot_t1(sel_vec)
  loc_env%m_h2o_drop = loc_env%m_h2o_drop(sel_vec)
  loc_env%d_m_h2o_evap = loc_env%d_m_h2o_evap(sel_vec)
  loc_env%d_E_h2o_evap = loc_env%d_E_h2o_evap(sel_vec)
  loc_env%T_sys_K = loc_env%T_sys_K(sel_vec)
  loc_env%T_sys_C = loc_env%T_sys_C(sel_vec)
  loc_env%Psat = loc_env%Psat(sel_vec)
  loc_env%Rh = loc_env%Rh(sel_vec)
  loc_env%mu_air = loc_env%mu_air(sel_vec)
  loc_env%rho_air = loc_env%rho_air(sel_vec)
  loc_env%T_wb_C = loc_env%T_wb_C(sel_vec)
  loc_env%T_f_C = loc_env%T_f_C(sel_vec)
  loc_env%Psat_f = loc_env%Psat_f(sel_vec)
  loc_env%mu_air_f = loc_env%mu_air_f(sel_vec)
  loc_env%rho_air_f = loc_env%rho_air_f(sel_vec)
  loc_env%rho_h2o_f = loc_env%rho_h2o_f(sel_vec)
  loc_env%D_h2o_f = loc_env%D_h2o_f(sel_vec)
  loc_env%Sc_h2o_f = loc_env%Sc_h2o_f(sel_vec)
  loc_env%dt = loc_env%dt(sel_vec)
  
  print *, 'droplet model run successful'
  
  end subroutine droplet_model