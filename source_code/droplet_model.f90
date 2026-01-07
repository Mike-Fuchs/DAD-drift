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
  
  real(i_kind) :: dt_min														![s]					|minimal time step width
  real(i_kind) :: dt_tmp														![s]					|temporal time step width
  real(i_kind) :: z_min															![s]					|threshold when to stop execution
  real(i_kind), dimension(:), allocatable :: vol_fact							![-]					|vector of volume factors for each droplet class
  real(i_kind), dimension(:), allocatable :: alpha								![-]					|sd vector for vertical puff dispersion
  real(i_kind), dimension(:), allocatable :: z									![m]					|vector of vertical position
  real(i_kind), dimension(:), allocatable :: U									![m/s]					|vector of wind speed
  real(i_kind) :: xmol															![-]					|molar fraction
  real(i_kind) :: mh2o															![-]					|mass of water
  real(i_kind) :: evap_fact														![-]					|evaporation fraction
  integer, dimension(:), allocatable :: sel_vec									![-]					|vector for data selection
  integer, dimension(:), allocatable :: iter_vec								![-]					|vector of maximal iteration steps
  integer :: s_max																![-]					|maximal step count
  integer :: n_dc																![-]					|number of droplet classes
  integer :: t																	![-]					|counter looping over time
  integer :: d																	![-]					|counter looping over droplet classes
  integer :: i																	![-]					|counter
  logical, dimension(:), allocatable :: ac										![-]					|logical vector if class is activly simulated, droplets z-position is above a given threshold
  logical, dimension(:), allocatable :: vv										![-]					|logical vector if velocity is deviating from ideal speed
  
  !set evap_fact
  evap_fact = 1.0_i_kind
  
  !set minimal time step and maximal step count
  dt_min = 0.001_i_kind
  s_max = 500000
  z_min = -10._i_kind
  
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
  
  !allocate interation count vector
  allocate(iter_vec(n_dc))
  iter_vec = s_max
  
  !define alpha vector
  alpha = real([(i,i=-10,10)],i_kind)/2._i_kind
  allocate(z(size(alpha)))
  allocate(U(size(alpha)))
  
  !allocate droplet model and loval environment
  call droplet_model_allo(s_max, n_dc,size(alpha))
  
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
  loc_env%dt(1) = dt_min
  
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
	drop_mod(d)%sigma_h(1) = 0._i_kind
	drop_mod(d)%sigma_v(1) = 0._i_kind
	!calculate distribution above ground
    drop_mod(d)%dag(1) = 1._i_kind
    !calculate distribution above ground for resolution
    drop_mod(d)%dag_res(1) = 1._i_kind
    !
    drop_mod(d)%dag_res_chg(1) = 0._i_kind
	!x_s
	drop_mod(d)%x_s(1,:) = 0._i_kind
	
    
	!!time step width estimation
	if (.not. abs(drop_mod(d)%V_z(1)) > 0._i_kind) then
	  dt_tmp = dt_min
	else
	  dt_tmp = (drop_mod(d)%V_z(1)/abs(drop_mod(d)%dV_z_dt(1)))/10._i_kind
	end if
	if(.not. abs(dt_tmp) > 0._i_kind) dt_tmp = dt_min
	drop_mod(d)%dt(1) = dt_tmp
	    
	!delta time
	if(drop_mod(d)%dt(1) < loc_env%dt(1))then
	  loc_env%dt(1) = drop_mod(d)%dt(1)
	end if
  end do
	
	
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
  
  
  !!calculate wind speed in dependency of z
  do d=1, n_dc
    !calculate wind speed
	if(drop_mod(d)%z(1) >= control_dat%z_0)then
	  !wind speed
	  if(drop_mod(d)%z(1) >= env_dat%Hv)then
	    drop_mod(d)%Uz(1) = (env_dat%U_fric/Kc)*log((drop_mod(d)%z(1)-env_dat%d)/env_dat%z0)
	  else
	    drop_mod(d)%Uz(1) = env_dat%Uh*exp((env_dat%LAI/2._i_kind)*((drop_mod(d)%z(1)/env_dat%Hv)-1._i_kind))
	  end if
	  !travelled distance
	  drop_mod(d)%Ud(1) = drop_mod(d)%Uz(1) * loc_env%dt(1)
	end if
  end do
  
  !calculate energy losst due to evaporation
  loc_env%d_E_h2o_evap(1) = -(loc_env%d_m_h2o_evap(1)*latent_heat_of_vaporization(loc_env%T_sys_K(1))) + ((loc_env%T_wb_C(1)+K_0)*loc_env%d_m_h2o_evap(1)*Cp_h2o)
  
  
  !!time loop
  t = 2
  do while (all((/any(ac),loc_env%time(t-1)<control_dat%max_time,t<=s_max/)))
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
	loc_env%dt(t) = dt_min
	
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
  	    !compute Uz
		if(drop_mod(d)%z(t) > 0._i_kind)then
		  drop_mod(d)%Uz(t) = wind_profile(drop_mod(d)%z(t))
		else
		  drop_mod(d)%Uz(t) = 0._i_kind
		end if
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
		
  	    if(vv(d))then
  	      !!time step width estimation
  	      if (.not. abs(drop_mod(d)%V_z(t)) > 0._i_kind) then
	        dt_tmp = dt_min
		  else
		    dt_tmp = (drop_mod(d)%V_z(t)/abs(drop_mod(d)%dV_z_dt(t)))/1000._i_kind
		  end if
  	      if(.not. abs(dt_tmp) > 0._i_kind) dt_tmp = dt_min
  	      drop_mod(d)%dt(t) = dt_tmp
		  
  	      !delta time
		  if(drop_mod(d)%dt(t) < loc_env%dt(t))then
  	        loc_env%dt(t) = drop_mod(d)%dt(t)
  	      end if
		  
		  !check if velocity is still different from direct solution
		  if(abs((drop_mod(d)%V_z(t)-g*drop_mod(d)%tau_m(t))/drop_mod(d)%V_z(t)) < 0.001_i_kind)then
		    vv(d) = .false.
		  end if
		end if
  	    
  	    !check if z-position if still above threshold
  	    if(all((/(drop_mod(d)%dag_res(t-1) < 0.0001_i_kind)/)))then
  	      ac(d) = .false.
  	  	  iter_vec(d) = t
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
	
	
	!!calculate wind speed in dependency of z
    do d=1, n_dc
      !calculate wind speed
	  if(drop_mod(d)%z(t) >= control_dat%z_0)then
	    !wind speed
	    if(drop_mod(d)%z(t) >= env_dat%Hv)then
	      drop_mod(d)%Uz(t) = (env_dat%U_fric/Kc)*log((drop_mod(d)%z(t)-env_dat%d)/env_dat%z0)
	    else
	      drop_mod(d)%Uz(t) = env_dat%Uh*exp((env_dat%LAI/2._i_kind)*((drop_mod(d)%z(t)/env_dat%Hv)-1._i_kind))
	    end if
	    !travelled distance
	    drop_mod(d)%Ud(t) = drop_mod(d)%Uz(t) * loc_env%dt(t)
	  end if
    end do
	
	!determine x dependent of z_0
	do d=1, n_dc
	  !reset vectors
	  z = 0._i_kind
	  U = 0._i_kind
	  
	  !determine z position
	  z = drop_mod(d)%z(t) + alpha * drop_mod(d)%sigma_v(t)
	  
	  !determine wind speed
	  do i=1, size(z)
	    if(z(i) > 0._i_kind)then
		  U(i) = wind_profile(z(i))
		else
		  U(i) = wind_profile(0._i_kind)
		end if
	  end do
	  
	  !compute x
	  drop_mod(d)%x_s(t,:) = drop_mod(d)%x_s(t-1,:) + U*loc_env%dt(t)
	  
	  !interpolate x_mean
	  if(any(z < control_dat%z_0) .and. any(z > control_dat%z_0))then
	    drop_mod(d)%x_mean(t) = interpolate(z,drop_mod(d)%x_s(t,:),control_dat%z_0)
	  else
	    drop_mod(d)%x_mean(t) = 0._i_kind
	  end if
	  !compute sigma_h
	  drop_mod(d)%sigma_h(t) = souton_sigma_h(env_dat%sig_h,loc_env%time(t),(loc_env%time(t)*env_dat%U_wind))
	end do	
	
	!calculate energy losst due to evaporation
	loc_env%d_E_h2o_evap(t) = -(loc_env%d_m_h2o_evap(t)*latent_heat_of_vaporization(loc_env%T_sys_K(t))) + ((loc_env%T_wb_C(t)+K_0)*loc_env%d_m_h2o_evap(t)*Cp_h2o)
	
	!increase time step width
	if(loc_env%time(t) >= 300._i_kind)then
	  dt_min = 1._i_kind
	else if(loc_env%time(t) >= 60._i_kind)then
	  dt_min = 0.1_i_kind
	else if(loc_env%time(t) >= 5._i_kind)then
	  dt_min = 0.01_i_kind
	end if	
	
	!update t
	t = t+1
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
  drop_mod(d)%Uz = drop_mod(d)%Uz(sel_vec)
  drop_mod(d)%Ud = drop_mod(d)%Ud(sel_vec) 
  drop_mod(d)%x_mean = drop_mod(d)%x_mean(sel_vec)
  drop_mod(d)%sigma_h = drop_mod(d)%sigma_h(sel_vec)
  drop_mod(d)%sigma_v = drop_mod(d)%sigma_v(sel_vec)
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