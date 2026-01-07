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

module data_module
  !! ~~~~ description ~~~~
  !! In this module various data frames are defined.
  use constants_module
    
  implicit none
  private
  public :: env_dat, drop_spec, control_dat, app_dat, drop_mod, loc_env, drift_curve_matrix
  
  type :: environmental_data
    real(i_kind) :: T_C																![°C]					|ambient temperature
	real(i_kind) :: Rh																![-]					|relative humidity
	real(i_kind) :: U_wind															![m/s]					|wind speed
    real(i_kind) :: D_wind															![°]					|wind direction
	real(i_kind) :: P																![Pa]					|ambient pressure
	real(i_kind) :: T_K																![K]					|ambient temperature
    real(i_kind) :: Rh_P															![%]					|relative humidity
	real(i_kind) :: mu_air															![kg/ms]				|dynamic viscosity of air at ambient temperature
	real(i_kind) :: Psat															![Pa]					|saturation vapor pressure at ambient temperature
	real(i_kind) :: rho_air															![kg/m³]				|air density
	real(i_kind) :: rho_h2o															![kg/m³]				|density of water at ambient temperature
	real(i_kind) :: ah																![kg/m³]				|absolute humidity at ambient temperature
	real(i_kind) :: sig_h															![m]					|horizontal eddy diffusivity
	real(i_kind) :: sig_v															![m]					|vertical eddy diffusivity
	real(i_kind) :: z0																![m]					|roughness height
	real(i_kind) :: z_wind															![m]					|reference height of mean wind velocity
	real(i_kind) :: U_fric															![m]					|friction velocity
	real(i_kind) :: Hv																![m]					|height of vegetation (canopy)
	real(i_kind) :: LAI																![m²/m²]				|leaf area index of vegetation
	real(i_kind) :: d																![m]					|zero plane displacement
	real(i_kind) :: Uh																![m/s]					|wind speed at the top of the canopy
	end type environmental_data
  type(environmental_data) :: env_dat
  
   
  type :: droplet_spectrum
    real(i_kind), dimension (:), allocatable :: ds									![m]					|droplet size
	real(i_kind), dimension (:), allocatable :: cf									![-]					|cummulative fraction of droplet spectrum
	real(i_kind), dimension (:), allocatable :: f									![-]					|fraction of droplet spectrum
  end type droplet_spectrum
  type(droplet_spectrum) :: drop_spec
  
  
  real(i_kind), dimension(:,:), allocatable :: drift_curve_matrix
  
  
  type control_data
    integer :: code1																!						|code controlling executed routines (1=drift curve, 0=drift raster)
	real(i_kind) :: z_0																![m]					|deposition height
	real(i_kind) :: max_dist														![m]					|maximum distance of drift calculation
	real(i_kind) :: max_time = 3600													![s]					|maximum distance of drift calculation
	integer :: max_n_time = 128														![-]					|maximum number of time steps
	integer :: field_count															![-]					|number of fields
	real(i_kind) :: cellsize														![m]					|cellsize of the raster
	real(i_kind) :: sd = 3._i_kind													![-]					|number of standard deviation relevant for pattern prediction
	logical :: dev_mode	= .false.													![-]					|flag if devmode should be active
	end type control_data
  type(control_data) control_dat
  
  
  type application_data
    real(i_kind) :: v_trac															![m/s]					|speed of tractor
	real(i_kind) :: b_width															![m]					|boom width
	real(i_kind) :: b_height														![m]					|boom height
	real(i_kind) :: V_initial														![m/s]					|initial droplet speed
	real(i_kind) :: V_horizontal													![m/s]					|initial droplet speed in horizontal direction
	real(i_kind) :: V_vertical = 0._i_kind											![m/s]					|initial droplet speed in vertical direction
	real(i_kind) :: app_pressure													![pa]					|application pressure
	real(i_kind) :: app_rate_mh														![m³/h]					|application rate
	real(i_kind) :: app_rate_mha													![m³/ha]				|application rate
	real(i_kind) :: app_rate														![kg/ha]				|application rate
	real(i_kind) :: nozzle_angle													![°]					|angle of the spray nozzle
    real(i_kind) :: c_solution														![kg/m³]				|concentration of the spray solution
    real(i_kind) :: rho_AI															![kg/m³]				|density of active ingredient
    real(i_kind) :: rho_solution													![kg/m³]				|density of the spray solution
	real(i_kind) :: mm_AI															![kg/mol]				|molar mass of active ingredient
    real(i_kind) :: P_AI															![Pa]					|vapor pressure of active ingredient
    integer :: n_pass																![-]					|number of passes during drift curve experiment
	real(i_kind) :: f_length														![m]					|length of the field
  end type application_data
  type(application_data) :: app_dat
  
  
  type droplet_model
    real(i_kind), dimension (:), allocatable :: time								![s]					|time
	real(i_kind), dimension (:), allocatable :: diameter							![m]					|droplet diameter
	real(i_kind), dimension (:), allocatable :: volume								![m³]					|droplet volume
	real(i_kind), dimension (:), allocatable :: m_h2o								![kg]					|mass of water
	real(i_kind), dimension (:), allocatable :: m_AI								![kg]					|mass of AI
	real(i_kind), dimension (:), allocatable :: fract_AI							![-]					|fraction of initial AI mass
	real(i_kind), dimension (:), allocatable :: x_h2o								![-]					|molar fraction of water
	real(i_kind), dimension (:), allocatable :: x									![m]					|x-position
	real(i_kind), dimension (:), allocatable :: z									![m]					|z-position
	real(i_kind), dimension (:), allocatable :: V_x									![m/s]					|x-velocity
	real(i_kind), dimension (:), allocatable :: V_z									![m/s]					|z-velocity
	real(i_kind), dimension (:), allocatable :: Re_s								![-]					|Reynolds number at system temperature
	real(i_kind), dimension (:), allocatable :: f									![-]					|Boothroyd's f
	real(i_kind), dimension (:), allocatable :: tau_m								![s]					|time constant
	real(i_kind), dimension (:), allocatable :: Re_f								![-]					|Reynolds number at film temperature
	real(i_kind), dimension (:), allocatable :: Sh_h2o_f							![-]					|Sherwood number of water at film temperature
	real(i_kind), dimension (:), allocatable :: F_h2o								![kg/m²s]				|Flux of water
	real(i_kind), dimension (:), allocatable :: F_AI								![kg/m²s]				|Flux of AI
	real(i_kind), dimension (:), allocatable :: dV_x_dt								![m/s²]					|rate of change of x-velocity
	real(i_kind), dimension (:), allocatable :: dV_z_dt								![m/s²]					|rate of change of z-velocity
	real(i_kind), dimension (:), allocatable :: dm_h2o_dt							![kg/s]					|rate of change of water mass
	real(i_kind), dimension (:), allocatable :: dm_AI_dt							![kg/s]					|rate of change of AI mass
	real(i_kind), dimension (:), allocatable :: dm_h2o								![kg]					|change of water mass
	real(i_kind), dimension (:), allocatable :: dm_AI								![kg]					|change of AI mass
	real(i_kind), dimension (:), allocatable :: AI_a								![kg]					|AI loss to atmosphere
	real(i_kind), dimension (:), allocatable :: dm_h2o_a							![kg]					|change of water mass relevant for ambient humidity
	real(i_kind), dimension (:), allocatable :: V_s									![m/s]					|settling velocity
	real(i_kind), dimension (:), allocatable :: edv									![m/s]					|effective deposition velocity
	real(i_kind), dimension (:), allocatable :: dag									![-]					|distribution above ground
	real(i_kind), dimension (:), allocatable :: dag_res								![-]					|distribution above ground
	real(i_kind), dimension (:), allocatable :: dag_res_chg							![-]					|distribution above ground
	real(i_kind), dimension (:), allocatable :: dt									![s]					|delta time
	real(i_kind), dimension (:), allocatable :: Uz									![s]					|wind speed depending on z
	real(i_kind), dimension (:), allocatable :: Ud									![s]					|travelled distance
	real(i_kind), dimension (:), allocatable :: x_mean								![s]					|mean x position of puff at deposition height
	real(i_kind), dimension (:), allocatable :: sigma_h								![s]					|horizontal puff dispersion coefficient
	real(i_kind), dimension (:), allocatable :: sigma_v								![s]					|vertical puff dispersion coefficient
	real(i_kind), dimension (:,:), allocatable :: x_s
  end type droplet_model
  type(droplet_model), dimension(:), allocatable :: drop_mod
  
  
  type local_environment
    real(i_kind), dimension (:), allocatable :: time								![s]					|time
	real(i_kind), dimension (:), allocatable :: volume								![m³]					|air volume
	real(i_kind), dimension (:), allocatable :: d_volume							![m³]					|change in air volume
	real(i_kind), dimension (:), allocatable :: m_air								![kg]					|mass of air
	real(i_kind), dimension (:), allocatable :: d_m_air								![kg]					|change in air mass
	real(i_kind), dimension (:), allocatable :: m_h2o								![kg]					|mass of water in air
	real(i_kind), dimension (:), allocatable :: d_m_h2o								![kg]					|change of water mass in air
	real(i_kind), dimension (:), allocatable :: d_E_air								![J]					|change of the internal energy of air
	real(i_kind), dimension (:), allocatable :: d_E_h2o								![J]					|change of the internal energy of water in air
	real(i_kind), dimension (:), allocatable :: E_tot_t1							![J]					|total internal energy of the volume before latent heat of vaporization is substracted 
	real(i_kind), dimension (:), allocatable :: m_h2o_drop							![kg]					|mass of mater in all droplets
	real(i_kind), dimension (:), allocatable :: d_m_h2o_evap						![kg]					|mass of water evaporating from droplet to athmosphere
	real(i_kind), dimension (:), allocatable :: d_E_h2o_evap						![J]					|change in energie by evaporation
	real(i_kind), dimension (:), allocatable :: T_sys_K								![K]					|temperature of system
	real(i_kind), dimension (:), allocatable :: T_sys_C								![°C]					|temperature of system
	real(i_kind), dimension (:), allocatable :: Psat								![Pa]					|saturation vapour pressure at system temperature
	real(i_kind), dimension (:), allocatable :: Rh									![-]					|realtive humidity
	real(i_kind), dimension (:), allocatable :: mu_air								![kg/ms]				|dynamic viscosity of air at system temperature
	real(i_kind), dimension (:), allocatable :: rho_air								![kg/m³]				|air density at system temperature
	real(i_kind), dimension (:), allocatable :: T_wb_C								![°C]					|wet-bulb temperature
	real(i_kind), dimension (:), allocatable :: T_f_C								![°C]					|film temperature
	real(i_kind), dimension (:), allocatable :: Psat_f								![Pa]					|saturation vapour pressure at film temperature
	real(i_kind), dimension (:), allocatable :: mu_air_f							![kg/ms]				|dynamic viscosity of air at film temperature
	real(i_kind), dimension (:), allocatable :: rho_air_f							![kg/m³]				|air density at film temperature
	real(i_kind), dimension (:), allocatable :: rho_h2o_f							![kg/m³]				|water density at film temperature
	real(i_kind), dimension (:), allocatable :: D_h2o_f								![m²/s]					|diffusion coefficient of water at film temperature
	real(i_kind), dimension (:), allocatable :: Sc_h2o_f							![-]					|Schmitt number of water at film temperature
	real(i_kind), dimension (:), allocatable :: dt									![s]					|time step width
  end type local_environment
  type(local_environment) :: loc_env
  
  end module data_module