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

module functions_module
  !! ~~~~ description ~~~~
  !! This module containes some functions
  
  use constants_module
  
  implicit none
  private
  public :: diffusion_model_xyz_integrated, x_of_ellipse, xyz_create, cell_budget, xy_res, time_res, interpolate, interpolate_idx, interpolate_2d_geo, expand_grid, settling_velocity, wind_profile, sutherland_formula, saturated_vapor_pressure, &
			wet_bulb_temperature, diffusion_coefficient_water, water_density, absolute_humidity, relative_humidity, air_density, latent_heat_of_vaporization, cdf_normal_sigma, souton_sigma_h, souton_sigma_v, souton_v_Kz, &
			unique, position_B_in_A, cbind, rbind
  
  contains
	function diffusion_model_xyz_integrated (x,y,z,M,z_off,x_off,V_z,sigma_h,t,res) result(f)
      !! ~~~~ description ~~~~
	  !! This function will compute the gaussian puff model
	  
	  use data_module
	  use constants_module
		
	  implicit none
	  real(i_kind), intent(in) :: x 					![m]					|x-coordinate
	  real(i_kind), intent(in) :: y 					![m]					|y-coordinate
	  real(i_kind), intent(in) :: z 					![m]					|z-coordinate
	  real(i_kind), intent(in) :: M						![kg]					|released mass
	  real(i_kind), intent(in) :: z_off					![m]					|offset in z-direction as function of t as calculated by droplet model
	  real(i_kind), intent(in) :: x_off					![m]					|mean x as calculated by droplet model
	  real(i_kind), intent(in) :: V_z					![m/s]					|effective deposition velocity
	  real(i_kind), intent(in) :: sigma_h				![m/s]					|effective deposition velocity
	  real(i_kind), intent(in) :: t						![s]					|time
	  real(i_kind), intent(in) :: res					![m]					|resolution of the matrix
	  real(i_kind) :: c									![kg/m³]				|concentration of AI in air
	  real(i_kind) :: f									![kg/m³s]				|deposition rate
	  real(i_kind) :: alpha								![-]					|alpha skewness factor
	  real(i_kind) :: norm_term							![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: x_term							![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: y_term							![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: z_term							![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: skew_term							![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: x1								![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: x2								![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: y1								![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: y2								![-]					|skewed term for implementation into gaussian model
	  real(i_kind) :: sigma_v							![m/s]					|effective deposition velocity
	  real(i_kind) :: Kz							![m/s]					|effective deposition velocity
	  real(i_kind) :: flux_term
	  
	  !determine x1 & x2
	  x1 = x - 0.5_i_kind*res
	  x2 = x + 0.5_i_kind*res
	  
	  !determine y1 & y2
	  y1 = y - 0.5_i_kind*res
	  y2 = y + 0.5_i_kind*res
	  
	  !calculate sigma_v
	  sigma_v = souton_sigma_v(env_dat%sig_v,t)
	  
	  !calculate sigma_v_Kz
	  !Kz = souton_v_Kz(env_dat%sig_v,t)
	  
	  !calculate flux term
	  flux_term = V_z
      
	  !alpha term for skewing
	  alpha = env_dat%k_skew * env_dat%U_wind * t**0.5_i_kind * (sigma_v / sigma_h)
	  
	  !mass normalization term
	  norm_term = M / ((2.0*pi)**1.5 * sigma_h**2 * sigma_v)
	  
	  !x-dimension term
	  x_term = (erf( (x2 - x_off) / (sqrt(2.0_i_kind)*sigma_h) ) - erf( (x1 - x_off) / (sqrt(2.0_i_kind)*sigma_h) )) * (sqrt(pi/2.0)*sigma_h)
	  
	  !y-dimension term
	  y_term = (erf(  y2 / (sqrt(2.0_i_kind)*sigma_h) ) - erf(  y1 / (sqrt(2.0_i_kind)*sigma_h) )) * (sqrt(pi/2.0)*sigma_h)
	  
	  !z-dimension term
	  z_term = exp(-((z-z_off)**2._i_kind)/(2._i_kind*((sigma_v)**2._i_kind)))
	  
	  !skewing term
	  skew_term = 1.0_i_kind + erf( alpha * (x  - x_off) / (sqrt(2.0_i_kind) * sigma_h) )
	  
	  !mass flux
	  f = norm_term * x_term * y_term * z_term * skew_term * flux_term
	  
    end function diffusion_model_xyz_integrated
	
	
	function x_of_ellipse (y,a,b) result(x)
	  !! ~~~~ description ~~~~
	  !! This function will return the two x coordinates of an ellips based on y,a,b
	  
	  implicit none
	  real(i_kind), intent(in) :: y									![-]					|
	  real(i_kind), intent(in) :: a									![-]					|
	  real(i_kind), intent(in) :: b									![-]					|
	  real(i_kind), dimension(2) :: x								![-]					|
	  real(i_kind) :: tmp											![-]					|
	  
	  tmp = (1._i_kind-((y**2._i_kind)/(b**2._i_kind)))**0.5_i_kind
	  
	  x(1) = -a*tmp
	  x(2) = +a*tmp
	  
	end function x_of_ellipse
	
	
	subroutine xyz_create (x_min, x_max, y_max, res, xyz, lat_coord, long_coord)
	  !! ~~~~ description ~~~~
	  !! This function will create the xyz table based on min&max of x&y and the wind direction
	  
	  use data_module
	  use constants_module
	  	  
	  implicit none
	  real(i_kind), intent(in) :: x_min													![-]					|droplet matrix
	  real(i_kind), intent(in) :: x_max													![-]					|droplet matrix
	  real(i_kind), intent(in) :: y_max													![-]					|droplet matrix
	  real(i_kind), intent(in) :: res												![m]					|resoultion
	  real(i_kind), dimension(:,:), allocatable, intent(out) :: xyz					![-]					|xyz table
	  real(i_kind), dimension(:), allocatable, intent(out) :: lat_coord				![m]					|x coordinate
      real(i_kind), dimension(:), allocatable, intent(out) :: long_coord			![m]					|y coordinate
      real(i_kind), dimension(:), allocatable :: x_coord							![m]					|x coordinate
      real(i_kind), dimension(:), allocatable :: y_coord							![m]					|y coordinate
      real(i_kind), dimension(:,:), allocatable :: temp_xy							![-]					|temporal xy matrix
      real(i_kind) :: angle															![°]					|angle of wind direction
      real(i_kind) :: lat_max														![m]					|maximal distance of lat
      real(i_kind) :: lat_min														![m]					|minimal distance of lat
      real(i_kind) :: long_max														![m]					|maximal distance of long
      real(i_kind) :: long_min														![m]					|minimal distance of long
	  real(i_kind) :: lat_max_adj													![m]					|maximal distance of lat adjusted for resolution
      real(i_kind) :: lat_min_adj													![m]					|minimal distance of lat adjusted for resolution
      real(i_kind) :: long_max_adj													![m]					|maximal distance of long adjusted for resolution
      real(i_kind) :: long_min_adj													![m]					|minimal distance of long adjusted for resolution
	  integer :: i																	![-]					|counter
	  
	  !initialize xy grid
      x_coord = real((/(i, i=int((x_min/res)),int((x_max/res)))/),i_kind)*res
      y_coord = real((/(i, i=int((-y_max/res)),int((y_max/res)))/),i_kind)*res
      allocate(temp_xy((size(x_coord)*size(y_coord)),4))
      temp_xy(:,1:2) = expand_grid(x_coord, y_coord)
      
      !rotate to lat long space
      angle = mod(env_dat%D_wind,360._i_kind)
      temp_xy(:,3) = -(cosd(angle)*temp_xy(:,1)+sind(angle)*temp_xy(:,2))
      temp_xy(:,4) = (-sind(angle)*temp_xy(:,1)+cosd(angle)*temp_xy(:,2))
      
      !find lat long upper bounds
      lat_max = real(ceiling(maxval(temp_xy(:,3))),i_kind)
      lat_min = real(floor(minval(temp_xy(:,3))),i_kind)
      long_max = real(ceiling(maxval(temp_xy(:,4))),i_kind)
      long_min = real(floor(minval(temp_xy(:,4))),i_kind)
      
	  !adjust by resolution
	  lat_max_adj = lat_max + 0.5_i_kind - res/2._i_kind
	  lat_min_adj = lat_min - 0.5_i_kind + res/2._i_kind
	  long_max_adj = long_max + 0.5_i_kind - res/2._i_kind
	  long_min_adj = long_min - 0.5_i_kind + res/2._i_kind
	  
	  !initialise lat long space	  
      lat_coord = real((/(i, i=int(lat_min_adj*2._i_kind/res), int(lat_max_adj*2._i_kind/res), 2)/),i_kind) * res / 2._i_kind
      long_coord = real((/(i, i=int(long_min_adj*2._i_kind/res),int(long_max_adj*2._i_kind/res), 2)/),i_kind) * res / 2._i_kind
	  if(res > 1._i_kind)then
	    lat_coord = lat_coord - 0.5_i_kind
		long_coord = long_coord - 0.5_i_kind
	  end if
      allocate(xyz((size(lat_coord)*size(long_coord)),5))
      xyz = 0._i_kind
      xyz(:,1:2) = expand_grid(lat_coord, long_coord)
      
      !rotate coordinate system to inverse lat and long to x and y needed for drift prediction
      xyz(:,3) = -(cosd(angle)*xyz(:,1) + sind(angle)*xyz(:,2))
      xyz(:,4) = (-sind(angle)*xyz(:,1) + cosd(angle)*xyz(:,2))
	  
	  !deallocate arrays
	  if(allocated(x_coord)) deallocate(x_coord)
	  if(allocated(y_coord)) deallocate(y_coord)
	  if(allocated(temp_xy)) deallocate(temp_xy)
	  
	end subroutine xyz_create
	
	function cell_budget(x_span, N_min, N_max, mu, sigma, alpha) result(N_cell)
	  !! ~~~~ description ~~~~
	  !! This function will compute the resolution in x-y space
	  
	  implicit none
	  real(i_kind), intent(in) :: x_span										![-]					|droplet matrix
	  real(i_kind), intent(in) :: N_min											![-]					|droplet matrix
	  real(i_kind), intent(in) :: N_max											![-]					|droplet matrix
	  real(i_kind), intent(in) :: mu											![-]					|droplet matrix
	  real(i_kind), intent(in) :: sigma											![-]					|droplet matrix
	  real(i_kind), intent(in) :: alpha											![-]					|droplet matrix
	  real(i_kind) :: N_cell														![m]					|resolution
	  real(i_kind) :: G, S, arg, N_real
	  
	  ! Gaussian core term
      G = exp( -((x_span - mu)**2) / (2.0_i_kind * sigma**2) )
	  
      ! Skewing factor
      arg = alpha * (x_span - mu) / (sigma * sqrt(2.0_i_kind))
      S   = 1.0_i_kind + erf(arg)
	  
      ! Final continuous N
      N_real = N_min + (N_max - N_min) * G * S
	  
	  ! Ensure output is a valid integer and never below N_min
      N_cell = max(N_real, N_min)
	  
	end function cell_budget
	
	function xy_res (x_min,x_max,y_max) result(res)
	  !! ~~~~ description ~~~~
	  !! This function will compute the resolution in x-y space
	  
	  use data_module
	   
	  implicit none
	  real(i_kind), intent(in) :: x_min,x_max,y_max								![-]					|droplet matrix
	  real(i_kind) :: res														![m]					|resolution
	  integer, dimension(31) :: ex												![-]					|exponent array
	  real(i_kind), dimension(31) :: dx, xmin, xmax, ymax, nx, ny, ncell
	  real(i_kind) :: n_allowed, x_span
	  
	  ex = (/-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20/)
	  
	  ! resolution candidates
      dx   = 1._i_kind / (2._i_kind ** real(ex, i_kind))
      
      ! x and y bounds
      xmin = floor(x_min / dx) * dx
      xmax = ceiling(x_max / dx) * dx
      ymax = ceiling(y_max / dx) * dx
	  
	  ! number of cells in x and y
      nx    = (xmax - xmin) / dx
      ny    = (2._i_kind * ymax) / dx
      
      ! total cell count
      ncell = nx * ny
	  
	  
	  !compute x_span
	  x_span = x_max - x_min
	  if(x_span < 1._i_kind) then
	    res = minval(pack(dx, mask = ncell < 9000._i_kind))
	  else
	    !calculate allowed cells
		n_allowed = cell_budget(x_span,30000._i_kind,180000._i_kind,700._i_kind,500._i_kind,1.3_i_kind)
      
	    !finalize resolution
	    res = minval(pack(dx, mask = ncell < n_allowed))
	  end if
	  
	end function xy_res
	
	
    function time_res(t1, t2) result(res)
      !! ~~~~ description ~~~~
      !! Compute time resolution from t1–t2 using control_dat%max_n_time
      
      use data_module
      
      implicit none
      real(i_kind), intent(in) :: t1, t2    ! [s] time bounds
      real(i_kind) :: res                   ! [s] time step width
      real(i_kind) :: N_t					! [-] number of time steps
      real(i_kind) :: tspan 				! [s] span min to max
	  
	  !calculate tspan
	  tspan = t2 - t1
	  
	  !calculate number of time steps
	  N_t = max(51._i_kind, + 23._i_kind*tspan**0.5_i_kind)
	  
	  !calculate resolution
	  res = tspan / N_t
    
    end function time_res
	
	
	function interpolate(arr_a, arr_b,a) result(b)
	  !! ~~~~ description ~~~~
	  !! This function will compute the value of b based on all values of a and b and the value of a
	  
	  implicit none
	  real(i_kind), dimension(:), intent(in) :: arr_a							![-]					|incomming array of variable a
	  real(i_kind), dimension(:), intent(in) :: arr_b							![-]					|incomming array of variable b
	  real(i_kind), intent(in) :: a												![-]					|incomming value of variable a
	  real(i_kind) :: b															![-]					|resulting value of variable b
	  real(i_kind) :: a_1, a_2													![-]					|temporal values
	  real(i_kind) :: b_1, b_2													![-]					|temporal values
	  integer :: loc															![-]					|location in array
	  
	  loc = maxloc(arr_a, dim = 1, mask = (arr_a<a))
	  
	  a_1 = arr_a(loc)
	  a_2 = arr_a(loc+1)
	  b_1 = arr_b(loc)
	  b_2 = arr_b(loc+1)
	  
	  b = b_1 + (a-a_1)*((b_2-b_1)/(a_2-a_1))
	
	end function interpolate
	
	
	function interpolate_idx(arr_a,arr_b,idx_1,idx_2,a) result(b)
	  !! ~~~~ description ~~~~
	  !! This function will compute the value of b based on all values of a and b, index 1&2, and the value of a
	  
	  implicit none
	  real(i_kind), dimension(:), intent(in) :: arr_a							![-]					|incomming array of variable a
	  real(i_kind), dimension(:), intent(in) :: arr_b							![-]					|incomming array of variable b
	  integer, intent(in) :: idx_1												![-]					|incomming index 1
	  integer, intent(in) :: idx_2												![-]					|incomming index 2
	  real(i_kind), intent(in) :: a												![-]					|incomming value of variable a
	  real(i_kind) :: b															![-]					|resulting value of variable b
	  real(i_kind) :: a_1, a_2													![-]					|temporal values
	  real(i_kind) :: b_1, b_2													![-]					|temporal values
	  
	  a_1 = arr_a(idx_1)
	  a_2 = arr_a(idx_2)
	  b_1 = arr_b(idx_1)
	  b_2 = arr_b(idx_2)
	  
	  b = b_1 + (a-a_1)*((b_2-b_1)/(a_2-a_1))
	
	end function interpolate_idx
	
	
    function interpolate_2d_geo(arr_t, arr_z, arr_mtx, t_idx_1, t_idx_2, z_idx_1_t1, z_idx_2_t1, z_idx_1_t2, z_idx_2_t2, t_val, z_val) result(b)
      !! ~~~~ description ~~~~
      !! Geometric bilinear interpolation on puff slice data arrays:
      !! - uses precomputed lower/upper indices in time (t_idx_1, t_idx_2)
      !!   and vertical space (z_idx_1, z_idx_2)
      !! - multiplicative blending in log-space for numerical robustness
    
      implicit none
      real(i_kind), dimension(:),   intent(in) :: arr_t          ! time array
      real(i_kind), dimension(:,:), intent(in) :: arr_z          ! z positions [time, N_slices]
      real(i_kind), dimension(:,:), intent(in) :: arr_mtx        ! quantity to interpolate [time, N_slices]
      integer, intent(in) :: t_idx_1, t_idx_2                    ! lower/upper time indices
      integer, intent(in) :: z_idx_1_t1, z_idx_2_t1              ! lower/upper vertical indices of t1
      integer, intent(in) :: z_idx_1_t2, z_idx_2_t2              ! lower/upper vertical indices of t2
      real(i_kind), intent(in) :: t_val, z_val                   ! target time and z position
      real(i_kind) :: b                                          ! resulting interpolated value
      real(i_kind) :: t1, t2, f_t
      real(i_kind) :: z1_t1, z2_t1, z1_t2, z2_t2
      real(i_kind) :: f_z_t1, f_z_t2
      real(i_kind) :: v11, v12, v21, v22
      real(i_kind) :: log_v11, log_v12, log_v21, log_v22
      real(i_kind) :: log_v_t1, log_v_t2
    
      ! --- time interpolation fraction ---
      t1  = arr_t(t_idx_1)
      t2  = arr_t(t_idx_2)
      f_t = (t_val - t1) / (t2 - t1)
      
      ! --- vertical geometry at t_idx_1 ---
      z1_t1 = arr_z(t_idx_1, z_idx_1_t1)
      z2_t1 = arr_z(t_idx_1, z_idx_2_t1)
      f_z_t1 = (z_val - z1_t1) / (z2_t1 - z1_t1)
      
      ! --- vertical geometry at t_idx_2 ---
      z1_t2 = arr_z(t_idx_2, z_idx_1_t2)
      z2_t2 = arr_z(t_idx_2, z_idx_2_t2)
      f_z_t2 = (z_val - z1_t2) / (z2_t2 - z1_t2)
      
      ! --- clamp fractions defensively ---
      f_z_t1 = max(0._i_kind, min(1._i_kind, f_z_t1))
      f_z_t2 = max(0._i_kind, min(1._i_kind, f_z_t2))
      f_t    = max(0._i_kind, min(1._i_kind, f_t))
      
      ! --- corner values ---
      v11 = max(arr_mtx(t_idx_1, z_idx_1_t1), 1.0e-12_i_kind)
      v12 = max(arr_mtx(t_idx_1, z_idx_2_t1), 1.0e-12_i_kind)
      v21 = max(arr_mtx(t_idx_2, z_idx_1_t2), 1.0e-12_i_kind)
      v22 = max(arr_mtx(t_idx_2, z_idx_2_t2), 1.0e-12_i_kind)
      
      ! --- log transform ---
      log_v11 = log(v11)
      log_v12 = log(v12)
      log_v21 = log(v21)
      log_v22 = log(v22)
      
      ! --- vertical interpolation in log space ---
      log_v_t1 = (1._i_kind - f_z_t1) * log_v11 + f_z_t1 * log_v12
      log_v_t2 = (1._i_kind - f_z_t2) * log_v21 + f_z_t2 * log_v22
      
      ! --- time interpolation in log space ---
      b = exp( (1._i_kind - f_t) * log_v_t1 + f_t * log_v_t2 )
    
    end function interpolate_2d_geo

	
	
	function expand_grid (a,b) result(c)
	  !! ~~~~ description ~~~~
	  !! This function will compute all combinatins of all values of the incomming arrays
	  
	  implicit none
	  real(i_kind), dimension(:), intent(in) :: a								![-]					|array a
	  real(i_kind), dimension(:), intent(in) :: b								![-]					|array b
	  real(i_kind), dimension(:,:), allocatable :: c							![-]					|arranged matrix
	  integer :: n_a															![-]					|length a
	  integer :: n_b															![-]					|length b
	  integer :: n_c 															![-]					|length c
	  integer :: n																![-]					|counter
	  integer :: i																![-]					|counter
	  integer :: j																![-]					|counter
	  
	  n_a = size(a)
	  n_b = size(b)
	  n_c = n_a*n_b
	  allocate(c(n_c,2))
	  
	  i = 1
	  j = 1
	  
	  do n=1, n_c
	  c(n,:) = (/a(i),b(j)/)
	  if(j==n_b)then
	  i = i+1
	  j = 1
	  else
	  j = j+1
	  end if	  
	  end do
	  	
	end function expand_grid
	
	
	function settling_velocity (D,V_in,rho_air,mu_air,rho_d) result (Vs)
	  !! ~~~~ description ~~~~
	  !! This function will compute the settling velocity iteratively 
	  
	  use constants_module
	  use data_module
	  
	  implicit none
	  real(i_kind), intent(in) :: D												![m]						|drop diameter
	  real(i_kind), intent(in) :: V_in											![m/s]						|initial guess of speed
	  real(i_kind), intent(in) :: rho_air										![kg/m³]					|density of air
	  real(i_kind), intent(in) :: mu_air										![kg/ms]					|dynamic viscosity of air
	  real(i_kind), intent(in) :: rho_d											![kg/m³]					|density of drop
	  real(i_kind) :: Vs														![m/s]						|settling velocity
	  real(i_kind), parameter :: a = 24._i_kind									![-]						|empirical parameter
	  real(i_kind), parameter :: b = 0.32_i_kind								![-]						|empirical parameter
	  real(i_kind), parameter :: c = 0.52_i_kind								![-]						|empirical parameter
	  real(i_kind) :: Re														![-]						|Reynolds number
	  real(i_kind) :: Cd														![-]						|drag coefficient
	  
	  
	  !!one itteration in the direction of the empirical settling velocity
	  !Reynolds number
	  Re = (D*V_in*rho_air)/mu_air
	  !drag coefficient
	  Cd = ((a/Re)**c+(b**c))**(1._i_kind/c)
	  !velocity
	  Vs = ((4._i_kind*rho_d*g*D)/(3._i_kind*rho_air*Cd))**0.5_i_kind
	  
	end function settling_velocity
	
	
	function wind_profile(z) result(U)
      !! ~~~~ description ~~~~
	  !! This function will compute the windspeed for a given height
	  
	  use data_module
	  
	  implicit none
	  real(i_kind) :: z			 												![m]					|vertical position
	  real(i_kind) :: U															![m/s]					|wind speed
	  
	  !wind speed
	    if(z >= env_dat%Hv)then
	      U = (env_dat%U_fric/Kc)*log((z-env_dat%d)/env_dat%z0)
	    else
	      U = env_dat%Uh*exp((env_dat%LAI/2._i_kind)*((z/env_dat%Hv)-1._i_kind))
	    end if
		
		! Cap wind speed to realistic upper limit
        U = min(U, env_dat%U_wind)
        U = max(U, 0._i_kind)
	  
	  end function wind_profile
	  
	
	function sutherland_formula(T_K) result(mu)
      !! ~~~~ description ~~~~
	  !! This function will compute the dynamic viscosity of air at a given temperature.
	  
	  implicit none
	  real(i_kind) :: T_0S = 273._i_kind 										![K]					|reference temperature of air
	  real(i_kind) :: mu_0S = 1.716E-5_i_kind									![kg/ms]				|reverence viscosity of air
	  real(i_kind) :: S_c = 111._i_kind											![K]					|sutherlnd constant of air
	  real(i_kind), intent(in) :: T_K											![K]					|air temperature
	  real(i_kind) :: mu														![kg/ms]				|dynamic viscosity of air
	  
	  mu = mu_0S*((T_0S+S_c)/(T_K+S_c))*((T_K/T_0S)**(3._i_kind/2._i_kind))
	  
	end function sutherland_formula
  
  
    function saturated_vapor_pressure(T_C) result(P_sat)
	  !! ~~~~ description ~~~~
	  !! This function will compute the saturated vapor pressure of water at a given temperature.
	  
	  implicit none
	  real(i_kind), intent(in) :: T_C											![°C]					|air temperature
	  real(i_kind) :: P_sat														![Pa]					|saturated vapor pressure
	  
	  P_sat = 610.78_i_kind*10._i_kind**((7.5_i_kind*T_C)/(T_C+237.3_i_kind))
	
	end function saturated_vapor_pressure
	   
	   
    function wet_bulb_temperature(T_C, Rh_P) result(T_wb)
	  !! ~~~~ description ~~~~
	  !! This function will compute the wet-bulb temperature at a given temperature.
	  
	  implicit none
	  real(i_kind), intent(in) :: T_C											![°C]					|air temperature
	  real(i_kind), intent(in) :: Rh_P											![%]					|relative humidity
	  real(i_kind) :: T_wb														![°C]					|wet-bulb temperature
	  
	  T_wb = T_C * atan(0.151977_i_kind*(Rh_P + 8.313659_i_kind)**0.5_i_kind)+atan(T_C + Rh_P)-atan(Rh_P-1.676331_i_kind)+0.00391838_i_kind*(Rh_P)**(3._i_kind/2._i_kind)*atan(0.023101_i_kind*Rh_P)-4.686035_i_kind
	
	end function wet_bulb_temperature
	  
	  
	function diffusion_coefficient_water(T_C) result(D_w)
	  !! ~~~~ description ~~~~
	  !! This function will compute the diffusion coefficient of water in air at a given temperature.
	  
	  implicit none
	  real(i_kind), intent(in) :: T_C											![°C]					|air temperature
	  real(i_kind) :: D_w														![m²/s]					|diffusion coefficient
	  
	  D_w = 21.2E-6_i_kind*(1._i_kind+0.0071_i_kind*T_C)
	  
	end function diffusion_coefficient_water
	  
	  
	function water_density(T_C) result(rho_h2o)
	  !! ~~~~ description ~~~~
	  !! This function will compute the density of water at a given temperature.
	  
	  implicit none
	  real(i_kind), intent(in) :: T_C											![°C]					|air temperature
	  real(i_kind) :: rho_h2o													![kg/m³]				|density of water
	  
	  rho_h2o = 999.85308_i_kind+6.32693E-2_i_kind*T_C-8.523829E-3_i_kind*T_C**2._i_kind+6.943248E-5_i_kind*T_C**3._i_kind-3.821216E-7_i_kind*T_C**4._i_kind
	
	end function water_density  
	
	
	function absolute_humidity(rh,p_sat,T) result(ah)
	  !! ~~~~ description ~~~~
	  !! This function will compute the absolute humidity
	  
	  implicit none
	  real(i_kind), intent(in) :: rh											![-]					|relative humidity
	  real(i_kind), intent(in) :: p_sat											![Pa]					|saturation vapour pressure
	  real(i_kind), intent(in) :: T												![K]					|air temperature
	  real(i_kind) :: ah														![kg/m³]				|absolute humidity
	  
	  ah = (rh*p_sat)/(R_h2o * T)
	  
	end function absolute_humidity  
	
	
	function relative_humidity(p_sat,m,T,V) result(rh)
	  !! ~~~~ description ~~~~
	  !! This function will compute the relative humidity
	  
	  implicit none
	  real(i_kind), intent(in) :: m												![kg]					|mass of water vapour
	  real(i_kind), intent(in) :: p_sat											![Pa]					|saturation vapour pressure
	  real(i_kind), intent(in) :: T												![K]					|air temperature
	  real(i_kind), intent(in) :: V												![m³]					|volume
	  real(i_kind) :: rh														![kg/m³]				|absolute humidity
	  
	  rh = (m*R_h2o*T)/(V*p_sat)
	
	end function relative_humidity
	
	
	function air_density(Psat,Rh,P,T) result(rho_air)
	  !! ~~~~ description ~~~~
	  !! This function will compute the air density
	  
	  implicit none
	  real(i_kind), intent(in) :: Rh											![-]					|relative humidity
	  real(i_kind), intent(in) :: Psat											![Pa]					|saturation vapour pressure
	  real(i_kind), intent(in) :: P												![PA]					|air pressure
	  real(i_kind), intent(in) :: T												![K]					|temperature
	  real(i_kind) :: rho_air													![kg/m³]				|air density
	  
	  rho_air = ((Psat * Rh * mm_h2o)+((P - (Psat * Rh))* mm_da))/(R * T)
	  
	end function air_density
	
	
	function latent_heat_of_vaporization(T) result(Hs)
	  !! ~~~~ description ~~~~
	  !! This function will compute the air density
	  
	  implicit none
	  real(i_kind), intent(in) :: T												![K]					|temperature
	  real(i_kind) :: Hs														![J/kg]					|
	  
	  Hs = (50.09_i_kind-0.9298_i_kind*T/1000._i_kind-65.19_i_kind*(T/1000._i_kind)**2._i_kind)*1000._i_kind/mm_h2o
	  
	end function latent_heat_of_vaporization

	
	function cdf_normal_sigma(sigma,mean,x) result(f)
	  !! ~~~~ description ~~~~
	  !! This function will estimate how much of an distribution is located above agiven z height
	  
	  implicit none
	  real(i_kind), intent(in) :: sigma											![m²/s]					|eddy difussion parameter
	  real(i_kind), intent(in) :: x												![m]					|value
	  real(i_kind), intent(in) :: mean											![m]					|mean value
	  real(i_kind) :: f															![-]					|fraction above cut height
	  
	  f = 0.5_i_kind*(1._i_kind+erf((x-mean)/((sigma)*((2._i_kind)**0.5_i_kind))))
	  
  	end function cdf_normal_sigma
	
	
	function souton_sigma_h(sig0,t,x) result(sigma)
	  !! ~~~~ description ~~~~
	  !! This function will calculate the horizontal dispersion coefficient 
	  
	  implicit none
	  real(i_kind), intent(in) :: sig0											![m²/s]					|eddy difussion parameter
	  real(i_kind), intent(in) :: x												![m]					|value
	  real(i_kind), intent(in) :: t												![m]					|mean value
	  real(i_kind) :: sigma														![-]					|fraction above cut height
	  
	  sigma = ((sig0*x) / (1._i_kind+(0.9_i_kind*((0.001_i_kind*t)**0.5_i_kind))))
	  
  	end function souton_sigma_h  
	
	
	function souton_sigma_v(sig0,t) result(sigma)
	  !! ~~~~ description ~~~~
	  !! This function will calculate the vertical dispersion coefficient
	  
	  implicit none
	  real(i_kind), intent(in) :: sig0											![m²/s]					|eddy difussion parameter
	  real(i_kind), intent(in) :: t												![m]					|value
	  real(i_kind) :: sigma														![-]					|fraction above cut height
	  
	  sigma = ((sig0*t) / (1._i_kind+(0.9_i_kind*((0.02_i_kind*t)**0.5_i_kind))))
	  
  	end function souton_sigma_v
	
	function souton_v_Kz(sig0, t) result(Kz)
      !! ~~~~ description ~~~~
	  !! This function will calculate the derivate of the vertical dispersion coefficient 
      implicit none
	  
      real(i_kind), intent(in) :: sig0   ! vertical eddy diffusion scale parameter
      real(i_kind), intent(in) :: t      ! time [s]
	  
      real(i_kind) :: Kz
      real(i_kind) :: d, sqrt_t, denom, factor
	  
      ! constant from Sutton parameterization
      d = 0.9_i_kind * sqrt(0.02_i_kind)
	  
      sqrt_t = sqrt(t)
      denom  = (1.0_i_kind + d * sqrt_t)**3
      factor = (1.0_i_kind + 1.5_i_kind * d * sqrt_t)
	  
      Kz = (sig0**2 * t / denom) * factor
    end function souton_v_Kz
	
	
	function unique(vec) result(vec_unique)  
      ! Return only the unique values from vec.
      
      implicit none
      real(i_kind), dimension(:), intent(in) :: vec								![-]					|array
      real(i_kind), dimension(:), allocatable :: vec_unique						![-]					|array
      integer :: i																![-]					|counter
	  integer :: num															![-]					|counter
      logical, dimension(size(vec)) :: mask										![-]					|mask
      
      mask = .false.
      do i=1,size(vec)
        !count the number of occurrences of this element:  
        num = count(vec(i) == vec)
        if (num == 1) then  
          !there is only one, flag it:  
          mask(i) = .true.  
        else  
          !flag this value only if it hasn't already been flagged:  
          if (.not. any(vec(i)==vec .and. mask)) mask(i) = .true.  
        end if
      end do
      
      !return only flagged elements:  
      allocate(vec_unique(count(mask)))  
      vec_unique = pack(vec, mask)
      
    end function unique
	
	
	function position_B_in_A (B,A) result(C)
	  !returns the positions of B values in A array
	  
	  implicit none
	  real(i_kind), dimension(:), intent(in) :: A								![-]					|array
	  real(i_kind), dimension(:), intent(in) :: B								![-]					|array
	  integer, dimension(:), allocatable :: C									![-]					|array
	  integer :: i																![-]					|counter
	  integer :: j																![-]					|counter
	  logical, dimension(:), allocatable :: mask								![-]					|maks
	  logical, dimension(:), allocatable :: tmp									![-]					|temporal array
	  
	  !allocate C, mask and tmp
	  allocate(mask(size(A)))
	  allocate(tmp(size(B)))
	  
	  !loop over all elements in A
	  do i=1,size(A)
	    tmp = .false.
	    !loop over all elements in B
		do j=1,size(B)
	      tmp(j) = abs(A(i)-B(j)) < 1.e-6_i_kind
	    end do
		mask(i) = any(tmp)
	  end do
	  C = pack([(i,i=1,size(A))],mask = mask)
	  
	  !deallocate
	  if(allocated(mask)) deallocate(mask)
	  if(allocated(tmp)) deallocate(tmp)
	  
	end function position_B_in_A
	
	
	function cbind (mat1,mat2) result (mat3)
    implicit none
    real(i_kind), dimension(:,:), intent(in) :: mat1							![-]					|matrix
	real(i_kind), dimension(:,:), intent(in) :: mat2							![-]					|matrix
	real(i_kind), dimension(:,:), allocatable :: mat3							![-]					|matrix
	integer :: nc1																![-]					|number of columns
	integer :: nr1																![-]					|number of rows
	integer :: nc2																![-]					|number of columns
	integer :: i																![-]					|counter
	
	!read dimesnions of incomming matrixes
	nc1 = size(mat1(1,:))
	nr1 = size(mat1(:,1))
	nc2 = size(mat2(1,:))
	
	!allocate outgoing matrix
	allocate(mat3(nr1,(nc1+nc2)))
	
	!add mat1
	mat3(:,[(i,i=1,nc1)]) = mat1
	
	!add mat2
	mat3(:,[(i,i=(nc1+1),(nc1+nc2))]) = mat2
	
  end function cbind
  
  
    function rbind (mat1,mat2) result (mat3)
      implicit none
      real(i_kind), dimension(:,:), intent(in) :: mat1							![-]					|matrix
  	  real(i_kind), dimension(:,:), intent(in) :: mat2							![-]					|matrix
  	  real(i_kind), dimension(:,:), allocatable :: mat3							![-]					|matrix
  	  integer :: nc1															![-]					|number of columns
  	  integer :: nr1															![-]					|number of rows
  	  integer :: nr2															![-]					|number of rows
  	  integer :: i																![-]					|counter
  	  
  	  !read dimesnions of incomming matrixes
  	  nc1 = size(mat1(1,:))
  	  nr1 = size(mat1(:,1))
  	  nr2 = size(mat2(:,1))
  	  
  	  !allocate outgoing matrix
  	  allocate(mat3((nr1+nr2),nc1))
  	  
  	  !add mat1
  	  mat3([(i,i=1,nr1)],:) = mat1
  	  
  	  !add mat2
  	  mat3([(i,i=(nr1+1),(nr1+nr2))],:) = mat2
  	  
    end function rbind
	
  end module functions_module