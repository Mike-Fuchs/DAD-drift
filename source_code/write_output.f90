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

subroutine write_output
  !! ~~~~ description ~~~~
  !! this routine will write the output
    
  use file_module
  use data_module
  use header_module
  use ascii_grid_module
  
  implicit none
  integer :: j															![-]					|counter
  integer :: i															![-]					|counter
  integer :: nrow_dev													![-]					|number of rows
  character(len=20) :: tmp_head											![-]					|temporal header
  
  print *, 'writing output'
  
  if(.not. control_dat%dev_mode)then
	!when controle code1 is 1 write drift curve to output file
	if(control_dat%code1 == 1) then
	  !write drift curve
	  open(unit=113, file=out_files%drift_curve_output, status='replace', action='write')
	  write(113,124) head_dc
	  write(113,124) units_dc
	  do j=1,size(drift_curve_matrix(:,1))
	    write(113,125) drift_curve_matrix(j,:)
	  end do
	  close(113)
	end if
	
	!when controle code1 is 0 write landscape drift raster
	if(control_dat%code1 == 0) then
	  call write_ascii_grid(landscape_drift_raster, out_files%landscape_drift)
	end if
	
  else
    !write drift pattern
	call write_ascii_grid(drift_pattern, out_files%drift_pattern)
	
    !write local environment file
	nrow_dev = size(loc_env%time)
	open(unit=150, file="output/loc_env.txt", status='replace', action='write')
	write(150,126) 'time', 'volume', 'd_volume', 'm_air', 'd_m_air', 'm_h2o', 'd_m_h2o', 'd_E_air', 'd_E_h2o', 'E_tot_t1', 'm_h2o_drop', 'd_m_h2o_evap', &
				   'd_E_h2o_evap', 'T_sys_K', 'T_sys_C', 'Psat', 'Rh', 'mu_air', 'rho_air', 'T_wb_C', 'T_f_C', 'Psat_f', 'mu_air_f', 'rho_air_f', 'rho_h2o_f', &
				   'D_h2o_f', 'Sc_h2o_f', 'dt'
	do i=1,nrow_dev,10
	  write(150,127) loc_env%time(i), loc_env%volume(i), loc_env%d_volume(i), loc_env%m_air(i), loc_env%d_m_air(i), loc_env%m_h2o(i), loc_env%d_m_h2o(i), &
					 loc_env%d_E_air(i), loc_env%d_E_h2o(i), loc_env%E_tot_t1(i), loc_env%m_h2o_drop(i), loc_env%d_m_h2o_evap(i), loc_env%d_E_h2o_evap(i), loc_env%T_sys_K(i), &
					 loc_env%T_sys_C(i), loc_env%Psat(i), loc_env%Rh(i), loc_env%mu_air(i), loc_env%rho_air(i), loc_env%T_wb_C(i), loc_env%T_f_C(i), &
					 loc_env%Psat_f(i), loc_env%mu_air_f(i), loc_env%rho_air_f(i), loc_env%rho_h2o_f(i), loc_env%D_h2o_f(i), loc_env%Sc_h2o_f(i), loc_env%dt(i)
	end do
	close(150)
	
	!write droplet model outputs
	do j=1, size(drop_mod)
	  nrow_dev =  size(drop_mod(j)%time)
	  write (tmp_head,128) j
	  open(unit=151, file='output/drop_mod_'//trim(tmp_head)//'.txt', status='replace', action='write')
	  write(151,129) 'time', 'diameter', 'volume', 'm_h2o', 'm_AI', 'fract_AI', 'x_h2o', 'x', 'z', 'V_x', 'V_z', 'Re_s', 'f', 'tau_m', 'Re_f', 'Sh_h2o_f', 'F_h2o', &
					 'F_AI', 'dV_x_dt', 'dV_z_dt', 'dm_h2o_dt', 'dm_AI_dt', 'dm_h2o', 'dm_AI', 'dm_h2o_a', 'V_s', 'edv', 'dag', 'dag_res', 'dag_res_chg', 'dt', 'sigma_v' 
	  do i=1,nrow_dev,10
	    write(151,130) drop_mod(j)%time(i), drop_mod(j)%diameter(i), drop_mod(j)%volume(i), drop_mod(j)%m_h2o(i), drop_mod(j)%m_AI(i), drop_mod(j)%fract_AI(i), drop_mod(j)%x_h2o(i), &
					   drop_mod(j)%x(i), drop_mod(j)%z(i), drop_mod(j)%V_x(i), drop_mod(j)%V_z(i), drop_mod(j)%Re_s(i), drop_mod(j)%f(i), drop_mod(j)%tau_m(i), &
					   drop_mod(j)%Re_f(i), drop_mod(j)%Sh_h2o_f(i), drop_mod(j)%F_h2o(i), drop_mod(j)%F_AI(i), drop_mod(j)%dV_x_dt(i), drop_mod(j)%dV_z_dt(i), drop_mod(j)%dm_h2o_dt(i), &
					   drop_mod(j)%dm_AI_dt(i), drop_mod(j)%dm_h2o(i), drop_mod(j)%dm_AI(i), drop_mod(j)%dm_h2o_a(i), drop_mod(j)%V_s(i), drop_mod(j)%edv(i), drop_mod(j)%dag(i), &
					   drop_mod(j)%dag_res(i), drop_mod(j)%dag_res_chg(i), drop_mod(j)%dt(i), drop_mod(j)%sigma_v(i)
	  end do
	  close(151)
	end do
	
	!when controle code1 is 1 write drift curve to output file
	if(control_dat%code1 == 1) then
	  !write drift curve
	  open(unit=113, file=out_files%drift_curve_output, status='replace', action='write')
	  write(113,124) head_dc
	  write(113,124) units_dc
	  do j=1,size(drift_curve_matrix(:,1))
	    write(113,125) drift_curve_matrix(j,:)
	  end do
	  close(113)
	end if
	
	!when controle code1 is 0 write landscape drift raster
	if(control_dat%code1 == 0) then
	  call write_ascii_grid(landscape_drift_raster, out_files%landscape_drift)
	end if
  end if
    
  !define formats
  !drift_curve
  124 format (2A20)
  125 format (1F20.2,1ES20.10E3)
  !dev out
  126 format (28A20)
  127 format (28E20.10)
  128 format (1I3.3)
  129 format (35A20)
  130 format (35E20.10)
  
end subroutine write_output