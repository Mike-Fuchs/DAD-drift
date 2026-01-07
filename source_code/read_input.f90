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

subroutine read_input
  !! ~~~~ description ~~~~
  !! this routine will read the input files
  
  use constants_module
  use file_module
  use data_module
  use ascii_grid_module
      
  implicit none
  integer :: i 														![-]					|counter
  real(i_kind), dimension(11) :: env_in								![-]					|temporal input array
  real(i_kind), dimension(14) :: app_in								![-]					|temporal input array
  character(len=60),dimension(2) :: con_in_cha1						![-]					|temporal input array
  character(len=60),dimension(:), allocatable :: con_in_cha2		![-]					|temporal input array
  real(i_kind), dimension(4) :: con_in								![-]					|temporal input array
  type(header) :: tmp_header										![-]					|temporal input header
  real(i_kind), dimension(:), allocatable :: cellsize				![m]					|raster cellsize
  integer, dimension(:), allocatable :: spec_vec					![-]					|vector of the non zero specrums
  integer :: nlines													![-]					|counter
    
  !!read controle input
  open(unit=104, file=in_files%controle_input, status='old', action='read')
  read(104,*) 
  read(104,*) con_in
  read(104,*) con_in_cha1
  if (con_in(1) == 0) then
    allocate(con_in_cha2(1+int(con_in(4))))
    read(104,*) con_in_cha2
  end if
  close(104)
  control_dat%code1 = int(con_in(1))
  control_dat%z_0 = con_in(2)
  control_dat%max_dist = con_in(3)
  control_dat%field_count = int(con_in(4))
  in_files%spec_input = 'input/' // trim(con_in_cha1(2))
  if (control_dat%code1 == 0) then
    if(control_dat%field_count > 0)then
      allocate(in_files%landscape_file_name(control_dat%field_count))
	  allocate(cellsize(control_dat%field_count))
      do i=1,control_dat%field_count
        in_files%landscape_file_name(i) = 'input/' // trim(con_in_cha2(i+1))
      end do
    end if
  end if
  
  !!reading environmental input file
  open(unit=100, file=in_files%env_input, status='old', action='read')
  read(100,*)
  read(100,*) 
  read(100,*) env_in
  close(100)
  !write data to env_dat
  env_dat%T_C = env_in(1)
  env_dat%Rh = env_in(2)
  env_dat%U_wind = env_in(3)
  env_dat%z_wind = env_in(4)
  env_dat%D_wind = env_in(5)
  env_dat%P = env_in(6)*1000._i_kind
  env_dat%sig_h = env_in(7)
  env_dat%sig_v = env_in(8)
  env_dat%z0 = env_in(9)
  env_dat%Hv = env_in(10)
  env_dat%LAI = env_in(11)
  
  !!read droplet distribution
  nlines = 0
  open(unit=101, file=in_files%spec_input, status='old', action='read')
  read(101,*)
  read(101,*)
  do
      read (101,*, END=10)
      nlines = nlines + 1
  end do
  10 close (101)
  allocate(drop_spec%ds(nlines))
  allocate(drop_spec%cf(nlines))
  allocate(drop_spec%f(nlines))
  open(unit=101, file=in_files%spec_input, status='old', action='read')
  read(101,*)
  read(101,*)
  do i=1,nlines
    read(101,*) drop_spec%ds(i), drop_spec%cf(i)
  end do
  close(101)
  !calculate fraction per droplet class
  drop_spec%f(1) = drop_spec%cf(1)
  do i=2,nlines
    drop_spec%f(i) = drop_spec%cf(i)-drop_spec%cf(i-1)
  end do
  
  !!droplet spectrum limitation
  !determine non zero droplet classes
  spec_vec = pack([(i,i=1,size(drop_spec%ds))],drop_spec%f > 0._i_kind)
  
  !limit droplet spectrum to non zero droplet classes
  drop_spec%ds = drop_spec%ds(spec_vec)
  drop_spec%cf = drop_spec%cf(spec_vec)
  drop_spec%f = drop_spec%f(spec_vec)
  
  !initialise grid list
  allocate(grid_list(size(drop_spec%ds)))
  
  !!read application input
  open(unit=102, file=in_files%app_input, status='old', action='read')
  read(102,*)
  read(102,*)
  read(102,*) app_in
  close(102)
  !write data to drop_dat and app_dat
  app_dat%v_trac = app_in(1)
  app_dat%b_width = app_in(2)
  app_dat%b_height = app_in(3)
  app_dat%nozzle_angle = app_in(4)
  app_dat%app_pressure = app_in(5)*1000._i_kind
  app_dat%app_rate_mh = app_in(6)
  app_dat%app_rate_mha = app_in(7)
  app_dat%app_rate = app_in(8)
  app_dat%c_solution = app_in(9)
  app_dat%rho_AI = app_in(10)
  app_dat%mm_AI = app_in(11)
  app_dat%P_AI = app_in(12)
  app_dat%n_pass = int(app_in(13))
  app_dat%f_length = real(ceiling(app_in(14)),i_kind)
  
  !!read landscape input raster if code is 0
  if (control_dat%code1 == 0) then
    !read headers
	do i=1,control_dat%field_count
	  tmp_header = read_ascii_header(in_files%landscape_file_name(i))
	  cellsize(i) = tmp_header%cellsize
	end do
	
	!define cell size
	control_dat%cellsize = cellsize(1)
	
	!check if all cellsizes are equal
	if(.not. all(control_dat%cellsize == cellsize))then
	  stop 'landscape raster cellsizes are unequal between different files'
	end if
	
	!check if landscape raster for resolution to be larger 1m
	if(control_dat%cellsize < 1)then
	  stop 'landscape raster cellsize must be equal or larger than 1m'
	end if
	
	!check if landscape raster for resolution to be a multiple of 1
	if(mod(real(control_dat%cellsize),1.) /= 0.)then
	  stop 'landscape raster cellsize must be a hole multiple of 1m'
	end if
  end if
  
  !deallocate arrays
  if(allocated(spec_vec)) deallocate(spec_vec)
  if(allocated(con_in_cha2)) deallocate(con_in_cha2)
  if(allocated(cellsize)) deallocate(cellsize)
  
  print *, 'input reading successful' 
  
  end subroutine read_input