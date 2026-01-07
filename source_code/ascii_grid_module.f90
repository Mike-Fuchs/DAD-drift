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

module ascii_grid_module
  !! ~~~~ description ~~~~
  !! This module defines the ascii grid data format, read, write, and data manipulation functions.
  use constants_module
  
  implicit none
  private
  public :: header, ascii_grid, drift_pattern, landscape_drift_raster, landscape_input_raster, grid_list, read_ascii_header, read_ascii_grid, write_ascii_grid, create_ascii_grid, llz_to_ascii_grid, &
			llz_to_ascii_grid_with_ll, ascii_grid_to_llz, add_ascii_grid, rescale_ascii_grid, rescale_ascii_grid_centered, add_ascii_grid_at, add_ascii_grid_at_fast, add_ascii_grid_fast, add_ascii_grid_landscape
  
  !define header
  type header
    integer :: ncols													![-]					|number of columns
	integer :: nrows													![-]					|number of rows
	real(i_kind) :: xllcorner											![m]					|x coordinate of lower left corner
	real(i_kind) :: yllcorner											![m]					|y coordinate of lower left corner
	real(i_kind) :: cellsize											![m]					|raster cell size
	integer :: nodata_value												![-]					|no data value
  end type header
  
  !define header names
  type header_n
    character(len = 13) :: ncols = 'ncols'
	character(len = 13) :: nrows = 'nrows'
	character(len = 13) :: xllcorner = 'xllcorner'
	character(len = 13) :: yllcorner = 'yllcorner'
	character(len = 13) :: cellsize = 'cellsize'
	character(len = 13) :: nodata_value = 'nodata_value'  
  end type header_n
  type(header_n) :: header_names
 
 !define ascii_grid_raster
  type :: ascii_grid
    type(header) :: header												![-]					|ratser header
	real(i_kind), dimension(:,:), allocatable :: data					![-]					|raster data matrix
  end type ascii_grid
  
  !define grid pattern raster
  type(ascii_grid) :: drift_pattern
  
  !define landscape drift raster
  type(ascii_grid) :: landscape_drift_raster 
  
  !define landscape raster
  type(ascii_grid) :: landscape_input_raster 
  
  !define grid list
  type(ascii_grid), dimension(:), allocatable :: grid_list
  
  contains
    function read_ascii_header (file_name) result(head)
	  !! ~~~~ description ~~~~
	  !! This function will read an ascii grid header
	  
	  implicit none	  
	  character(*), intent(in) :: file_name								![-]					|file name
	  type(header) :: head												![-]					|raster header
	  character(len=20) :: dump											![-]					|temporal dump
	  
	  open(unit=1000, file=file_name, status='old', action='read')
	  read(1000,*) dump, head%nrows
	  read(1000,*) dump, head%ncols
	  read(1000,*) dump, head%xllcorner
	  read(1000,*) dump, head%yllcorner
	  read(1000,*) dump, head%cellsize
	  read(1000,*) dump, head%nodata_value
	  close(1000)
	    
	end function read_ascii_header
	
	
	function read_ascii_grid (file_name) result(raster)
	  !! ~~~~ description ~~~~
	  !! This function will read an ascii grid raster
	  
	  implicit none	  
	  character(*), intent(in) :: file_name								![-]					|file name
	  type(ascii_grid) :: raster										![-]					|raster
	  integer :: i														![-]					|counter
	  character(len=20) :: dump											![-]					|temporal dump
	  
	  open(unit=1000, file=file_name, status='old', action='read')
	  read(1000,*) dump, raster%header%ncols
	  read(1000,*) dump, raster%header%nrows
	  read(1000,*) dump, raster%header%xllcorner
	  read(1000,*) dump, raster%header%yllcorner
	  read(1000,*) dump, raster%header%cellsize
	  read(1000,*) dump, raster%header%nodata_value
	  
	  !correcting xllcorner and yllcorner to match with resolution
	  raster%header%xllcorner = real(nint(raster%header%xllcorner - raster%header%cellsize/2._i_kind),i_kind) + raster%header%cellsize/2._i_kind
	  raster%header%yllcorner = real(nint(raster%header%yllcorner - raster%header%cellsize/2._i_kind),i_kind) + raster%header%cellsize/2._i_kind
	  
	  allocate(raster%data(raster%header%nrows,raster%header%ncols))
	  do i = raster%header%nrows, 1, -1
        read(1000,*) raster%data(i,:)
      end do
	  close(1000)
	    
	end function read_ascii_grid
	
	
	subroutine write_ascii_grid (raster, file_name)	  
	  !! ~~~~ description ~~~~
	  !! This function will write an ascii grid raster
	  
	  implicit none	  
	  character(*), intent(in) :: file_name								![-]					|file name
	  type(ascii_grid), intent(in) :: raster							![-]					|raster
	  integer :: i														![-]					|counter
	  character(len=20) :: dummy										![-]					|temporal dump
	   
	  open(unit=1001, file=trim(file_name), status='replace', action='write')
	  write(dummy,'(I20)') 			raster%header%ncols
	  write(1001,'(A13,A20)') 		header_names%ncols, adjustl(dummy)
	  write(dummy,'(I20)') 			raster%header%nrows
	  write(1001,'(A13,A20)') 		header_names%nrows, adjustl(dummy)
	  write(dummy,'(F20.10)') 		raster%header%xllcorner
	  write(1001,'(A13,A20)') 		header_names%xllcorner, adjustl(dummy)
	  write(dummy,'(F20.10)') 		raster%header%yllcorner
	  write(1001,'(A13,A20)') 		header_names%yllcorner, adjustl(dummy)
	  write(dummy,'(F20.5)') 		raster%header%cellsize
	  write(1001,'(A13,A20)')	 	header_names%cellsize, adjustl(dummy)
	  write(dummy,'(I20)') 			raster%header%nodata_value
	  write(1001,'(A13,A20)') 		header_names%nodata_value, adjustl(dummy)
	  
	  if(size(raster%data) > 0)then
	    do i = raster%header%nrows, 1, -1
          write(1001,*) raster%data(i,:)
        end do
	  end if
	  close(1001)
	  	    
	end subroutine write_ascii_grid
	
	
	function create_ascii_grid (ncols, nrows, xllcorner, yllcorner, cellsize) result(raster)
	  !! ~~~~ description ~~~~
	  !! This function will create an empty ascii grid raster
	  
	  implicit none
	  integer, intent(in) :: ncols										![-]					|number of columns
	  integer, intent(in) :: nrows										![-]					|number of rows
	  real(i_kind), intent(in) :: xllcorner								![m]					|x coordinate of lower left corner
	  real(i_kind), intent(in) :: yllcorner								![m]					|y coordinate of lower left corner
	  real(i_kind), intent(in) :: cellsize								![m]					|raster cell size
	  type(ascii_grid) :: raster										![-]					|raster
	  
	  !define raster header
	  raster%header%ncols = ncols
	  raster%header%nrows = nrows
	  raster%header%xllcorner = xllcorner
	  raster%header%yllcorner = yllcorner
	  raster%header%cellsize = cellsize
	  raster%header%nodata_value = -99
	  
	  !allocate raster data
	  allocate(raster%data(raster%header%nrows,raster%header%ncols))
	  raster%data = 0._i_kind
	  
	end function create_ascii_grid
	
	
	function llz_to_ascii_grid (lat_long_z, res) result(raster)
	  !! ~~~~ description ~~~~
	  !! This function will convert lat-long matrix into a ascii grid raster
	  
	  use functions_module
	  
      implicit none
      real(i_kind), dimension(:,:), intent(in) :: lat_long_z			![-]					|lat-long matrix
	  real(i_kind), intent(in) :: res									![m]					|cellsize
	  type(ascii_grid) :: raster										![-]					|raster
	  integer :: i														![-]					|counter
	  integer :: j														![-]					|counter
	  integer :: n														![-]					|counter
	  
	  !write head elements
	  raster%header%ncols = size(unique(lat_long_z(:,2)))
	  raster%header%nrows = size(unique(lat_long_z(:,1)))
	  raster%header%xllcorner = minval(lat_long_z(:,2)) - (res/2._i_kind)
	  raster%header%yllcorner = minval(lat_long_z(:,1)) - (res/2._i_kind)
	  raster%header%cellsize = res
	  raster%header%nodata_value = -99
	
	  !converte z to matrix
	  allocate(raster%data(raster%header%nrows,raster%header%ncols))
	  
	  !nested do
	  n = 1
	  do i=raster%header%nrows,1,-1
	    do j=1,raster%header%ncols
		  raster%data(i,j) = lat_long_z(n,3)
		  n = n + 1
		end do	    
	  end do
	  
	end function llz_to_ascii_grid
	
	
	function llz_to_ascii_grid_with_ll (lat_long_z, res, lat_coord, long_coord) result(raster)
	  !! ~~~~ description ~~~~
	  !! This function will convert lat-long matrix into a ascii grid raster, with lateral and longitudinal arrays as additional input
	  
	  use functions_module
	  
      implicit none
      real(i_kind), dimension(:,:), intent(in) :: lat_long_z			![-]					|lat-long matrix
	  real(i_kind), intent(in) :: res									![m]					|cellsize
	  real(i_kind), dimension(:), intent(in) :: lat_coord				![m]					|lateral coordinates
	  real(i_kind), dimension(:), intent(in) :: long_coord  			![m]					|longitudinal coordinates
	  type(ascii_grid) :: raster										![-]					|raster
	  integer :: i,j,n
	  
	  !write head elements
	  raster%header%ncols = size(long_coord)
	  raster%header%nrows = size(lat_coord)
	  raster%header%xllcorner = minval(long_coord) - res/2._i_kind
	  raster%header%yllcorner = minval(lat_coord) - res/2._i_kind
	  raster%header%cellsize = res
	  raster%header%nodata_value = -99
	
	  !converte z to matrix
	  allocate(raster%data(raster%header%nrows,raster%header%ncols))
	  
	  !nested do
	  n = 1
	  do i=raster%header%nrows,1,-1
	    do j=1,raster%header%ncols
		  raster%data(i,j) = lat_long_z(n,3)
		  n = n + 1
		end do	    
	  end do
	  
	end function llz_to_ascii_grid_with_ll
	
	
	function ascii_grid_to_llz (raster) result(lat_long_z)
	  !! ~~~~ description ~~~~
	  !! 
	  
	  use functions_module
	  
	  implicit none
	  type(ascii_grid), intent(in) :: raster							![-]					|raster
	  real(i_kind), dimension(:,:), allocatable :: lat_long_z			![-]					|lat-long matrix
	  real(i_kind), dimension(:), allocatable :: lat					![m]					|lateral coordinates
	  real(i_kind), dimension(:), allocatable :: long					![m]					|longitudinal coordinates
	  integer :: i														![-]					|counter
	  integer :: j														![-]					|counter
	  integer :: n														![-]					|counter
	  
	  !create lat and long array
	  lat = (real([(i,i=1,raster%header%nrows,1)],i_kind) * raster%header%cellsize) + raster%header%yllcorner - (raster%header%cellsize/2._i_kind)
	  long = (real([(i,i=1,raster%header%ncols,1)],i_kind) * raster%header%cellsize) + raster%header%xllcorner - (raster%header%cellsize/2._i_kind)
	  
	  !create lat_long_z matrix
	  allocate(lat_long_z(raster%header%nrows*raster%header%ncols,3))
	  lat_long_z(:,1:2) = expand_grid(lat,long)
	  
	  !convert matrix to z
	  i = 1
	  j = raster%header%nrows
	  do n=1,size(lat_long_z(:,1))
	    lat_long_z(n,3) = raster%data(j,i)
		if(i < raster%header%ncols)then
		  i = i + 1
		else
		  j = j - 1
		  i = 1
		end if	    
	  end do
	  
	end function ascii_grid_to_llz
	
	
	function add_ascii_grid (raster1, raster2) result(raster3)
	  !! ~~~~ description ~~~~
	  !! This function will add two ascii grid rasters
	  
	  use functions_module
	
	  implicit none
	  type(ascii_grid), intent(in) :: raster1							![-]					|raster
	  type(ascii_grid), intent(in) :: raster2							![-]					|raster
	  type(ascii_grid) :: raster3										![-]					|raster
	  real(i_kind), dimension(:), allocatable :: arr_cols1				![-]					|column array
	  real(i_kind), dimension(:), allocatable :: arr_rows1				![-]					|row array
	  real(i_kind), dimension(:), allocatable :: arr_cols2				![-]					|column array
	  real(i_kind), dimension(:), allocatable :: arr_rows2				![-]					|row array
	  real(i_kind), dimension(:), allocatable :: arr_cols3				![-]					|column array
	  real(i_kind), dimension(:), allocatable :: arr_rows3				![-]					|row array
	  integer, dimension(:), allocatable :: col_pos1					![-]					|column position array
	  integer, dimension(:), allocatable :: row_pos1					![-]					|row position array
	  integer, dimension(:), allocatable :: col_pos2					![-]					|column position array
	  integer, dimension(:), allocatable :: row_pos2					![-]					|row position array
	  real(i_kind) :: min_cols											![-]					|minimum value column
	  real(i_kind) :: min_rows											![-]					|minimum value row
	  real(i_kind) :: max_cols											![-]					|maximum value column
	  real(i_kind) :: max_rows											![-]					|maximum value row
	  real(i_kind) :: xllcorner_raster1									![-]					|x coordinate of lower left corner
	  real(i_kind) :: yllcorner_raster1									![-]					|y coordinate of lower left corner
	  real(i_kind) :: xllcorner_raster2									![-]					|x coordinate of lower left corner
	  real(i_kind) :: yllcorner_raster2									![-]					|y coordinate of lower left corner
	  real(i_kind) :: cellsize											![-]					|cellsize
	  integer :: i														![-]					|counter
	  
	  !if resolutions are unequal exit function
	  if(raster1%header%cellsize /= raster2%header%cellsize) goto 20
	  
	  !when raster header are equal, add matrix values
	  if(raster1%header%ncols == raster2%header%ncols .and. raster1%header%nrows == raster2%header%nrows .and. raster1%header%xllcorner == raster2%header%xllcorner .and. raster1%header%yllcorner == raster2%header%yllcorner .and. raster1%header%cellsize == raster2%header%cellsize) then
	    raster3%header = raster1%header
		raster3%data = raster1%data + raster2%data
	  else
	    !define cellsize
		cellsize = raster1%header%cellsize
		
		!round xllcorner and yllcorner for raster1 and raster2
		xllcorner_raster1 = real(nint(raster1%header%xllcorner*10._i_kind),i_kind)/10._i_kind
		yllcorner_raster1 = real(nint(raster1%header%yllcorner*10._i_kind),i_kind)/10._i_kind
		xllcorner_raster2 = real(nint(raster2%header%xllcorner*10._i_kind),i_kind)/10._i_kind
		yllcorner_raster2 = real(nint(raster2%header%yllcorner*10._i_kind),i_kind)/10._i_kind
		
		!construct cols and rows arrays for raster 1 and 2
	    arr_cols1 = (real([(i,i=0,raster1%header%ncols-1)],i_kind) * cellsize) + xllcorner_raster1 + (cellsize/2._i_kind)
	    arr_rows1 = (real([(i,i=0,raster1%header%nrows-1)],i_kind) * cellsize) + yllcorner_raster1 + (cellsize/2._i_kind)
	    arr_cols2 = (real([(i,i=0,raster2%header%ncols-1)],i_kind) * cellsize) + xllcorner_raster2 + (cellsize/2._i_kind)
	    arr_rows2 = (real([(i,i=0,raster2%header%nrows-1)],i_kind) * cellsize) + yllcorner_raster2 + (cellsize/2._i_kind)
	    
		!find min and max cols and rows for result raster
		min_cols = minval((/arr_cols1,arr_cols2/))
		max_cols = maxval((/arr_cols1,arr_cols2/))
		min_rows = minval((/arr_rows1,arr_rows2/))
		max_rows = maxval((/arr_rows1,arr_rows2/))
		
		!define resolution
		raster3%header%cellsize = cellsize
		
		!construct cols and rows for raster 3
		arr_cols3 = real([(i,i=int(min_cols/raster3%header%cellsize),int(max_cols/raster3%header%cellsize),1)],i_kind) * raster3%header%cellsize
	    arr_rows3 = real([(i,i=int(max_rows/raster3%header%cellsize),int(min_rows/raster3%header%cellsize),-1)],i_kind) * raster3%header%cellsize
	    
		!define header
		raster3%header%ncols = size(arr_cols3)
		raster3%header%nrows = size(arr_rows3)
		raster3%header%xllcorner = min_cols - (raster3%header%cellsize/2._i_kind)
		raster3%header%yllcorner = min_rows - (raster3%header%cellsize/2._i_kind)
		raster3%header%nodata_value = raster1%header%nodata_value
		
		!define data matrix
		allocate(raster3%data(raster3%header%nrows,raster3%header%ncols))
		raster3%data = 0._i_kind
		
		!add raster1 to raster3 at correct position
		col_pos1 = position_B_in_A(arr_cols1, arr_cols3)
		row_pos1 = position_B_in_A(arr_rows1, arr_rows3)
		raster3%data(row_pos1,col_pos1) = raster3%data(row_pos1,col_pos1) + raster1%data
		
		!add raster2 to raster3 at correct position
		col_pos2 = position_B_in_A(arr_cols2, arr_cols3)
		row_pos2 = position_B_in_A(arr_rows2, arr_rows3)
		raster3%data(row_pos2,col_pos2) = raster3%data(row_pos2,col_pos2) + raster2%data  
	  end if
	  
	  !deallocate arrays
	  if(allocated(arr_cols1)) deallocate(arr_cols1)
	  if(allocated(arr_rows1)) deallocate(arr_rows1)
	  if(allocated(arr_cols2)) deallocate(arr_cols2)
	  if(allocated(arr_rows2)) deallocate(arr_rows2)
	  if(allocated(arr_cols3)) deallocate(arr_cols3)
	  if(allocated(arr_rows3)) deallocate(arr_rows3)
	  if(allocated(col_pos1)) deallocate(col_pos1)
	  if(allocated(row_pos1)) deallocate(row_pos1)
	  if(allocated(col_pos2)) deallocate(col_pos2)
	  if(allocated(row_pos2)) deallocate(row_pos2)
	  
20	end function add_ascii_grid
	
    function add_ascii_grid_fast (raster1, raster2) result(raster3)
      !! Corrected to match codebase convention:
      !!  - (0,0) = lower-left corner
      !!  - row index grows northward (up)
      !!  - consistent with read_ascii_grid centering logic
    
      implicit none
    
      type(ascii_grid), intent(in) :: raster1
      type(ascii_grid), intent(in) :: raster2
      type(ascii_grid) :: raster3
    
      real(i_kind) :: cell
      real(i_kind) :: x_min, y_min, x_max, y_max
      integer :: nrows_out, ncols_out
      integer :: row_off1, col_off1, row_off2, col_off2
      integer :: i1_start, i1_end, j1_start, j1_end
      integer :: i2_start, i2_end, j2_start, j2_end
    
      if (abs(raster1%header%cellsize - raster2%header%cellsize) > 1.e-6_i_kind) then
         print *, '⚠️ Warning: cellsize mismatch in add_ascii_grid_fast'
      end if
    
      cell = raster1%header%cellsize
    
      !–– Determine full extent (cell centers)
      x_min = min(raster1%header%xllcorner, raster2%header%xllcorner)
      y_min = min(raster1%header%yllcorner, raster2%header%yllcorner)
      x_max = max(raster1%header%xllcorner + (raster1%header%ncols - 1)*cell, &
                  raster2%header%xllcorner + (raster2%header%ncols - 1)*cell)
      y_max = max(raster1%header%yllcorner + (raster1%header%nrows - 1)*cell, &
                  raster2%header%yllcorner + (raster2%header%nrows - 1)*cell)
    
      ncols_out = int((x_max - x_min)/cell + 1.5_i_kind)
      nrows_out = int((y_max - y_min)/cell + 1.5_i_kind)
    
      raster3%header%ncols = ncols_out
      raster3%header%nrows = nrows_out
      raster3%header%xllcorner = x_min
      raster3%header%yllcorner = y_min
      raster3%header%cellsize  = cell
      raster3%header%nodata_value = raster1%header%nodata_value
    
      allocate(raster3%data(nrows_out, ncols_out))
      raster3%data = 0._i_kind
    
      !–– Offsets using bottom-left origin (no inversion)
      row_off1 = int((raster1%header%yllcorner - raster3%header%yllcorner)/cell + 0.5_i_kind)
      col_off1 = int((raster1%header%xllcorner - raster3%header%xllcorner)/cell + 0.5_i_kind)
      row_off2 = int((raster2%header%yllcorner - raster3%header%yllcorner)/cell + 0.5_i_kind)
      col_off2 = int((raster2%header%xllcorner - raster3%header%xllcorner)/cell + 0.5_i_kind)
    
      !–– insert raster1
      i1_start = row_off1 + 1
      i1_end   = i1_start + raster1%header%nrows - 1
      j1_start = col_off1 + 1
      j1_end   = j1_start + raster1%header%ncols - 1
      raster3%data(i1_start:i1_end, j1_start:j1_end) = raster3%data(i1_start:i1_end, j1_start:j1_end) + raster1%data
    
      !–– insert raster2
      i2_start = row_off2 + 1
      i2_end   = i2_start + raster2%header%nrows - 1
      j2_start = col_off2 + 1
      j2_end   = j2_start + raster2%header%ncols - 1
      raster3%data(i2_start:i2_end, j2_start:j2_end) = raster3%data(i2_start:i2_end, j2_start:j2_end) + raster2%data
    
    end function add_ascii_grid_fast



	function rescale_ascii_grid (raster1, target_res, methode) result(raster2)
	  !! ~~~~ description ~~~~
	  !! This function will rescale a ascii grid raster
	  
	  implicit none
	  type(ascii_grid), intent(in) :: raster1							![-]					|raster
	  real(i_kind), intent(in) :: target_res							![m]					|target cellsize
	  character(*), intent(in) :: methode								![-]					|methode of scaling
	  type(ascii_grid) :: raster2										![-]					|raster
	  real(i_kind), dimension(:,:), allocatable :: tmp					![-]					|temporal matrix
	  integer :: i														![-]					|counter
	  integer :: k														![-]					|counter
	  integer :: j														![-]					|counter
	  integer :: f														![-]					|counter
	  
	  if(target_res > raster1%header%cellsize) then
	    !determine scaling factor
	    f = nint(target_res/raster1%header%cellsize)
	    
		!start writing header information
	    raster2%header%ncols = raster1%header%ncols/f
	    raster2%header%nrows = raster1%header%nrows/f
	    raster2%header%xllcorner = raster1%header%xllcorner
	    raster2%header%yllcorner = raster1%header%yllcorner
	    raster2%header%cellsize = target_res
	    raster2%header%nodata_value = -99
	    
	    !allocate data matrix size
	    allocate(raster2%data(raster2%header%nrows,raster2%header%ncols))
	    
	    !do loop
	    do i=1,raster2%header%nrows
	      do j=1,raster2%header%ncols
		    tmp = raster1%data([(k,k=(f*i-(f-1)),(f*i),1)],[(k,k=(f*j-(f-1)),(f*j),1)])
		    !calculate sum
		    if(methode == 'sum') raster2%data(i,j) = sum(tmp)
		    !calculate mena
		    if(methode == 'mean') raster2%data(i,j) = sum(tmp)/(real(size(tmp(:,1)),i_kind)**2._i_kind)		  
		  end do	    
	    end do
	  end if
	  
	  if(target_res < raster1%header%cellsize) then
	    !determine scaling factor
	    f = nint(raster1%header%cellsize/target_res)
	    
	    !start writing header information
	    raster2%header%ncols = raster1%header%ncols*f
	    raster2%header%nrows = raster1%header%nrows*f
	    raster2%header%xllcorner = raster1%header%xllcorner
	    raster2%header%yllcorner = raster1%header%yllcorner
	    raster2%header%cellsize = target_res
	    raster2%header%nodata_value = -99
	    
		!allocate data matrix size
	    allocate(raster2%data(raster2%header%nrows,raster2%header%ncols))
	    raster2%data = 0._i_kind
				
	    !do loop
	    do i=1,raster1%header%nrows
	      do j=1,raster1%header%ncols
		    !calculate reverse sum
			if(methode == 'sum')   raster2%data([(k,k=(f*i-(f-1)),(f*i),1)],[(k,k=(f*j-(f-1)),(f*j),1)]) = raster1%data(i,j)/real((f**2),i_kind)
		    !calculate reverse mean
		    if(methode == 'mean')  raster2%data([(k,k=(f*i-(f-1)),(f*i),1)],[(k,k=(f*j-(f-1)),(f*j),1)]) = raster1%data(i,j)
		  end do	    
	    end do
	  end if
	  
	  !deallocate arrays
	  if(allocated(tmp)) deallocate(tmp)
	  
	end function rescale_ascii_grid
	
	
	function rescale_ascii_grid_centered (raster_in, target_res, methode) result(raster_out)
	  !! ~~~~ description ~~~~
	  !! This function will rescale a ascii grid raster centered at the 0,0 coordinates
	  
	  use functions_module
	
	  implicit none
	  type(ascii_grid), intent(in) :: raster_in							![-]					|raster
	  real(i_kind), intent(in) :: target_res							![m]					|target cellsize
	  character(*), intent(in) :: methode								![-]					|methode of scaling
	  integer :: target_res_int											![m]					|target cellsize
	  type(ascii_grid) :: raster_out									![-]					|raster
	  type(ascii_grid) :: raster_tmp									![-]					|raster
	  real(i_kind), dimension(:,:), allocatable :: tmp1					![-]					|temporal matrix
	  real(i_kind), dimension(:,:), allocatable :: tmp2					![-]					|temporal matrix
	  real(i_kind), dimension(:,:), allocatable :: tmp					![-]					|temporal matrix
	  integer :: i														![-]					|counter
	  integer :: j														![-]					|counter
	  integer :: k														![-]					|counter
	  integer :: f														![-]					|counter
	  real(i_kind), dimension(:), allocatable :: arr_corner_out			![m]					|array of possible corners
	  real(i_kind), dimension(:), allocatable :: arr_corner_in			![m]					|array of possible corners
	  real(i_kind) :: x_corner_opt										![m]					|optimal x coordinate of lower left corner
	  real(i_kind) :: y_corner_opt										![m]					|optimal y coordinate of lower left corner
	  integer :: n_left													![-]					|counter
	  integer :: n_right												![-]					|counter
	  integer :: n_top													![-]					|counter
	  integer :: n_bottom												![-]					|counter
	  integer :: n_row													![-]					|counter
	  integer :: n_col													![-]					|counter
	  
	  !convert target res to int
	  target_res_int = int(target_res)
	  
	  !set raster out
	  raster_tmp = raster_in
	  
	  if(target_res > raster_tmp%header%cellsize) then
	    !determine scaling factor
	    f = nint(target_res/raster_tmp%header%cellsize)
	    
		!generate array of possible xllcorner and yllcorner based on output and input resolution
		arr_corner_out = - real(floor(target_res/2._i_kind),i_kind) - 1._i_kind/2._i_kind*target_res - real([(i,i=0,ceiling(max(abs(raster_tmp%header%xllcorner),abs(raster_tmp%header%yllcorner))/target_res))],i_kind)*target_res
		arr_corner_in = arr_corner_out + 0.5_i_kind*target_res - 0.5_i_kind*raster_in%header%cellsize
		
		
		!! lateral correction
		if(.not. any(abs(arr_corner_in - raster_tmp%header%yllcorner) < 1.e-5_i_kind))then
		  !determine ideal minlat
		  y_corner_opt = maxval(pack(arr_corner_in ,mask = arr_corner_in < raster_tmp%header%yllcorner))
		  n_left = int(raster_tmp%header%yllcorner-y_corner_opt)
		  
		  !determine row count of matrix
		  n_row = size(raster_tmp%data(:,1))
		  
		  !allocate filler matrix
		  if(allocated(tmp1)) deallocate(tmp1)
		  allocate(tmp1(n_row,n_left))
		  tmp1 = 0._i_kind
		  
		  !bind filler matrix on left side of matrix
		  tmp2 = cbind(tmp1,raster_tmp%data)
		  
		  !deallocate and redefine ma1
		  deallocate(raster_tmp%data)
		  raster_tmp%data = tmp2
		else
		  y_corner_opt = raster_tmp%header%yllcorner
		end if
		
		!right fill if necceressary
		n_right = target_res_int - mod(size(raster_tmp%data(1,:)),target_res_int)
		if(n_right /= 0)then
		  !determine row count of matrix
		  n_row = size(raster_tmp%data(:,1))
		  
		  !allocate filler matrix
		  if(allocated(tmp1)) deallocate(tmp1)
		  allocate(tmp1(n_row,n_right))
		  tmp1 = 0._i_kind
		  
		  !bind filler matrix on right side of matrix
		  tmp2 = cbind(raster_tmp%data,tmp1)
		  
		  !deallocate and redefine ma1
		  deallocate(raster_tmp%data)
		  raster_tmp%data = tmp2
		end if
		
		!! longitudinal correction
		if(.not. any(abs(arr_corner_in - raster_tmp%header%xllcorner) < 1.e-5_i_kind))then
		  !determine ideal minlat
		  x_corner_opt = maxval(pack(arr_corner_in ,mask = arr_corner_in < raster_tmp%header%xllcorner))
		  n_bottom = int(raster_tmp%header%xllcorner-x_corner_opt)
		  
		  !determine row count of matrix
		  n_col = size(raster_tmp%data(1,:))
		  
		  !allocate filler matrix
		  if(allocated(tmp1)) deallocate(tmp1)
		  allocate(tmp1(n_bottom,n_col))
		  tmp1 = 0._i_kind
		  
		  !bind filler matrix on bottom of matrix
		  tmp2 = rbind(raster_tmp%data,tmp1)
		  
		  !deallocate and redefine ma1
		  deallocate(raster_tmp%data)
		  raster_tmp%data = tmp2
		else
		  x_corner_opt = raster_tmp%header%xllcorner
		end if
		
		!top fill if necceressary
		n_top = target_res_int - mod(size(raster_tmp%data(:,1)),target_res_int)
		if(n_top /= 0)then
		  !determine row count of matrix
		  n_col = size(raster_tmp%data(1,:))
		  
		  !allocate filler matrix
		  if(allocated(tmp1)) deallocate(tmp1)
		  allocate(tmp1(n_top,n_col))
		  tmp1 = 0._i_kind
		  
		  !bind filler matrix on top of matrix
		  tmp2 = rbind(tmp1,raster_tmp%data)
		  
		  !deallocate and redefine ma1
		  deallocate(raster_tmp%data)
		  raster_tmp%data = tmp2
		end if
		
		! writing header information
	    raster_out%header%ncols = size(raster_tmp%data(1,:))/target_res_int
	    raster_out%header%nrows = size(raster_tmp%data(:,1))/target_res_int
	    raster_out%header%xllcorner = x_corner_opt
	    raster_out%header%yllcorner = y_corner_opt
	    raster_out%header%cellsize = target_res
	    raster_out%header%nodata_value = -99
		
		!allocate data matrix size
	    allocate(raster_out%data(raster_out%header%nrows,raster_out%header%ncols))
	    
	    !do loop
	    do i=1,raster_out%header%nrows
	      do j=1,raster_out%header%ncols
		    tmp = raster_tmp%data([(k,k=(f*i-(f-1)),(f*i),1)],[(k,k=(f*j-(f-1)),(f*j),1)])
		    !calculate sum
		    if(methode == 'sum') raster_out%data(i,j) = sum(tmp)
		    !calculate mena
		    if(methode == 'mean') raster_out%data(i,j) = sum(tmp)/(real(size(tmp(:,1)),i_kind)**2._i_kind)		  
		  end do	    
	    end do
	  end if
	  
	  
	end function rescale_ascii_grid_centered
	
	
	function add_ascii_grid_at (raster_in, matrix_in) result(raster_out)
	  !! ~~~~ description ~~~~
	  !! This function will add one ascii grid raster to another ascii grid raster at a defined position
	  
	  use functions_module
	
	  implicit none
	  !input-output
	  type(ascii_grid), intent(in) :: raster_in							![-]					|raster
	  real(i_kind), dimension(:,:), intent(in) :: matrix_in				![-]					|matrix
	  type(ascii_grid) :: raster_out									![-]					|raster
	  integer :: ncol_in												![-]					|counter
	  integer :: nrow_in												![-]					|counter
	  real(i_kind) :: min_row_in										![-]					|counter
	  real(i_kind) :: min_col_in										![-]					|counter
	  real(i_kind), dimension(:), allocatable :: arr_cols_out			![-]					|column array
	  real(i_kind), dimension(:), allocatable :: arr_rows_out			![-]					|row array
	  real(i_kind), dimension(:), allocatable :: arr_cols_in			![-]					|column array
	  real(i_kind), dimension(:), allocatable :: arr_rows_in			![-]					|row array
	  integer, dimension(:), allocatable :: col_pos						![-]					|position array
	  integer, dimension(:), allocatable :: row_pos						![-]					|position array
	  integer, dimension(:), allocatable :: col_off						![-]					|offset array
	  integer, dimension(:), allocatable :: row_off						![-]					|offset array
	  integer :: i														![-]					|counter
	  
	  !read matrix_in rows and cols
	  nrow_in = int((maxval(matrix_in(:,1))/ raster_in%header%cellsize) - (minval(matrix_in(:,1))/ raster_in%header%cellsize)) + 1
	  ncol_in = int((maxval(matrix_in(:,2))/ raster_in%header%cellsize) - (minval(matrix_in(:,2))/ raster_in%header%cellsize)) + 1
	  min_row_in = minval(matrix_in(:,1))
	  min_col_in = minval(matrix_in(:,2))
	  
	  !generate output raster
	  raster_out = create_ascii_grid( &
	  (raster_in%header%ncols + ncol_in), &
	  (raster_in%header%nrows + nrow_in), &
	  (raster_in%header%xllcorner + min_col_in), &
	  (raster_in%header%yllcorner + min_row_in), &
	  raster_in%header%cellsize)
	  
	  !construct output raster arrays
	  arr_cols_out = (real([(i,i=0,raster_out%header%ncols-1)],i_kind) * raster_out%header%cellsize) + raster_out%header%xllcorner + (raster_out%header%cellsize/2._i_kind)
	  arr_rows_out = (real([(i,i=raster_out%header%nrows-1,0,-1)],i_kind) * raster_out%header%cellsize) + raster_out%header%yllcorner + (raster_out%header%cellsize/2._i_kind)
	  
	  !construct input raster arrays
	  arr_cols_in = (real([(i,i=0,raster_in%header%ncols-1)],i_kind) * raster_in%header%cellsize) + raster_in%header%xllcorner + (raster_in%header%cellsize/2._i_kind) + min_col_in
	  arr_rows_in = (real([(i,i=0,raster_in%header%nrows-1)],i_kind) * raster_in%header%cellsize) + raster_in%header%yllcorner + (raster_in%header%cellsize/2._i_kind) + min_row_in
	  
	  !position of input arrays in output arrays
	  col_pos = position_B_in_A(arr_cols_in, arr_cols_out)
	  row_pos = position_B_in_A(arr_rows_in, arr_rows_out)
	  
	  !compute offset arrays
	  row_off = (int(matrix_in(:,1)/ raster_in%header%cellsize) - int(min_row_in/ raster_in%header%cellsize))
	  col_off = (int(matrix_in(:,2)/ raster_in%header%cellsize) - int(min_col_in/ raster_in%header%cellsize))
	  
	  !loop over landscape data matrix
	  do i=1, size(matrix_in(:,1))
		raster_out%data((row_pos - row_off(i)),(col_pos + col_off(i))) = &
		raster_out%data((row_pos - row_off(i)),(col_pos + col_off(i))) + raster_in%data
	  end do
	  
	  !deallocate
	  if(allocated(arr_cols_out)) deallocate(arr_cols_out)
	  if(allocated(arr_rows_out)) deallocate(arr_rows_out)
	  if(allocated(arr_cols_in)) deallocate(arr_cols_in)
	  if(allocated(arr_rows_in)) deallocate(arr_rows_in)
	  if(allocated(col_pos)) deallocate(col_pos)
	  if(allocated(row_pos)) deallocate(row_pos)
	  if(allocated(row_off)) deallocate(row_off)
	  if(allocated(col_off)) deallocate(col_off)
		
	end function add_ascii_grid_at
	
	
    function add_ascii_grid_at_fast (raster_in, matrix_in) result(matrix_out)
       !! --------------------------------------------------------------------
       !! Optimized version of add_ascii_grid_at
       !! Performs same operation as original but:
       !!   – avoids temporary allocations
       !!   – uses direct integer index arithmetic
       !!   – returns matrix_out directly (lat, lon, value)
       !! Works under IFX and GFORTRAN, single-thread safe.
       !! --------------------------------------------------------------------
    
       implicit none
    
       !–– input / output ––
       type(ascii_grid), intent(in) :: raster_in          				! deposition pattern
       real(i_kind), dimension(:,:), intent(in) :: matrix_in 			! placement coordinates
       real(i_kind), dimension(:,:), allocatable :: matrix_out			! flattened result
    
       !–– locals ––
       integer :: nrows_in, ncols_in, nrows_out, ncols_out
       real(i_kind) :: x_min, y_min, x_max, y_max, cell
       real(i_kind) :: xll_out, yll_out
       integer :: i, r_idx, c_idx, npoints
       real(i_kind), allocatable :: data_out(:,:)
    
       !–– grid setup ––
       cell     = raster_in%header%cellsize
       nrows_in = raster_in%header%nrows
       ncols_in = raster_in%header%ncols
    
       !–– bounds of placement matrix ––
       y_min = minval(matrix_in(:,1))
       y_max = maxval(matrix_in(:,1))
       x_min = minval(matrix_in(:,2))
       x_max = maxval(matrix_in(:,2))
    
       !–– determine output grid geometry
       nrows_out = nrows_in + int((y_max - y_min)/cell) + 1
       ncols_out = ncols_in + int((x_max - x_min)/cell) + 1
       yll_out   = raster_in%header%yllcorner + y_min
       xll_out   = raster_in%header%xllcorner + x_min
    
       !–– allocate accumulation field
       allocate(data_out(nrows_out, ncols_out))
       data_out = 0._i_kind
    
       !–– overlay each placement
       do i = 1, size(matrix_in,1)
          ! Compute target indices (y reversed like ASCII raster convention)
          r_idx = nrows_out - int((matrix_in(i,1) - y_min) / cell) - nrows_in + 1
          c_idx = int((matrix_in(i,2) - x_min) / cell) + 1
    
          ! Bounds check – ensure full patch fits in output grid
          if (r_idx >= 1 .and. r_idx + nrows_in - 1 <= nrows_out .and. &
              c_idx >= 1 .and. c_idx + ncols_in - 1 <= ncols_out) then
             data_out(r_idx:r_idx+nrows_in-1, c_idx:c_idx+ncols_in-1) = &
                  data_out(r_idx:r_idx+nrows_in-1, c_idx:c_idx+ncols_in-1) + raster_in%data
          end if
       end do
    
       !–– flatten to (y, x, value)
       npoints = nrows_out * ncols_out
       allocate(matrix_out(npoints, 3))
    
       do i = 1, npoints
          matrix_out(i,1) = yll_out + ((nrows_out - 1) - (i - 1)/ncols_out) * cell + cell/2._i_kind
          matrix_out(i,2) = xll_out + mod(i - 1, ncols_out) * cell + cell/2._i_kind
          matrix_out(i,3) = data_out((i - 1)/ncols_out + 1, mod(i - 1, ncols_out) + 1)
       end do
    
       deallocate(data_out)
    
    end function add_ascii_grid_at_fast



    function add_ascii_grid_landscape(landscape_in, drift_pattern) result(raster_out)
      !! ~~~~ description ~~~~
      !! Combines raster selection and stitching into one pass.
      !! Finds all cells with value ≈ 1 in landscape_in and
      !! adds drift_pattern at each cell position.
      !!
      !! Output: full stitched ascii_grid (same cellsize)
      !! Compatible with new bottom-up (yllcorner-based) logic.
    
      implicit none
    
      !–– input / output ––
      type(ascii_grid), intent(in) :: landscape_in
      type(ascii_grid), intent(in) :: drift_pattern
      type(ascii_grid)             :: raster_out
    
      !–– locals ––
      real(i_kind) :: cell
      integer :: nrows_land, ncols_land
      integer :: nrows_pat,  ncols_pat
      real(i_kind) :: xll_min, yll_min, xll_max, yll_max
      real(i_kind) :: xpat_min, ypat_min
	  real(i_kind) :: x_pos, y_pos
      integer :: nrows_out, ncols_out
      integer :: r, c, i_add
      integer :: r_start, r_end, c_start, c_end
      integer :: hits
      
      !–– setup
      cell       = landscape_in%header%cellsize
      nrows_land = landscape_in%header%nrows
      ncols_land = landscape_in%header%ncols
      nrows_pat  = drift_pattern%header%nrows
      ncols_pat  = drift_pattern%header%ncols
    
      ! --- compute drift pattern 0,0 offset
	  xpat_min = drift_pattern%header%xllcorner + 0.5_i_kind*cell
	  ypat_min = drift_pattern%header%yllcorner + 0.5_i_kind*cell
	  
	  !–– compute total bounds (extend by pattern size)
      xll_min = landscape_in%header%xllcorner + xpat_min
      yll_min = landscape_in%header%yllcorner + ypat_min
      xll_max = xll_min + (ncols_land+ncols_pat)*cell
      yll_max = yll_min + (nrows_land+nrows_pat)*cell
    
      ncols_out = int((xll_max - xll_min)/cell + 1.5_i_kind)
      nrows_out = int((yll_max - yll_min)/cell + 1.5_i_kind)
    
      raster_out%header%ncols = ncols_out
      raster_out%header%nrows = nrows_out
      raster_out%header%xllcorner = xll_min
      raster_out%header%yllcorner = yll_min
      raster_out%header%cellsize  = cell
      raster_out%header%nodata_value = landscape_in%header%nodata_value
    
      allocate(raster_out%data(nrows_out, ncols_out))
      raster_out%data = 0._i_kind
    
      !–– main loop: place pattern wherever landscape cell ≈ 1
      do r = 1, nrows_land
         do c = 1, ncols_land
            if (abs(landscape_in%data(r,c) - 1._i_kind) < 1.e-5_i_kind) then
               ! compute physical coordinates of this cell’s center
               y_pos = landscape_in%header%yllcorner + (r - 0.5_i_kind) * cell
               x_pos = landscape_in%header%xllcorner + (c - 0.5_i_kind) * cell
      
               ! translate to output indices
               r_start = int(((y_pos-(yll_min+0.5_i_kind*cell))/cell) + (ypat_min / cell)) + 1
               c_start = int(((x_pos-(xll_min+0.5_i_kind*cell))/cell) + (xpat_min / cell)) + 1
			   r_end   = r_start + nrows_pat - 1
               c_end   = c_start + ncols_pat - 1
      
               if (r_start >= 1 .and. r_end <= nrows_out .and. c_start >= 1 .and. c_end <= ncols_out) then
                  raster_out%data(r_start:r_end, c_start:c_end) = raster_out%data(r_start:r_end, c_start:c_end) + drift_pattern%data
               end if
            end if
         end do
      end do
    
    end function add_ascii_grid_landscape

!    function rescale_drift_pattern_fast(pattern_in, target_cellsize) result(pattern_out)
!      !! ---------------------------------------------------------------------------
!      !! Fast rescaler for drift pattern rasters
!      !! Replaces the old "create_ascii_grid → add → rescale" pipeline.
!      !! Works directly in raster space using integer index arithmetic.
!      !! ---------------------------------------------------------------------------
!    
!      implicit none
!      type(ascii_grid), intent(in) :: pattern_in
!      real(i_kind), intent(in)     :: target_cellsize
!      type(ascii_grid)             :: pattern_out
!    
!      ! local
!      integer :: scale_factor, i, j
!      integer :: new_nrows, new_ncols
!      real(i_kind), allocatable :: tmp_data(:,:)
!      real(i_kind) :: ratio
!    
!      !–– if already correct resolution, just copy ––
!      if (abs(pattern_in%header%cellsize - target_cellsize) < 1.e-6_i_kind) then
!         pattern_out = pattern_in
!         return
!      end if
!    
!      !–– scaling ratio ––
!      ratio = pattern_in%header%cellsize / target_cellsize
!      if (ratio < 1._i_kind) ratio = 1._i_kind / ratio    ! ensure >= 1
!    
!      !–– new raster size ––
!      new_nrows = int(real(pattern_in%header%nrows, i_kind) * ratio)
!      new_ncols = int(real(pattern_in%header%ncols, i_kind) * ratio)
!    
!      allocate(tmp_data(new_nrows, new_ncols))
!      tmp_data = 0._i_kind
!    
!      !–– aggregate / distribute deposition values ––
!      ! use conservative "sum" scaling — same as rescale_ascii_grid_centered
!      do i = 1, pattern_in%header%nrows
!         do j = 1, pattern_in%header%ncols
!            ! compute target block bounds in output grid
!            integer :: r1, r2, c1, c2
!            r1 = int((i-1) * ratio) + 1
!            r2 = min(int(i * ratio), new_nrows)
!            c1 = int((j-1) * ratio) + 1
!            c2 = min(int(j * ratio), new_ncols)
!            tmp_data(r1:r2, c1:c2) = tmp_data(r1:r2, c1:c2) + pattern_in%data(i,j)
!         end do
!      end do
!    
!      !–– assemble header ––
!      pattern_out%header = pattern_in%header
!      pattern_out%header%cellsize = target_cellsize
!      pattern_out%header%nrows = new_nrows
!      pattern_out%header%ncols = new_ncols
!      allocate(pattern_out%data(new_nrows, new_ncols))
!      pattern_out%data = tmp_data
!    
!      deallocate(tmp_data)
!    
!    end function rescale_drift_pattern_fast

	
end module ascii_grid_module