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

program main
  use data_module
  
  implicit none
  
  !! Copyright disclaimer 
  write(*,'(A)') ""
  write(*,'(A)') "Copyright (C) 2024 BASF SE"
  write(*,'(A)') ""
  write(*,'(A)') "This program is free software; you can redistribute it and/or"
  write(*,'(A)') "modify it under the terms of the GNU General Public License"
  write(*,'(A)') "as published by the Free Software Foundation; version 2."
  write(*,'(A)') ""
  write(*,'(A)') "This program is distributed in the hope that it will be useful,"
  write(*,'(A)') "but WITHOUT ANY WARRANTY; without even the implied warranty of"
  write(*,'(A)') "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
  write(*,'(A)') "GNU General Public License for more details."
  write(*,'(A)') ""
  write(*,'(A)') "You should have received a copy of the GNU General Public License"
  write(*,'(A)') "along with this program; If not, see <https://www.gnu.org/licenses/>."
  write(*,'(A)') ""
  
  !reading input files
  call read_input
  
  !initialisation of environmental parameters
  call environment_initialisation
  
  !initialisation of application parameters
  call application_initialisation
  
  !run droplet model
  call droplet_model
  
  !predict drift pattern
  call predict_drift_pattern

  !predict landscape drift
  if (control_dat%code1 == 0) then
    call landscape_drift
  end if
  
  !predict drift curve
  if (control_dat%code1 == 1) then
    call drift_curve
  end if
  
  !write outputs
  call write_output
  
end program main