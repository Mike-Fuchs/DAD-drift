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

module constants_module
  !! ~~~~ description ~~~~
  !! In this module universal constants are defined.
    
  implicit none
  integer, parameter :: i_kind = 		selected_real_kind(15,307)		!						|precision of real(i_kind)s
  real(i_kind), parameter :: R	= 		8.31446261815324_i_kind			![J/(mol*K)]			|ideal (universal) gas constant
  real(i_kind), parameter :: g	= 		9.80665_i_kind					![m/s²]					|earth-surface gravitational acceleration
  real(i_kind), parameter :: K_0 = 		273.15_i_kind					![K]					|0°C in kelvin
  real(i_kind), parameter :: mm_da = 	28.9652E-3_i_kind				![kg/mol]				|molar mass of dry air
  real(i_kind), parameter :: mm_h2o = 	18.016E-3_i_kind				![kg/mol]				|molar mass of water
  real(i_kind), parameter :: pi	= 		4._i_kind*ATAN(1._i_kind)		![-]					|value pi
  real(i_kind), parameter :: R_h2o = 	R/mm_h2o						![J/(kg*K)]				|specific gas constant for water
  real(i_kind), parameter :: R_da = 	R/mm_da							![J/(kg*K)]				|specific gas constant for dry air
  real(i_kind), parameter :: Cp_air = 	1003.5_i_kind					![J/(kg*K)]				|specific heat capacity of air
  real(i_kind), parameter :: Cp_h2o = 	4181.3_i_kind					![J/(kg*K)]				|specific heat capacity of water
  real(i_kind), parameter :: kB = 		1.380649E-23_i_kind				![J/K]					|Boltzmann constant
  real(i_kind), parameter :: Kc = 		0.4_i_kind						![-]					|von Karman constant
  
  end module constants_module