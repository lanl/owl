!
! © 2025. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!

module elastic_vhtiort_2d

    use libflit
    use mod_fdcoef
    use mod_su
    use mod_grid
    use mod_source_receiver, only: source_receiver_geometry

    implicit none

    type wave_solver_elastic_vhtiort_2d

        integer :: nx = 0, nz = 0
        real :: dx = 0, dz = 0
        real :: ox = 0, oz = 0
        real :: tmax = 1.0
        real :: dt = 0
        real :: data_tmax = 1.0
        real :: data_dt = 0
        integer :: pml = 15
        logical :: free_surface = .false.
        real, allocatable, dimension(:, :) :: vp, vs, rho
        real, allocatable, dimension(:, :) :: tieps, tidel, tithe, tieta
        real, allocatable, dimension(:, :) :: c11, c13, c33, c55
        character(len=1024) :: dir_synthetic = './'
        character(len=1024) :: dir_adjoint = './'
        character(len=1024) :: dir_snapshot = './'
        character(len=1024) :: dir_working = './'
        integer :: cc_step_interval = 1
        logical :: verbose = .false.
        real, allocatable, dimension(:) :: snaps
        type(source_receiver_geometry) :: gmtr
        logical :: energy_precond = .false.
        logical :: reconstruct = .false.

        character(len=12) :: anisotropy_type = 'iso'
        character(len=24) :: kernel_v = 'full'
        character(len=24) :: kernel_a = 'full'

        type(su) :: seis_vx, seis_vz, seis_vxr, seis_vzr

        logical :: compx = .true.
        logical :: compz = .true.

        real :: free_surface_dz_refine = 4.0
        real :: dz_max = 0

        integer :: nc_mt = 6
        real, allocatable, dimension(:) :: mt
        real, allocatable, dimension(:) :: stf

        logical :: yn_update_medium = .true.
        logical :: yn_update_source = .false.

    contains

        procedure, public :: forward => compute_forward
#ifdef _fwi_
        procedure, public :: adjoint => compute_adjoint
#endif
    end type


    interface
        module subroutine compute_forward(this)
            class(wave_solver_elastic_vhtiort_2d), intent(inout) :: this
        end subroutine
#ifdef _fwi_
        module subroutine compute_adjoint(this)
            class(wave_solver_elastic_vhtiort_2d), intent(inout) :: this
        end subroutine
#endif
    end interface


    private
    public :: wave_solver_elastic_vhtiort_2d

end module elastic_vhtiort_2d
