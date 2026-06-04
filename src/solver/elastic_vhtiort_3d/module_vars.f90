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


module elastic_vhtiort_3d_vars

    use libflit
    use mod_fdcoef
    use mod_su
    use mod_source_receiver, only: source_receiver_geometry
    use mod_utility, only: check_dt_f0
    use mod_source_receiver, only : nkw
    use mod_anisotropy, only: thomsen_to_cij, alkhalifah_tsvankin_to_cij

    use elastic_vhtiort_3d

    implicit none

    real, allocatable, dimension(:, :, :) :: vp, vs, rho, tieps, tidel, tigam, tithe, tiphi, tieta

    real, allocatable, dimension(:, :, :) :: vx, vy, vz, stressxx, stressyy, stresszz, stressyz, stressxz, stressxy
    real, allocatable, dimension(:, :, :) :: memory_pdxxx, memory_pdyxy, memory_pdzxz
    real, allocatable, dimension(:, :, :) :: memory_pdxxy, memory_pdyyy, memory_pdzyz
    real, allocatable, dimension(:, :, :) :: memory_pdxxz, memory_pdyyz, memory_pdzzz
    real, allocatable, dimension(:, :, :) :: memory_pdxvx, memory_pdyvx, memory_pdzvx
    real, allocatable, dimension(:, :, :) :: memory_pdxvy, memory_pdyvy, memory_pdzvy
    real, allocatable, dimension(:, :, :) :: memory_pdxvz, memory_pdyvz, memory_pdzvz

    real, allocatable, dimension(:, :, :) :: vxr, vyr, vzr, stressxxr, stressyyr, stresszzr, stressyzr, stressxzr, stressxyr
    real, allocatable, dimension(:, :, :) :: memory_pdxxxr, memory_pdyxyr, memory_pdzxzr
    real, allocatable, dimension(:, :, :) :: memory_pdxxyr, memory_pdyyyr, memory_pdzyzr
    real, allocatable, dimension(:, :, :) :: memory_pdxxzr, memory_pdyyzr, memory_pdzzzr
    real, allocatable, dimension(:, :, :) :: memory_pdxvxr, memory_pdyvxr, memory_pdzvxr
    real, allocatable, dimension(:, :, :) :: memory_pdxvyr, memory_pdyvyr, memory_pdzvyr
    real, allocatable, dimension(:, :, :) :: memory_pdxvzr, memory_pdyvzr, memory_pdzvzr

    real, allocatable, dimension(:, :, :) :: energy_src_v, energy_rec_v
    real, allocatable, dimension(:, :, :) :: energy_src_a, energy_rec_a

    real, allocatable, dimension(:, :, :) :: snapvx, snapvy, snapvz
    real, allocatable, dimension(:, :, :) :: snapvxr, snapvyr, snapvzr

    real, allocatable, dimension(:, :, :) :: prev_vx, prev_vy, prev_vz
    real, allocatable, dimension(:, :, :) :: prev_stressxx, prev_stressyy, prev_stresszz
    real, allocatable, dimension(:, :, :) :: prev_stressyz, prev_stressxz, prev_stressxy

    real, allocatable, dimension(:, :, :) :: prev_pr, prev_stressxxr, prev_stressyyr, prev_stresszzr
    real, allocatable, dimension(:, :, :) :: prev_stressyzr, prev_stressxzr, prev_stressxyr

    real, allocatable, dimension(:, :, :) :: strainxx, strainyy, strainzz, strainxz, strainyz, strainxy
    real, allocatable, dimension(:, :, :) :: strainxxr, strainyyr, strainzzr, strainxzr, strainyzr, strainxyr
    real, allocatable, dimension(:, :, :) :: src_vx, src_vy, src_vz, rec_vx, rec_vy, rec_vz

    real, allocatable, dimension(:, :, :) :: strainxx_hilbert, strainyy_hilbert, strainzz_hilbert, strainxz_hilbert, strainyz_hilbert, strainxy_hilbert
    real, allocatable, dimension(:, :, :) :: strainxxr_hilbert, strainyyr_hilbert, strainzzr_hilbert, strainxzr_hilbert, strainyzr_hilbert, strainxyr_hilbert
    real, allocatable, dimension(:, :, :) :: src_hilbert, rec_hilbert

    real, allocatable, dimension(:, :, :) :: s11, s12, s13, s22, s23, s33, s44, s55, s66

    real, allocatable, dimension(:, :, :) :: c11, c12, c13, c14, c15, c16
    real, allocatable, dimension(:, :, :) :: c22, c23, c24, c25, c26
    real, allocatable, dimension(:, :, :) :: c33, c34, c35, c36
    real, allocatable, dimension(:, :, :) :: c44, c45, c46
    real, allocatable, dimension(:, :, :) :: c55, c56
    real, allocatable, dimension(:, :, :) :: c66

    real, allocatable, dimension(:, :, :) :: grad_c11, grad_c12, grad_c13
    real, allocatable, dimension(:, :, :) :: grad_c22, grad_c23, grad_c33
    real, allocatable, dimension(:, :, :) :: grad_c44, grad_c55, grad_c66
    real, allocatable, dimension(:, :, :) :: grad_vp, grad_vs, grad_eps, grad_del, grad_gam, grad_rho, grad_eta
    real, allocatable, dimension(:, :, :) :: energy_src, energy_rec

    real, allocatable, dimension(:) :: snaps

    integer :: np

    real :: idx, idy, idz

    integer :: nx, ny, nz
    real :: dx, dy, dz
    real :: ox, oy, oz

    integer :: nt
    real :: dt
    real :: tmax

    integer :: data_nt
    real :: data_dt
    real :: data_tmax

    integer :: pml

    logical :: yn_reconstruct = .false.
    logical :: yn_free_surface = .false.
    type(source_receiver_geometry) :: sgmtr

    integer :: cc_step_interval
    logical :: verbose = .false.

    logical :: yn_energy_precond = .true.

    character(len=1024) :: dir_synthetic, dir_snapshot, dir_working, dir_adjoint

    logical :: yn_compx = .true.
    logical :: yn_compy = .true.
    logical :: yn_compz = .true.

    character(len=12) :: aniso_param = 'iso'

    integer :: nx1, nx2, ny1, ny2, nz1, nz2
    integer :: nx1_interior, nx2_interior, ny1_interior, ny2_interior, nz1_interior, nz2_interior

    real :: pmlvp

    real, allocatable, dimension(:) :: zz_i, zz_h, dz_i, dz_h
    real, allocatable, dimension(:) :: dz_scaling_i, dz_scaling_h

    integer :: htlen = 30

    integer, allocatable, dimension(:, :) :: trace_range

    integer :: nc_mt
    real, allocatable, dimension(:) :: grad_mt
    logical :: yn_update_medium = .true.
    logical :: yn_update_source = .false.
    real, allocatable, dimension(:, :) :: dstf_dt

contains

    !
    !> Map the model on the original regular mesh to the depth-varying mesh
    !
    subroutine map_regular_to_irregular(v, this, rr)

        real, allocatable, dimension(:, :, :), intent(inout) :: v
        type(wave_solver_elastic_vhtiort_3d), intent(in) :: this
        integer, dimension(1:6), intent(in) :: rr

        real, allocatable, dimension(:, :, :) :: m
        real, allocatable, dimension(:) :: z
        integer, allocatable, dimension(:) :: r
        integer :: i, j

        r = zeros(6)
        r(1) = lbound(v, dim=1)
        r(2) = ubound(v, dim=1)
        r(3) = lbound(v, dim=2)
        r(4) = ubound(v, dim=2)
        r(5) = lbound(v, dim=3)
        r(6) = ubound(v, dim=3)

        ! Depth coordinates of regular-mesh array
        z = regspace(r(5) - 1.0, 1.0, r(6) - 1.0)*this%dz

        ! Interpolate to irregular-mesh array
        call alloc_array(m, rr)
        !$omp parallel do private(i, j) collapse(2)
        do j = max(r(3), rr(3)), min(r(4), rr(4))
            do i = max(r(1), rr(1)), min(r(2), rr(2))
                m(i, j, :) = ginterp(z, v(i, j, r(5):r(6)), zz_i, 'linear')
            end do
        end do
        !$omp end parallel do

        ! Overwrite original
        v = m

    end subroutine

    !
    !> Map the model on the depth-varying mesh to the original regular mesh
    !
    subroutine map_irregular_to_regular(v, this, rr)

        real, allocatable, dimension(:, :, :), intent(inout) :: v
        type(wave_solver_elastic_vhtiort_3d), intent(in) :: this
        integer, dimension(1:6), intent(in) :: rr

        real, allocatable, dimension(:, :, :) :: m
        real, allocatable, dimension(:) :: zz
        integer, allocatable, dimension(:) :: r
        integer :: i, j

        r = zeros(6)
        r(1) = lbound(v, dim=1)
        r(2) = ubound(v, dim=1)
        r(3) = lbound(v, dim=2)
        r(4) = ubound(v, dim=2)
        r(5) = lbound(v, dim=3)
        r(6) = ubound(v, dim=3)

        zz = regspace(rr(5) - 1.0, 1.0, rr(6) - 1.0)*this%dz

        call alloc_array(m, rr)
        !$omp parallel do private(i, j) collapse(2)
        do j = max(r(3), rr(3)), min(r(4), rr(4))
            do i = max(r(1), rr(1)), min(r(2), rr(2))
                m(i, j, :) = ginterp(zz_i(r(5):r(6)), v(i, j, r(5):r(6)), zz, 'cubic')
            end do
        end do
        !$omp end parallel do

        v = m

    end subroutine

    !
    !> Prepare for modeling
    !
    subroutine prepare_modeling(this)

        type(wave_solver_elastic_vhtiort_3d), intent(in) :: this

        real :: minv, maxv, dtstable, f0clean, temp
        integer :: i, j, k, l, refine_nz
        real :: rr, dz_max, rayleigh_wavelength
        real, allocatable, dimension(:, :, :) :: mdispersion, mstability
        integer :: nbeg, nend

        nx = this%nx
        ny = this%ny
        nz = this%nz
        dx = this%dx
        dy = this%dy
        dz = this%dz
        ox = this%ox
        oy = this%oy
        oz = this%oz

        idx = 1.0d0/dx
        idy = 1.0d0/dy
        idz = 1.0d0/dz

        pml = this%pml

        ! Allocate memory for Cij models
        c11 = zeros(nx, ny, nz)
        c12 = zeros(nx, ny, nz)
        c13 = zeros(nx, ny, nz)
        c14 = zeros(nx, ny, nz)
        c15 = zeros(nx, ny, nz)
        c16 = zeros(nx, ny, nz)
        c22 = zeros(nx, ny, nz)
        c23 = zeros(nx, ny, nz)
        c24 = zeros(nx, ny, nz)
        c25 = zeros(nx, ny, nz)
        c26 = zeros(nx, ny, nz)
        c33 = zeros(nx, ny, nz)
        c34 = zeros(nx, ny, nz)
        c35 = zeros(nx, ny, nz)
        c36 = zeros(nx, ny, nz)
        c44 = zeros(nx, ny, nz)
        c45 = zeros(nx, ny, nz)
        c46 = zeros(nx, ny, nz)
        c55 = zeros(nx, ny, nz)
        c56 = zeros(nx, ny, nz)
        c66 = zeros(nx, ny, nz)
        rho = zeros(nx, ny, nz)

        aniso_param = this%anisotropy_type
        select case (aniso_param)

            case ('iso')
                vp = permute(this%vp, 321)
                vs = permute(this%vs, 321)
                rho = permute(this%rho, 321)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                c11 = rho*vp**2
                c22 = c11
                c33 = c11
                c44 = rho*vs**2
                c55 = c44
                c66 = c44
                c12 = c11 - 2*c55
                c13 = c12
                c23 = c12

                minv = minval(vs)
                maxv = maxval(vp)

            case ('thomsen')
                vp = permute(this%vp, 321)
                vs = permute(this%vs, 321)
                tieps = permute(this%tieps, 321)
                tidel = permute(this%tidel, 321)
                tigam = permute(this%tigam, 321)
                tithe = permute(this%tithe, 321)
                tiphi = permute(this%tiphi, 321)
                rho = permute(this%rho, 321)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                !$omp parallel do private(i, j, k, temp) collapse(3)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            call thomsen_to_cij(vp(i, j, k), vs(i, j, k), rho(i, j, k), &
                                tieps(i, j, k), tidel(i, j, k), tigam(i, j, k), tithe(i, j, k), tiphi(i, j, k), &
                                c11(i, j, k), c12(i, j, k), c13(i, j, k), temp, temp, temp, &
                                c22(i, j, k), c23(i, j, k), temp, temp, temp, &
                                c33(i, j, k), temp, temp, temp, &
                                c44(i, j, k), temp, temp, &
                                c55(i, j, k), temp, &
                                c66(i, j, k))
                        end do
                    end do
                end do
                !$omp end parallel do

            case ('a-t')
                vp = permute(this%vp, 321)
                vs = permute(this%vs, 321)
                tieps = permute(this%tieps, 321)
                tieta = permute(this%tieta, 321)
                tigam = permute(this%tigam, 321)
                tithe = permute(this%tithe, 321)
                tiphi = permute(this%tiphi, 321)
                rho = permute(this%rho, 321)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                !$omp parallel do private(i, j, k, temp) collapse(3)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            call alkhalifah_tsvankin_to_cij(vp(i, j, k), vs(i, j, k), rho(i, j, k), &
                                tieps(i, j, k), tieta(i, j, k), tigam(i, j, k), tithe(i, j, k), tiphi(i, j, k), &
                                c11(i, j, k), c12(i, j, k), c13(i, j, k), temp, temp, temp, &
                                c22(i, j, k), c23(i, j, k), temp, temp, temp, &
                                c33(i, j, k), temp, temp, temp, &
                                c44(i, j, k), temp, temp, &
                                c55(i, j, k), temp, &
                                c66(i, j, k))
                        end do
                    end do
                end do
                !$omp end parallel do

            case ('cij')
                c11 = permute(this%c11, 321)
                c12 = permute(this%c12, 321)
                c13 = permute(this%c13, 321)
                c22 = permute(this%c22, 321)
                c23 = permute(this%c23, 321)
                c33 = permute(this%c33, 321)
                c44 = permute(this%c44, 321)
                c55 = permute(this%c55, 321)
                c66 = permute(this%c66, 321)
                rho = permute(this%rho, 321)

                call assert(maxval(abs(c11)) > 0, ' Error: c11 is all zero.')
                call assert(maxval(abs(c22)) > 0, ' Error: c22 is all zero.')
                call assert(maxval(abs(c33)) > 0, ' Error: c33 is all zero.')
                call assert(maxval(abs(c44)) > 0, ' Error: c44 is all zero.')
                call assert(maxval(abs(c55)) > 0, ' Error: c55 is all zero.')
                call assert(maxval(abs(c66)) > 0, ' Error: c66 is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

        end select

        ! PML Vp
        pmlvp = max( &
            maxval(c11(1, :, :)/rho(1, :, :)), &
            maxval(c11(nx, :, :)/rho(nx, :, :)), &
            maxval(c11(:, 1, :)/rho(:, 1, :)), &
            maxval(c11(:, ny, :)/rho(:, ny, :)), &
            maxval(c11(:, :, 1)/rho(:, :, 1)), &
            maxval(c11(:, :, nz)/rho(:, :, nz)))
        pmlvp = max(pmlvp, &
            maxval(c22(1, :, :)/rho(1, :, :)), &
            maxval(c22(nx, :, :)/rho(nx, :, :)), &
            maxval(c22(:, 1, :)/rho(:, 1, :)), &
            maxval(c22(:, ny, :)/rho(:, ny, :)), &
            maxval(c22(:, :, 1)/rho(:, :, 1)), &
            maxval(c22(:, :, nz)/rho(:, :, nz)))
        pmlvp = max(pmlvp, &
            maxval(c33(1, :, :)/rho(1, :, :)), &
            maxval(c33(nx, :, :)/rho(nx, :, :)), &
            maxval(c33(:, 1, :)/rho(:, 1, :)), &
            maxval(c33(:, ny, :)/rho(:, ny, :)), &
            maxval(c33(:, :, 1)/rho(:, :, 1)), &
            maxval(c33(:, :, nz)/rho(:, :, nz)))
        pmlvp = sqrt(pmlvp)

        ! Pad models with PML
        call pad_array(c11, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c12, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c13, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c22, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c23, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c33, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c44, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c55, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c66, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(rho, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])

        ! Refine near-surface mesh along z when doing free-surface modeling
        yn_free_surface = this%free_surface
        if (yn_free_surface) then

            ! Approximate Rayleigh wavelength at the free surface
            rayleigh_wavelength = sqrt(max(mean(c44(:, :, 1)/rho(:, :, 1)), &
                mean(c55(:, :, 1)/rho(:, :, 1)), &
                mean(c66(:, :, 1)/rho(:, :, 1))))/minval(this%gmtr%srcr(:)%f0)
            refine_nz = ceiling(0.5*rayleigh_wavelength/dz)

            ! Compute integer grid point positions
            ! Initial dz = dz/dz_refine, and it gradually increases after refine_nz, by 2.5% per grid
            ! 5% can cause very weak artificial reflections, while 1% is just too costly.
            rr = this%free_surface_dz_refine
            dz_max = dz/rr
            zz_i = regspace(0.0, dz_max, (refine_nz - 1)*dz)
            l = size(zz_i)
            do while (zz_i(l) < (this%nz - 1)*dz)
                dz_max = min(dz_max*1.025, this%dz_max)
                zz_i = [zz_i, zz_i(l) + dz_max]
                l = l + 1
            end do

            ! New nz is the length of zz_i
            nz = size(zz_i)

            ! Half integer grid points are interpolated
            zz_h = zeros(nz)
            zz_h(1) = -0.5*dz/rr
            do i = 1, nz - 1
                zz_h(i + 1) = 0.5*(zz_i(i) + zz_i(i + 1))
            end do

            ! Depth varying grid size
            dz_i = deriv(zz_i, method='center')
            dz_scaling_i = dz/dz_i

            dz_h = deriv(zz_h, method='center')
            dz_scaling_h = dz/dz_h

            call pad_array(dz_i, [pml + 1, pml + 1])
            call pad_array(dz_h, [pml + 1, pml + 1])
            call pad_array(dz_scaling_i, [pml + 1, pml + 1])
            call pad_array(dz_scaling_h, [pml + 1, pml + 1])
            call pad_array(zz_i, [pml + 1, pml + 1])
            call pad_array(zz_h, [pml + 1, pml + 1])

            do i = 1, pml + 1
                zz_i(1 - i) = zz_i(1) - i*dz/rr
                zz_i(nz + i) = zz_i(nz) + i*dz_max
                zz_h(1 - i) = zz_h(1) - i*dz/rr
                zz_h(nz + i) = zz_h(nz) + i*dz_max
            end do

            ! Interpolate medium parameter models
            call map_regular_to_irregular(c11, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c12, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c13, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c22, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c23, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c33, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c44, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c55, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c66, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(rho, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])

            if (rankid_group == 0) then
                call warn(date_time_compact()//' Original nz = '//num2str(this%nz)//', adjusted nz = '//num2str(nz))
                call warn(date_time_compact()//' min(dz) = '//num2str(min(minval(dz_i), minval(dz_h)), '(es)') &
                    //', max(dz) = '//num2str(max(maxval(dz_i), maxval(dz_h)), '(es)'))
            end if

        end if

        call mpibarrier_group

        ! Check stability and dispersion
        dt = this%dt
        tmax = this%tmax
        temp = sum(abs(fdcoefs))

        if (yn_free_surface) then
            ! For depth-varying mesh, the stability & dispersion criteria vary in space as well

            mdispersion = zeros(nx, ny, nz)
            mstability = zeros(nx, ny, nz)

            !$omp parallel do private(i, j, k, minv, maxv) collapse(3)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx

                        minv = min(sqrt(c44(i, j, k)/rho(i, j, k)), sqrt(c55(i, j, k)/rho(i, j, k)), sqrt(c66(i, j, k)/rho(i, j, k)))
                        maxv = max(sqrt(c11(i, j, k)/rho(i, j, k)), sqrt(c22(i, j, k)/rho(i, j, k)), sqrt(c33(i, j, k)/rho(i, j, k)))

                        mstability(i, j, k) = 1.0/(temp*maxv*sqrt(1.0/dx**2 + 1.0/dy**2 + 1.0/min(dz_i(k), dz_h(k))**2))
                        ! Empirical relation to ensure sufficient accuracy for surface waves
                        if (sum(dz_i(1:k)) <= 0.5*rayleigh_wavelength) then
                            mdispersion(i, j, k) = min(0.9*minv/7.0/max(dx, dy), 0.2*0.9*minv/10.0/max(dz_i(k), dz_h(k)))
                        else
                            mdispersion(i, j, k) = minv/7.0/max(dx, dy, max(dz_i(k), dz_h(k)))
                        end if

                    end do
                end do
            end do
            !$omp end parallel do

            dtstable = minval(mstability)
            f0clean = minval(mdispersion)

        else
            ! Fixed mesh

            minv = min(minval(sqrt(c44/rho)), minval(sqrt(c55/rho)), minval(sqrt(c66/rho)))
            maxv = max(maxval(sqrt(c11/rho)), maxval(sqrt(c22/rho)), maxval(sqrt(c33/rho)))

            dtstable = 1.0/(temp*maxv*sqrt(1.0/dx**2 + 1.0/dy**2 + 1.0/dz**2))
            f0clean = minv/max(dx, dy, dz)/7.0

        end if

        if (rankid_group == 0) then
            call warn(date_time_compact()//' Stable dt = '//num2str(dtstable, '(es)') &
                //' s, clean f0 = '//num2str(f0clean, '(es)')//' Hz')
        end if

        call check_dt_f0(dt, dtstable, maxval(this%gmtr%srcr(:)%f0), f0clean)

        nt = nint(tmax/dt + 1)

        ! Domain decomposition
        if (yn_free_surface) then

            call domain_decomp_regular_group(nx + 2*pml, ny + 2*pml, nz + pml, nx1, nx2, ny1, ny2, nz1, nz2, &
                weights1=[ones(pml)*1.25, ones(nx), ones(pml)*1.25], &
                weights2=[ones(pml)*1.25, ones(ny), ones(pml)*1.25], &
                weights3=[ones(nz), ones(pml)*1.25])

            nx1 = nx1 - pml
            nx2 = nx2 - pml
            ny1 = ny1 - pml
            ny2 = ny2 - pml

        else

            call domain_decomp_regular_group(nx + 2*pml, ny + 2*pml, nz + 2*pml, nx1, nx2, ny1, ny2, nz1, nz2, &
                weights1=[ones(pml)*1.25, ones(nx), ones(pml)*1.25], &
                weights2=[ones(pml)*1.25, ones(ny), ones(pml)*1.25], &
                weights3=[ones(pml)*1.25, ones(nz), ones(pml)*1.25])

            nx1 = nx1 - pml
            nx2 = nx2 - pml
            ny1 = ny1 - pml
            ny2 = ny2 - pml
            nz1 = nz1 - pml
            nz2 = nz2 - pml

        end if
        nx1_interior = max(1, nx1)
        nx2_interior = min(nx2, nx)
        ny1_interior = max(1, ny1)
        ny2_interior = min(ny2, ny)
        nz1_interior = max(1, nz1)
        nz2_interior = min(nz2, nz)

        ! Crop models
        call alloc_array(c11, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c11(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c12, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c12(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c13, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c13(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c22, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c22(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c23, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c23(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c33, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c33(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c44, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c44(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c55, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c55(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c66, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=c66(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(rho, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=rho(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))

        ! Prepare geometry
        sgmtr = this%gmtr
        sgmtr%nx = nx
        sgmtr%ny = ny
        sgmtr%nz = this%nz
        sgmtr%dx = dx
        sgmtr%dy = dy
        sgmtr%dz = this%dz
        sgmtr%ox = ox
        sgmtr%oy = oy
        sgmtr%oz = oz
        sgmtr%xmin = max(ox, this%gmtr%xmin)
        sgmtr%xmax = min((nx - 1)*dx + ox, this%gmtr%xmax)
        sgmtr%ymin = max(oy, this%gmtr%ymin)
        sgmtr%ymax = min((ny - 1)*dy + oy, this%gmtr%ymax)
        sgmtr%zmin = max(oz, this%gmtr%zmin)
        sgmtr%zmax = min((this%nz - 1)*this%dz + oz, this%gmtr%zmax)
        sgmtr%sxmin = max(ox, this%gmtr%sxmin)
        sgmtr%sxmax = min((nx - 1)*dx + ox, this%gmtr%sxmax)
        sgmtr%symin = max(oy, this%gmtr%symin)
        sgmtr%symax = min((ny - 1)*dy + oy, this%gmtr%symax)
        sgmtr%szmin = max(oz, this%gmtr%szmin)
        sgmtr%szmax = min((this%nz - 1)*this%dz + oz, this%gmtr%szmax)
        sgmtr%rxmin = max(ox, this%gmtr%rxmin)
        sgmtr%rxmax = min((nx - 1)*dx + ox, this%gmtr%rxmax)
        sgmtr%rymin = max(oy, this%gmtr%rymin)
        sgmtr%rymax = min((ny - 1)*dy + oy, this%gmtr%rymax)
        sgmtr%rzmin = max(oz, this%gmtr%rzmin)
        sgmtr%rzmax = min((this%nz - 1)*this%dz + oz, this%gmtr%rzmax)

        if (yn_free_surface) then

            call alloc_array(sgmtr%z_i, [1, 1, 1, 1, -pml, nz + pml + 1])
            call alloc_array(sgmtr%z_h, [1, 1, 1, 1, -pml, nz + pml + 1])
            call alloc_array(sgmtr%dz_i, [1, 1, 1, 1, -pml, nz + pml + 1])
            call alloc_array(sgmtr%dz_h, [1, 1, 1, 1, -pml, nz + pml + 1])
            sgmtr%z_i(1, 1, :) = zz_i
            sgmtr%z_h(1, 1, :) = zz_h
            sgmtr%dz_i(1, 1, :) = dz_i
            sgmtr%dz_h(1, 1, :) = dz_h

            ! Avoid placing source on the free surface as otherwise the resulting amplitude
            ! can be problematic
            where (sgmtr%srcr(:)%z < dz_i(1))
                sgmtr%srcr(:)%z = dz_i(1)
                sgmtr%srcr(:)%amp = sgmtr%srcr(:)%amp*1.2
            end where

        end if

        call sgmtr%prepare_geometry

        sgmtr%dt = dt
        sgmtr%nt = nt
        do i = 1, sgmtr%ns
            if (sgmtr%srcr(i)%mechanism == 'explosion' .or. sgmtr%srcr(i)%mechanism == 'mt') then
                sgmtr%srcr(i)%time_integration = -1
            end if
        end do
        call sgmtr%prepare_stf

        data_dt = this%data_dt
        data_tmax = this%data_tmax
        data_nt = nint(data_tmax/data_dt + 1)

        cc_step_interval = this%cc_step_interval
        verbose = this%verbose
        yn_reconstruct = this%reconstruct

        dir_synthetic = this%dir_synthetic
        dir_adjoint = this%dir_adjoint
        dir_working = this%dir_working

        if (rankid_group == 0) then
            call make_directory(dir_synthetic)
        end if

        snaps = this%snaps
        np = size(snaps)
        if (np > 0) then
            dir_snapshot = this%dir_snapshot
            if (rankid_group == 0) then
                call make_directory(dir_snapshot)
            end if
        end if

        yn_compx = this%compx
        yn_compy = this%compy
        yn_compz = this%compz

        ! If mt inversion is required
        yn_update_medium = this%yn_update_medium
        yn_update_source = this%yn_update_source

        nc_mt = this%nc_mt
        if (norm2(this%mt) > 0) then

            dstf_dt = zeros(nt, sgmtr%ns)

            do i = 1, sgmtr%ns
                sgmtr%srcr(i)%mechanism = 'mt'
                sgmtr%srcr(i)%moment_tensor(1, 1) = this%mt(1)
                sgmtr%srcr(i)%moment_tensor(2, 2) = this%mt(2)
                sgmtr%srcr(i)%moment_tensor(3, 3) = this%mt(3)
                sgmtr%srcr(i)%moment_tensor(1, 2) = this%mt(4)
                sgmtr%srcr(i)%moment_tensor(1, 3) = this%mt(5)
                sgmtr%srcr(i)%moment_tensor(2, 3) = this%mt(6)

                nbeg = nint(sgmtr%srcr(i)%t0/dt) + 1
                nend = nbeg + sgmtr%srcr(i)%nt - 1
                dstf_dt(nbeg:nend, i) = sgmtr%srcr(i)%stf
                dstf_dt(:, i) = deriv(dstf_dt(:, i))

            end do

        end if

        ! Split traces into group processes
        call alloc_array(trace_range, [0, nrank_group - 1, 1, 2])
        call cut(1, sgmtr%nr, nrank_group, trace_range)

        call mpibarrier_group

    end subroutine prepare_modeling

    !
    !> Allocate memory
    !
    subroutine alloc_forward_wavefield

        call alloc_array(stressxx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxxx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvz, [nx1, nx2, ny1, ny2, nz1, nz2])

    end subroutine alloc_forward_wavefield

    subroutine alloc_adjoint_wavefield

        call alloc_array(stressxxr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyyr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszzr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxzr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyzr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxyr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vxr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vyr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vzr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxxxr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxyr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvxr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvyr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxyr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyyr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvxr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvyr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvxr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvyr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvzr, [nx1, nx2, ny1, ny2, nz1, nz2])

    end subroutine alloc_adjoint_wavefield

    !
    !> Check if a point is in a domain decomposition block
    !
    pure function is_in_block(i, j, k) result(y)

        integer, intent(in) :: i, j, k
        logical :: y

        y = (i >= nx1_interior .and. i <= nx2_interior &
            .and. j >= ny1_interior .and. j <= ny2_interior &
            .and. k >= nz1_interior .and. k <= nz2_interior)

    end function

end module
