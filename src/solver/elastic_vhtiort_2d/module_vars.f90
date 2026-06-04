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

module elastic_vhtiort_2d_vars

    use libflit
    use mod_fdcoef
    use mod_su
    use mod_source_receiver, only: source_receiver_geometry
    use mod_utility, only: check_dt_f0
    use mod_source_receiver, only: nkw
    use mod_anisotropy, only: thomsen_to_cij, alkhalifah_tsvankin_to_cij

    use elastic_vhtiort_2d

    implicit none

    real, allocatable, dimension(:, :) :: vx, vz, stressxx, stresszz, stressxz
    real, allocatable, dimension(:, :) :: memory_pdxxx, memory_pdxxz
    real, allocatable, dimension(:, :) :: memory_pdzxz, memory_pdzzz
    real, allocatable, dimension(:, :) :: memory_pdxvx, memory_pdzvx
    real, allocatable, dimension(:, :) :: memory_pdxvz, memory_pdzvz

    real, allocatable, dimension(:, :) :: vxr, vzr, stressxxr, stresszzr, stressxzr
    real, allocatable, dimension(:, :) :: memory_pdxxxr, memory_pdxxzr
    real, allocatable, dimension(:, :) :: memory_pdzxzr, memory_pdzzzr
    real, allocatable, dimension(:, :) :: memory_pdxvxr, memory_pdzvxr
    real, allocatable, dimension(:, :) :: memory_pdxvzr, memory_pdzvzr

    real, allocatable, dimension(:, :) :: snapvx, snapvz
    real, allocatable, dimension(:, :) :: snapvxr, snapvzr

    real, allocatable, dimension(:, :) :: prev_p, prev_stressxx, prev_stresszz, prev_stressxz
    real, allocatable, dimension(:, :) :: strain, strainxx, strainzz, strainxz

    real, allocatable, dimension(:, :) :: prev_pr, prev_stressxxr, prev_stresszzr, prev_stressxzr
    real, allocatable, dimension(:, :) :: strainr, strainxxr, strainzzr, strainxzr

    real, allocatable, dimension(:, :) :: prev_vx, prev_vz
    real, allocatable, dimension(:, :) :: strainvx, strainvz, strainvxr, strainvzr

    real, allocatable, dimension(:, :) :: src_vx, src_vz, rec_vx, rec_vz

    real, allocatable, dimension(:, :) :: energy_src_v, energy_rec_v
    real, allocatable, dimension(:, :) :: energy_src_a, energy_rec_a

    real, allocatable, dimension(:) :: strainxx_udsh, strainxx_udrh, strainxx_udarh
    real, allocatable, dimension(:) :: strainzz_udsh, strainzz_udrh, strainzz_udarh
    real, allocatable, dimension(:) :: strainxz_udsh, strainxz_udrh, strainxz_udarh
    real, allocatable, dimension(:) :: strainxx_lrsh, strainxx_lrrh, strainxx_lrarh
    real, allocatable, dimension(:) :: strainzz_lrsh, strainzz_lrrh, strainzz_lrarh
    real, allocatable, dimension(:) :: strainxz_lrsh, strainxz_lrrh, strainxz_lrarh

    real, allocatable, dimension(:) :: strainvx_udsh, strainvx_udrh
    real, allocatable, dimension(:) :: strainvz_udsh, strainvz_udrh
    real, allocatable, dimension(:) :: strainvx_lrsh, strainvx_lrrh
    real, allocatable, dimension(:) :: strainvz_lrsh, strainvz_lrrh

    real, allocatable, dimension(:) :: p_lrsh, p_lrrh, p_lrarh
    real, allocatable, dimension(:) :: p_udsh, p_udrh, p_udarh

    real, allocatable, dimension(:, :) :: vp, vs, tieps, tidel, tithe, rho, tieta

    real, allocatable, dimension(:, :) :: s11, s13, s33, s55
    real, allocatable, dimension(:, :) :: c11, c13, c15, c33, c35, c55

    real, allocatable, dimension(:, :) :: grad_vp, grad_vs, grad_rho
    real, allocatable, dimension(:, :) :: grad_epsilon, grad_delta, grad_eta
    real, allocatable, dimension(:, :) :: grad_c11, grad_c33, grad_c13, grad_c55
    real, allocatable, dimension(:) :: grad_mt

    real, allocatable, dimension(:) :: snaps

    integer :: np

    real :: idx, idz

    integer :: nx, nz
    real :: dx, dz
    real :: ox, oz

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
    logical :: yn_compz = .true.

    character(len=12) :: aniso_param = 'iso'

    real, allocatable, dimension(:) :: zz_i, zz_h, dz_i, dz_h
    real, allocatable, dimension(:) :: dz_scaling_i, dz_scaling_h

    real :: pmlvp

    integer :: nc_mt

    logical :: yn_update_medium = .true.
    logical :: yn_update_source = .false.
    real, allocatable, dimension(:, :) :: dstf_dt

contains

    !
    !> Map the model on the original regular mesh to the depth-varying mesh
    !
    subroutine map_regular_to_irregular(v, this, rr)

        real, allocatable, dimension(:, :), intent(inout) :: v
        type(wave_solver_elastic_vhtiort_2d), intent(in) :: this
        integer, dimension(1:4), intent(in) :: rr

        real, allocatable, dimension(:, :) :: m
        real, allocatable, dimension(:) :: z
        integer, allocatable, dimension(:) :: r
        integer :: i

        r = zeros(4)
        r(1) = lbound(v, dim=1)
        r(2) = ubound(v, dim=1)
        r(3) = lbound(v, dim=2)
        r(4) = ubound(v, dim=2)

        ! Depth coordinates of regular-mesh array
        z = regspace(r(3) - 1.0, 1.0, r(4) - 1.0)*this%dz

        ! Interpolate to irregular-mesh array
        call alloc_array(m, rr)
        !$omp parallel do private(i)
        do i = max(r(1), rr(1)), min(r(2), rr(2))
            m(i, :) = ginterp(z, v(i, r(3):r(4)), zz_i, 'linear')
        end do
        !$omp end parallel do

        ! Overwrite original
        v = m

    end subroutine

    !
    !> Map the model on the depth-varying mesh to the original regular mesh
    !
    subroutine map_irregular_to_regular(v, this, rr)

        real, allocatable, dimension(:, :), intent(inout) :: v
        type(wave_solver_elastic_vhtiort_2d), intent(in) :: this
        integer, dimension(1:4), intent(in) :: rr

        real, allocatable, dimension(:, :) :: m
        real, allocatable, dimension(:) :: zz
        integer, allocatable, dimension(:) :: r
        integer :: i

        r = zeros(4)
        r(1) = lbound(v, dim=1)
        r(2) = ubound(v, dim=1)
        r(3) = lbound(v, dim=2)
        r(4) = ubound(v, dim=2)

        zz = regspace(rr(3) - 1.0, 1.0, rr(4) - 1.0)*this%dz

        call alloc_array(m, rr)
        !$omp parallel do private(i)
        do i = max(r(1), rr(1)), min(r(2), rr(2))
            m(i, :) = ginterp(zz_i(r(3):r(4)), v(i, r(3):r(4)), zz, 'cubic')
        end do
        !$omp end parallel do

        v = m

    end subroutine

    !
    !> Prepare for simulation
    !
    subroutine prepare_modeling(this)

        type(wave_solver_elastic_vhtiort_2d), intent(in) :: this

        real :: minv, maxv, dtstable, f0clean, temp
        integer :: i, j, l, refine_nz
        real :: rr, dz_max, rayleigh_wavelength
        real, allocatable, dimension(:, :) :: mdispersion, mstability
        integer :: nbeg, nend

        ! Set model
        nx = this%nx
        nz = this%nz
        dx = this%dx
        dz = this%dz
        ox = this%ox
        oz = this%oz

        idx = 1.0d0/dx
        idz = 1.0d0/dz

        pml = this%pml

        call alloc_array(c11, [1, nx, 1, nz])
        call alloc_array(c13, [1, nx, 1, nz])
        call alloc_array(c15, [1, nx, 1, nz])
        call alloc_array(c33, [1, nx, 1, nz])
        call alloc_array(c35, [1, nx, 1, nz])
        call alloc_array(c55, [1, nx, 1, nz])
        call alloc_array(rho, [1, nx, 1, nz])

        aniso_param = this%anisotropy_type
        select case (aniso_param)

            case ('iso')
                vp = transpose(this%vp)
                vs = transpose(this%vs)
                rho = transpose(this%rho)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                c33 = rho*vp**2
                c55 = rho*vs**2
                c11 = c33
                c13 = c33 - 2*c55

            case ('thomsen')
                vp = transpose(this%vp)
                vs = transpose(this%vs)
                tieps = transpose(this%tieps)
                tidel = transpose(this%tidel)
                tithe = transpose(this%tithe)
                rho = transpose(this%rho)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                call thomsen_to_cij(vp, vs, rho, tieps, tidel, tithe, c11, c13, c15, c33, c35, c55)

            case ('a-t')
                vp = transpose(this%vp)
                vs = transpose(this%vs)
                tieps = transpose(this%tieps)
                tieta = transpose(this%tieta)
                tithe = transpose(this%tithe)
                rho = transpose(this%rho)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                call alkhalifah_tsvankin_to_cij(vp, vs, rho, tieps, tieta, tithe, c11, c13, c15, c33, c35, c55)

            case ('cij')
                c11 = transpose(this%c11)
                c13 = transpose(this%c13)
                c33 = transpose(this%c33)
                c55 = transpose(this%c55)
                rho = transpose(this%rho)

                call assert(maxval(abs(c11)) > 0, ' Error: c11 is all zero.')
                call assert(maxval(abs(c33)) > 0, ' Error: c33 is all zero.')
                call assert(maxval(abs(c55)) > 0, ' Error: c55 is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

        end select

        pmlvp = max(maxval(c11(1, :)/rho(1, :)), maxval(c11(nx, :)/rho(nx, :)), &
            maxval(c11(:, 1)/rho(:, 1)), maxval(c11(:, nz)/rho(:, nz)))
        pmlvp = max(pmlvp, max(maxval(c33(1, :)/rho(1, :)), maxval(c33(nx, :)/rho(nx, :)), &
            maxval(c33(:, 1)/rho(:, 1)), maxval(c33(:, nz)/rho(:, nz))))
        pmlvp = sqrt(pmlvp)

        ! Pad models with PML
        call pad_array(c11, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c13, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c33, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c55, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(rho, [pml + 1, pml + 1, pml + 1, pml + 1])

        ! Refine near-surface mesh along z when doing free-surface modeling
        yn_free_surface = this%free_surface
        if (yn_free_surface) then

            ! Approximate Rayleigh wavelength at the free surface
            rayleigh_wavelength = sqrt(mean(c55(:, 1)/rho(:, 1)))/minval(this%gmtr%srcr(:)%f0)
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
            call map_regular_to_irregular(c11, this, [-pml, nx + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c13, this, [-pml, nx + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c15, this, [-pml, nx + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c33, this, [-pml, nx + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c35, this, [-pml, nx + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c55, this, [-pml, nx + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(rho, this, [-pml, nx + pml + 1, -pml, nz + pml + 1])

            call warn(date_time_compact()//' Original nz = '//num2str(this%nz)//', adjusted nz = '//num2str(nz))
            call warn(date_time_compact()//' min(dz) = '//num2str(min(minval(dz_i), minval(dz_h)), '(es)') &
                //', max(dz) = '//num2str(max(maxval(dz_i), maxval(dz_h)), '(es)'))

        end if

        ! Check stability and dispersion
        dt = this%dt
        tmax = this%tmax
        temp = sum(abs(fdcoefs))

        if (yn_free_surface) then
            ! For depth-varying mesh, the stability & dispersion criteria vary in space as well

            mdispersion = zeros(nx, nz)
            mstability = zeros(nx, nz)

            !$omp parallel do private(i, j, minv, maxv) collapse(2)
            do j = 1, nz
                do i = 1, nx

                    minv = sqrt(c55(i, j)/rho(i, j))
                    maxv = max(sqrt(c11(i, j)/rho(i, j)), sqrt(c33(i, j)/rho(i, j)))

                    mstability(i, j) = 1.0/(temp*maxv*sqrt(1.0/dx**2 + 1.0/min(dz_i(j), dz_h(j))**2))
                    ! Empirical relation to ensure sufficient accuracy for surface waves
                    if (sum(dz_i(1:j)) <= 0.5*rayleigh_wavelength) then
                        mdispersion(i, j) = min(0.9*minv/7.0/dx, 0.2*0.9*minv/10.0/max(dz_i(j), dz_h(j)))
                    else
                        mdispersion(i, j) = minv/7.0/max(dx, max(dz_i(j), dz_h(j)))
                    end if

                end do
            end do
            !$omp end parallel do

            dtstable = minval(mstability)
            f0clean = minval(mdispersion)

        else
            ! Fixed mesh

            minv = minval(sqrt(c55/rho))
            maxv = max(maxval(sqrt(c11/rho)), maxval(sqrt(c33/rho)))

            dtstable = 1.0/(temp*maxv*sqrt(1.0/dx**2 + 1.0/dz**2))
            f0clean = minv/max(dx, dz)/7.0

        end if

        call warn(date_time_compact()//' Stable dt = '//num2str(dtstable, '(es)') &
            //' s, clean f0 = '//num2str(f0clean, '(es)')//' Hz')

        call check_dt_f0(dt, dtstable, maxval(this%gmtr%srcr(:)%f0), f0clean)

        nt = nint(tmax/dt + 1)

        ! Prepare geometry
        sgmtr = this%gmtr
        sgmtr%nx = nx
        sgmtr%nz = this%nz
        sgmtr%dx = dx
        sgmtr%dz = this%dz
        sgmtr%ox = ox
        sgmtr%oz = oz
        sgmtr%xmin = max(ox, this%gmtr%xmin)
        sgmtr%xmax = min((nx - 1)*dx + ox, this%gmtr%xmax)
        sgmtr%zmin = max(oz, this%gmtr%zmin)
        sgmtr%zmax = min((this%nz - 1)*this%dz + oz, this%gmtr%zmax)
        sgmtr%sxmin = max(ox, this%gmtr%sxmin)
        sgmtr%sxmax = min((nx - 1)*dx + ox, this%gmtr%sxmax)
        sgmtr%szmin = max(oz, this%gmtr%szmin)
        sgmtr%szmax = min((this%nz - 1)*this%dz + oz, this%gmtr%szmax)
        sgmtr%rxmin = max(ox, this%gmtr%rxmin)
        sgmtr%rxmax = min((nx - 1)*dx + ox, this%gmtr%rxmax)
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

        call make_directory(dir_synthetic)

        snaps = this%snaps
        np = size(snaps)
        if (np > 0) then
            dir_snapshot = this%dir_snapshot
            call make_directory(dir_snapshot)
        end if

        yn_compx = this%compx
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
                sgmtr%srcr(i)%moment_tensor(3, 3) = this%mt(3)
                sgmtr%srcr(i)%moment_tensor(1, 3) = this%mt(5)

                nbeg = nint(sgmtr%srcr(i)%t0/dt) + 1
                nend = nbeg + sgmtr%srcr(i)%nt - 1
                dstf_dt(nbeg:nend, i) = sgmtr%srcr(i)%stf
                dstf_dt(:, i) = deriv(dstf_dt(:, i))

            end do

        end if

    end subroutine prepare_modeling

    !
    !> Allocate memory for forward wavefields
    !
    subroutine alloc_forward_wavefield

        call alloc_array(vx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(vz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(stressxx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(stresszz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(stressxz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxvx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzvx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxvz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzvz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxxx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzxz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxxz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzzz, [1, nx, 1, nz], pad=pml + fdhalf)

    end subroutine alloc_forward_wavefield

    !
    !> Allocate memory for adjoint wavefields
    !
    subroutine alloc_adjoint_wavefield

        call alloc_array(vxr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(vzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(stressxxr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(stresszzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(stressxzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxvxr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzvxr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxvzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzvzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxxxr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzxzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxxzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzzzr, [1, nx, 1, nz], pad=pml + fdhalf)

    end subroutine alloc_adjoint_wavefield

end module
