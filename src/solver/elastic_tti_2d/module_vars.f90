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

module elastic_tti_2d_vars

    use libflit
    use mod_fdcoef
    use mod_su
    use mod_source_receiver, only: source_receiver_geometry
    use mod_utility, only: check_dt_f0
    use mod_source_receiver, only: nkw
    use mod_anisotropy, only: thomsen_to_cij, alkhalifah_tsvankin_to_cij, min_max_phase_velocity_2d

    use elastic_tti_2d

    implicit none

    real, allocatable, dimension(:, :) :: vx_hxiz, vx_ixhz
    real, allocatable, dimension(:, :) :: vz_hxiz, vz_ixhz
    real, allocatable, dimension(:, :) :: stressxx_hxhz, stressxx_ixiz
    real, allocatable, dimension(:, :) :: stresszz_hxhz, stresszz_ixiz
    real, allocatable, dimension(:, :) :: stressxz_hxhz, stressxz_ixiz
    real, allocatable, dimension(:, :) :: memory_pdxvx_hxhz, memory_pdxvz_hxhz
    real, allocatable, dimension(:, :) :: memory_pdxvx_ixiz, memory_pdxvz_ixiz
    real, allocatable, dimension(:, :) :: memory_pdxxx_hxiz, memory_pdxxz_hxiz
    real, allocatable, dimension(:, :) :: memory_pdxxx_ixhz, memory_pdxxz_ixhz
    real, allocatable, dimension(:, :) :: memory_pdzvx_hxhz, memory_pdzvz_hxhz
    real, allocatable, dimension(:, :) :: memory_pdzvx_ixiz, memory_pdzvz_ixiz
    real, allocatable, dimension(:, :) :: memory_pdzxx_hxiz, memory_pdzxz_hxiz, memory_pdzzz_hxiz
    real, allocatable, dimension(:, :) :: memory_pdzxx_ixhz, memory_pdzxz_ixhz, memory_pdzzz_ixhz

    real, allocatable, dimension(:, :) :: vxr_hxiz, vxr_ixhz
    real, allocatable, dimension(:, :) :: vzr_hxiz, vzr_ixhz
    real, allocatable, dimension(:, :) :: stressxxr_hxhz, stressxxr_ixiz
    real, allocatable, dimension(:, :) :: stresszzr_hxhz, stresszzr_ixiz
    real, allocatable, dimension(:, :) :: stressxzr_hxhz, stressxzr_ixiz
    real, allocatable, dimension(:, :) :: memory_pdxvxr_hxhz, memory_pdxvzr_hxhz
    real, allocatable, dimension(:, :) :: memory_pdxvxr_ixiz, memory_pdxvzr_ixiz
    real, allocatable, dimension(:, :) :: memory_pdxxxr_hxiz, memory_pdxxzr_hxiz
    real, allocatable, dimension(:, :) :: memory_pdxxxr_ixhz, memory_pdxxzr_ixhz
    real, allocatable, dimension(:, :) :: memory_pdzvxr_hxhz, memory_pdzvzr_hxhz
    real, allocatable, dimension(:, :) :: memory_pdzvxr_ixiz, memory_pdzvzr_ixiz
    real, allocatable, dimension(:, :) :: memory_pdzxxr_hxiz, memory_pdzxzr_hxiz, memory_pdzzzr_hxiz
    real, allocatable, dimension(:, :) :: memory_pdzxxr_ixhz, memory_pdzxzr_ixhz, memory_pdzzzr_ixhz

    real, allocatable, dimension(:, :) :: prev_vx_hxiz, prev_vx_ixhz
    real, allocatable, dimension(:, :) :: prev_vz_hxiz, prev_vz_ixhz
    real, allocatable, dimension(:, :) :: prev_stressxx_hxhz, prev_stressxx_ixiz
    real, allocatable, dimension(:, :) :: prev_stresszz_hxhz, prev_stresszz_ixiz
    real, allocatable, dimension(:, :) :: prev_stressxz_hxhz, prev_stressxz_ixiz

    real, allocatable, dimension(:, :) :: strain, strainxx, strainzz, strainxz
    real, allocatable, dimension(:, :) :: strainr, strainxxr, strainzzr, strainxzr
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

    real, allocatable, dimension(:, :) :: vp, vs, tieps, tidel, tithe, rho
    real, allocatable, dimension(:, :) :: tieta

    real, allocatable, dimension(:, :) :: s11, s13, s15, s33, s35, s55
    real, allocatable, dimension(:, :) :: c11, c13, c15, c33, c35, c55

    real, allocatable, dimension(:, :) :: grad_c11, grad_c13, grad_c15, grad_c33, grad_c35, grad_c55
    real, allocatable, dimension(:, :) :: grad_vp, grad_vs, grad_epsilon, grad_delta, grad_rho
    real, allocatable, dimension(:, :) :: grad_eta

    real, allocatable, dimension(:, :) :: snapvx, snapvz

    integer :: nx1, nx2
    integer :: nz1, nz2

    real, allocatable, dimension(:) :: snaps

    integer :: np

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

    character(len=24) :: kernel_v = 'full'
    character(len=24) :: kernel_a = 'full'

    character(len=1024) :: dir_synthetic, dir_snapshot, dir_working, dir_adjoint

    logical :: yn_compx = .true.
    logical :: yn_compz = .true.

    character(len=12) :: aniso_param = 'iso'

    real, allocatable, dimension(:, :) :: zz_i, zz_h, dz_i, dz_h
    real, allocatable, dimension(:, :) :: dz_scaling_i, dz_scaling_h

    real, allocatable, dimension(:) :: topo_i, slopex_i
    real, allocatable, dimension(:) :: topo_h, slopex_h

    real :: depth_max, eta_max

    real, allocatable, dimension(:) :: eta_zz_i, eta_zz_h
    real, allocatable, dimension(:) :: eta_dz_i, eta_dz_h
    real, allocatable, dimension(:) :: eta_dz_scaling_i, eta_dz_scaling_h

    real, allocatable, dimension(:) :: topox

    real :: topo_max
    real :: pmlvp

    integer :: nc_mt

    logical :: yn_update_medium = .true.
    logical :: yn_update_source = .false.
    real, allocatable, dimension(:, :) :: dstf_dt
    real, allocatable, dimension(:) :: grad_mt

contains

    !
    !> Map the model on the original regular mesh to the depth-varying mesh
    !
    subroutine map_regular_to_irregular(v, this, rr)

        real, allocatable, dimension(:, :), intent(inout) :: v
        type(wave_solver_elastic_tti_2d), intent(in) :: this
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
        z = regspace(r(3) - 1.0, 1.0, r(4) - 1.0)*this%dz - maxval(topo_i)

        ! Interpolate to irregular-mesh array
        call alloc_array(m, rr)

        !$omp parallel do private(i)
        do i = max(r(1), rr(1)), min(r(2), rr(2))

            m(i, :) = ginterp(z, v(i, r(3):r(4)), zz_i(i, :), 'linear')

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
        type(wave_solver_elastic_tti_2d), intent(in) :: this
        integer, dimension(1:4), intent(in) :: rr

        real, allocatable, dimension(:, :) :: m
        real, allocatable, dimension(:) :: zz
        integer, allocatable, dimension(:) :: r
        integer :: i, j
        real :: max_topo

        r = zeros(4)
        r(1) = lbound(v, dim=1)
        r(2) = ubound(v, dim=1)
        r(3) = lbound(v, dim=2)
        r(4) = ubound(v, dim=2)

        zz = regspace(rr(3) - 1.0, 1.0, rr(4) - 1.0)*this%dz

        call alloc_array(m, rr)

        max_topo = maxval(topo_i)

        !$omp parallel do private(i, j)
        do i = max(r(1), rr(1)), min(r(2), rr(2))

            m(i, :) = ginterp(zz_i(i, r(3):r(4)), v(i, r(3):r(4)), zz - max_topo, 'cubic')

            ! Mask out points above the free surface
            do j = 1, size(zz)
                if (zz(j) < max_topo - topo_i(i)) then
                    m(i, j) = 0
                end if
            end do

        end do
        !$omp end parallel do

        v = m

    end subroutine

    !
    !> Prepare modeling-related arrays and parameters
    !
    subroutine prepare_modeling(this)

        type(wave_solver_elastic_tti_2d), intent(in) :: this

        real :: minv, maxv, dtstable, f0clean, temp
        integer :: i, j, l, refine_nz, ntp
        real :: rr, dz_max, rayleigh_wavelength
        real :: tmpmin, tmpmax
        real, allocatable, dimension(:, :) :: mdispersion, mstability, topo
        real, allocatable, dimension(:) :: stp, rtp
        real :: dz0
        integer :: nbeg, nend

        nx = this%nx
        nz = this%nz
        dx = this%dx
        dz = this%dz
        ox = this%ox
        oz = this%oz

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
                c15 = transpose(this%c15)
                c33 = transpose(this%c33)
                c35 = transpose(this%c35)
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
        call pad_array(c15, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c33, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c35, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c55, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(rho, [pml + 1, pml + 1, pml + 1, pml + 1])

        ! Refine near-surface mesh along z when doing free-surface modeling
        yn_free_surface = this%free_surface
        if (yn_free_surface) then

            ! =================================================================
            ! Setup topography
            if (this%file_topo /= '') then

                ntp = count_nonempty_lines(this%file_topo)
                topo = load(this%file_topo, ntp, 2, ascii=.true.)

                topox = regspace(ox - (pml + 1)*dx, dx, ox + (nx + pml)*dx)

                call alloc_array(topo_i, [-pml, nx + pml + 1], &
                    source=ginterp(topo(:, 1), topo(:, 2), topox, method=this%topo_interp))
                call alloc_array(slopex_i, [-pml, nx + pml + 1], &
                    source=deriv(topo_i, method='center')/dx)

                call alloc_array(topo_h, [-pml, nx + pml + 1], &
                    source=ginterp(topox, topo_i, topox - 0.5*dx, method=this%topo_interp))
                call alloc_array(slopex_h, [-pml, nx + pml + 1], &
                    source=deriv(topo_h, method='center')/dx)

                call warn(date_time_compact()//' Number of topography control points = '//num2str(ntp))
                call warn(date_time_compact()//' Max elevation above z=0 = '//num2str(maxval(topo_i), '(es)'))
                call warn(date_time_compact()//' Slope range = '//num2str(minval(slopex_i), '(es)') &
                    //' ~ '//num2str(maxval(slopex_i), '(es)'))

            else

                call alloc_array(topo_i, [-pml, nx + pml + 1])
                call alloc_array(slopex_i, [-pml, nx + pml + 1])

                call alloc_array(topo_h, [-pml, nx + pml + 1])
                call alloc_array(slopex_h, [-pml, nx + pml + 1])

            end if

            ! =================================================================
            ! Setup mesh

            ! Maximum depth in physical and virtual domains
            topo_max = maxval(topo_i)
            depth_max = (this%nz - 1)*this%dz - topo_max
            eta_max = depth_max

            ! Check stability of mesh
            if (any((topo_i + depth_max)/eta_max < max(1.0, abs(slopex_i)))) then
                eta_max = 0.975*minval((topo_i + depth_max)/max(1.0, abs(slopex_i)))
                call warn(date_time_compact()//' Set eta_max = '//num2str(eta_max, '(es)')//' to ensure stability. ')
            end if

            dz = eta_max/(nz - 1)

            ! Approximate Rayleigh wavelength at the free surface
            rayleigh_wavelength = sqrt(mean(c55(:, 1)/rho(:, 1)))/minval(this%gmtr%srcr(:)%f0)
            refine_nz = ceiling(0.5*rayleigh_wavelength/dz)

            ! Compute integer grid point positions
            ! Initial dz = dz/dz_refine, and it gradually increases after refine_nz, by 2.5% per grid
            ! 5% can cause very weak artificial reflections, while 1% is just too costly.
            rr = this%free_surface_dz_refine
            dz_max = min(this%dz_max*eta_max/(depth_max + topo_max), dz/rr)

            ! First compute a standard depth varying 1D mesh
            eta_zz_i = regspace(0.0, dz_max, (refine_nz - 1)*dz)
            l = size(eta_zz_i)
            do while (eta_zz_i(l) < eta_max)
                dz_max = min(dz_max*1.025, this%dz_max*eta_max/(depth_max + topo_max))
                eta_zz_i = [eta_zz_i, eta_zz_i(l) + dz_max]
                l = l + 1
            end do

            ! New nz is the length of eta_zz_i
            nz = size(eta_zz_i)

            ! Half integer grid points are interpolated
            eta_zz_h = zeros(nz)
            eta_zz_h(1) = -0.5*dz/rr
            do i = 1, nz - 1
                eta_zz_h(i + 1) = 0.5*(eta_zz_i(i) + eta_zz_i(i + 1))
            end do

            call pad_array(eta_zz_i, [pml + 1, pml + 1])
            call pad_array(eta_zz_h, [pml + 1, pml + 1])

            ! Note that after extension, max(eta_zz_i) may >, =, or < eta_max
            ! To map the irregular wavefield to the regular mesh, the
            do i = 1, pml + 1
                eta_zz_i(1 - i) = eta_zz_i(1) - i*dz/rr
                eta_zz_i(nz + i) = eta_zz_i(nz) + i*dz_max
                eta_zz_h(1 - i) = eta_zz_h(1) - i*dz/rr
                eta_zz_h(nz + i) = eta_zz_h(nz) + i*dz_max
            end do

            call alloc_array(eta_dz_i, [1, nz], pad=pml + 1, source=deriv(eta_zz_i, method='center'))
            call alloc_array(eta_dz_h, [1, nz], pad=pml + 1, source=deriv(eta_zz_h, method='center'))

            call alloc_array(eta_dz_scaling_i, [1, nz], pad=pml + 1, source=dz/eta_dz_i)
            call alloc_array(eta_dz_scaling_h, [1, nz], pad=pml + 1, source=dz/eta_dz_h)

            ! Then scale the standard irregular mesh to the entire domain
            call alloc_array(zz_i, [1, nx, 1, nz], pad=pml + 1)
            call alloc_array(zz_h, [1, nx, 1, nz], pad=pml + 1)

            ! Then set the mesh above the bottom PML interface
            !$omp parallel do private(i, j)
            do j = -pml, nx + pml + 1

                zz_i(j, 1:nz) = rescale(eta_zz_i(1:nz), [-topo_i(j), depth_max])

            end do
            !$omp end parallel do

            dz_max = mean(zz_i(:, nz) - zz_i(:, nz - 1))

            ! Then set the mesh below the bottom PML interface
            !$omp parallel do private(i, j)
            do j = -pml, nx + pml + 1

                do i = 1, pml + 1
                    zz_i(j, nz + i) = zz_i(j, nz) + i*dz_max
                end do

                zz_h(j, 1) = -topo_i(j) - 0.5*(zz_i(j, 2) - zz_i(j, 1))
                do i = 1, nz + pml
                    zz_h(j, i + 1) = 0.5*(zz_i(j, i) + zz_i(j, i + 1))
                end do

                do i = 1, pml + 1
                    zz_i(j, 1 - i) = zz_i(j, 1) - i*(zz_i(j, 2) - zz_i(j, 1))
                    zz_h(j, 1 - i) = zz_h(j, 1) - i*(zz_h(j, 2) - zz_h(j, 1))
                end do

            end do
            !$omp end parallel do

            ! Depth varying grid size
            call alloc_array(dz_i, [1, nx, 1, nz], pad=pml + 1, source=deriv(zz_i, method='center', dim=2))
            call alloc_array(dz_h, [1, nx, 1, nz], pad=pml + 1, source=deriv(zz_h, method='center', dim=2))

            ! =================================================================
            ! Map regular-mesh models to irregular mesh
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

            ! Clip virtual mesh to make sure it does not exceed eta_max
            eta_zz_i = clip(eta_zz_i, -pml*dz, eta_max)
            eta_zz_h = clip(eta_zz_h, -pml*dz, eta_max)

        end if

        ! Domain decomposition -- must be placed here since mesh adjusting changed nz
        nx1 = -pml + 1
        nx2 = nx + pml
        nz1 = -pml + 1
        nz2 = nz + pml

        ! Check stability and dispersion
        dt = this%dt
        tmax = this%tmax
        temp = sum(abs(fdcoefs))

        if (yn_free_surface) then
            ! Irregular mesh

            call alloc_array(mdispersion, [1, nx, 1, nz], pad=pml + 1)
            call alloc_array(mstability, [1, nx, 1, nz], pad=pml + 1)
            mdispersion = float_huge
            mstability = float_huge

            !$omp parallel do private(i, j, minv, maxv, dz0) collapse(2)
            do j = 1, nz + pml + 1
                do i = -pml + 1, nx + pml

                    select case (aniso_param)
                        case ('iso')
                            call min_max_phase_velocity_2d(c11(i, j), c13(i, j), c15(i, j), c33(i, j), c35(i, j), c55(i, j), &
                                rho(i, j), 1, maxv, minv)
                        case default
                            ! Get min/max phase velocities for each point; for each point, the velocities are computed at
                            ! 36 polar angles for ensure accuracy.
                            call min_max_phase_velocity_2d(c11(i, j), c13(i, j), c15(i, j), c33(i, j), c35(i, j), c55(i, j), &
                                rho(i, j), 36, maxv, minv)
                    end select

                    dz0 = sqrt((1.0/dx + max(0.0, eta_max - eta_zz_i(j))/eta_max*slopex_i(i)*1.0/eta_dz_i(j))**2 &
                        + (1.0/eta_dz_i(j))**2)

                    mstability(i, j) = 1.0/(temp*maxv*dz0)
                    ! Empirical relation to ensure sufficient accuracy for surface waves
                    if (sum(dz_i(i, 1:j)) <= 0.5*rayleigh_wavelength) then
                        mdispersion(i, j) = min(0.9*minv/7.0/dx, 0.2*0.9*minv/10.0/max(dz_i(i, j), dz_h(i, j)))
                    else
                        mdispersion(i, j) = minv/7.0/max(dx, max(dz_i(i, j), dz_h(i, j)))
                    end if

                end do
            end do
            !$omp end parallel do

            dtstable = minval(mstability)
            f0clean = minval(mdispersion)

        else
            ! Regular mesh

            select case (aniso_param)
                case ('iso')
                    minv = minval(vs)
                    maxv = maxval(vp)
                case ('thomsen')
                    minv = minval(vs)
                    maxv = maxval(vp*sqrt(1.0 + 2.0*tieps))
                case ('a-t')
                    minv = minval(vs)
                    maxv = maxval(vp)
                case ('cij')
                    minv = float_huge
                    maxv = -float_huge
                    !$omp parallel do private(i, j) collapse(2) reduction(min:minv) reduction(max:maxv)
                    do j = 1, nz
                        do i = 1, nx

                            ! Get min/max phase velocities for each point; for each point, the velocities are computed at
                            ! 36 polar angles for ensure accuracy.
                            call min_max_phase_velocity_2d(c11(i, j), c13(i, j), c15(i, j), c33(i, j), c35(i, j), c55(i, j), &
                                rho(i, j), 36, tmpmax, tmpmin)

                            minv = min(minv, tmpmin)
                            maxv = max(maxv, tmpmax)

                        end do
                    end do
                    !$omp end parallel do
            end select

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

            ! If measure from the free surface then + topo (here, topo has already been sign-corrected)
            ! e.g., sz = 400, topo = 100, then real depth = 500
            topox = regspace(0.0, dx, (nx + 2*pml + 1)*dx) - (pml + 1)*dx

            stp = ginterp(topox, zz_i(:, 1), sgmtr%srcr(:)%x, this%topo_interp)
            rtp = ginterp(topox, zz_i(:, 1), sgmtr%recr(:)%x, this%topo_interp)

            if (this%measure_source_depth_from_surface) then
                sgmtr%srcr(:)%z = sgmtr%srcr(:)%z + stp
            end if
            if (this%measure_receiver_depth_from_surface) then
                sgmtr%recr(:)%z = sgmtr%recr(:)%z + rtp
            end if

            ! Map real location to virtual mesh
            sgmtr%srcr(:)%z = (sgmtr%srcr(:)%z - stp)/(-stp + depth_max)*eta_max
            sgmtr%recr(:)%z = (sgmtr%recr(:)%z - rtp)/(-rtp + depth_max)*eta_max

            ! For FSG with topographic free surface, if a source is at the surface (depth = 0),
            ! then must adjust it to the first grid below the free surface; because
            ! unlike in SSG + horizontal free surface, here we cannot set sigmaxz(j = 1) = sigmazz(j = 1) = 0
            where (sgmtr%srcr(:)%z < eta_dz_i(1))
                sgmtr%srcr(:)%z = eta_dz_i(1)
                sgmtr%srcr(:)%amp = sgmtr%srcr(:)%amp*1.2
            end where

            call alloc_array(sgmtr%z_i, [-pml, nx + pml + 1, 1, 1, -pml, nz + pml + 1])
            call alloc_array(sgmtr%z_h, [-pml, nx + pml + 1, 1, 1, -pml, nz + pml + 1])
            call alloc_array(sgmtr%dz_i, [-pml, nx + pml + 1, 1, 1, -pml, nz + pml + 1])
            call alloc_array(sgmtr%dz_h, [-pml, nx + pml + 1, 1, 1, -pml, nz + pml + 1])
            !$omp parallel do private(i)
            do i = -pml, nx + pml + 1
                sgmtr%z_i(i, 1, :) = eta_zz_i
                sgmtr%z_h(i, 1, :) = eta_zz_h
                sgmtr%dz_i(i, 1, :) = eta_dz_i
                sgmtr%dz_h(i, 1, :) = eta_dz_h
            end do
            !$omp end parallel do

        end if

        call sgmtr%prepare_geometry

        ! Prepare source wavelet
        sgmtr%dt = dt
        sgmtr%nt = nt
        do i = 1, sgmtr%ns
            if (sgmtr%srcr(i)%mechanism == 'explosion' .or. sgmtr%srcr(i)%mechanism == 'mt') then
                sgmtr%srcr(i)%time_integration = -1
            end if
        end do
        call sgmtr%prepare_stf

        ! Other parameters
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

        ! Save generated mesh
        if (yn_free_surface .and. this%save_mesh) then

            open(3, file=tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_mesh.txt')
            do i = -pml + 1, nx + pml
                do j = 1, nz + pml
                    write(3, '(4es)') (i - 1)*dx + ox, zz_i(i, j), c11(i, j), mstability(i, j)
                end do
            end do
            close(3)

        end if

        !        ! Save source interpolation function for debugging
        !        open(3, file=tidy(dir_synthetic)//'/source_interp.txt')
        !        do i = -nkw, nkw
        !            do j = -nkw, nkw
        !                write(3, *) i*1.0, eta_zz_i(sgmtr%srcr(1)%gz + j), &
            !                   sgmtr%srcr(1)%interp_ix(i)*sgmtr%srcr(1)%interp_iz(j)
        !            end do
        !        end do
        !        do i = -nkw, nkw
        !            do j = -nkw, nkw
        !                write(3, *) i*1.0 + 10.0, eta_zz_h(sgmtr%srcr(1)%hz + j), &
            !                   sgmtr%srcr(1)%interp_hx(i)*sgmtr%srcr(1)%interp_hz(j)
        !            end do
        !        end do
        !        close(3)

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

    end subroutine

    !
    !> Allocate memory for forward wavefields
    !
    subroutine alloc_forward_wavefield

        call alloc_array(vx_hxiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(vx_ixhz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(vz_hxiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(vz_ixhz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxx_ixiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxx_hxhz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszz_ixiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszz_hxhz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxz_ixiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxz_hxhz, [nx1, nx2, nz1, nz2], pad=fdhalf)

        if (yn_free_surface) then

            call alloc_array(memory_pdxvx_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvz_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvx_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvz_ixiz, [nx1, nx2, nz1, nz2])

            call alloc_array(memory_pdzvx_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvz_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvx_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvz_ixiz, [nx1, nx2, nz1, nz2])

            call alloc_array(memory_pdxxx_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxz_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxx_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxz_ixhz, [nx1, nx2, nz1, nz2])

            call alloc_array(memory_pdzxx_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxz_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzzz_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxx_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxz_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzzz_ixhz, [nx1, nx2, nz1, nz2])

        else
            call alloc_array(memory_pdxvx_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvz_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvx_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvz_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxx_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxz_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxx_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxz_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvx_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvz_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvx_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvz_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxz_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzzz_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxz_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzzz_ixhz, [nx1, nx2, nz1, nz2])
        end if

    end subroutine

    !
    !> Allocate memory for adjoint wavefields
    !
    subroutine alloc_adjoint_wavefield

        call alloc_array(vxr_hxiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(vxr_ixhz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(vzr_hxiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(vzr_ixhz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxxr_ixiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxxr_hxhz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszzr_ixiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszzr_hxhz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxzr_ixiz, [nx1, nx2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxzr_hxhz, [nx1, nx2, nz1, nz2], pad=fdhalf)

        if (yn_free_surface) then

            call alloc_array(memory_pdxvxr_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvzr_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvxr_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvzr_ixiz, [nx1, nx2, nz1, nz2])

            call alloc_array(memory_pdzvxr_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvzr_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvxr_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvzr_ixiz, [nx1, nx2, nz1, nz2])

            call alloc_array(memory_pdxxxr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxzr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxxr_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxzr_ixhz, [nx1, nx2, nz1, nz2])

            call alloc_array(memory_pdzxxr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxzr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzzzr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxxr_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxzr_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzzzr_ixhz, [nx1, nx2, nz1, nz2])

        else

            call alloc_array(memory_pdxvxr_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvzr_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvxr_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxvzr_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxxr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxzr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxxr_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdxxzr_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvxr_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvzr_hxhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvxr_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzvzr_ixiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxzr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzzzr_hxiz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzxzr_ixhz, [nx1, nx2, nz1, nz2])
            call alloc_array(memory_pdzzzr_ixhz, [nx1, nx2, nz1, nz2])

        end if

    end subroutine

end module
