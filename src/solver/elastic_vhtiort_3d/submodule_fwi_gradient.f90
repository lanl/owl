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

submodule(elastic_vhtiort_3d) elastic_vhtiort_3d_fwi_gradient

    use libflit
    use elastic_vhtiort_3d_vars
    use elastic_vhtiort_3d_boundary_saving
    use elastic_vhtiort_3d_wavefield

    implicit none

#include '../../lib/macro_thomsen_3d.f90'
#include '../../lib/macro_alkhalifah_tsvankin_3d.f90'

contains

    module subroutine compute_adjoint(this)

        class(wave_solver_elastic_vhtiort_3d), intent(inout) :: this

        integer :: l, ir, irx, iry, irz, rgx, rgy, rgz
        type(grid3) :: grd
        integer :: i, j, k, t
        real, allocatable, dimension(:, :, :) :: grad
        real :: wmin1, wmin2, wmin3
        real :: wmax1, wmax2, wmax3
        logical :: wnan

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield
        call alloc_adjoint_wavefield

        yn_energy_precond = this%energy_precond

        ! Adjoint source
        if (yn_compx) then
            call this%seis_vxr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_x.su')
            call this%seis_vxr%zero_foreign_rank_traces_group
            call this%seis_vxr%resamp(nnt=nt, ddt=dt)
            call this%seis_vxr%collect_group
        end if

        if (yn_compy) then
            call this%seis_vyr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_y.su')
            call this%seis_vyr%zero_foreign_rank_traces_group
            call this%seis_vyr%resamp(nnt=nt, ddt=dt)
            call this%seis_vyr%collect_group
        end if

        if (yn_compz) then
            call this%seis_vzr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_z.su')
            call this%seis_vzr%zero_foreign_rank_traces_group
            call this%seis_vzr%resamp(nnt=nt, ddt=dt)
            call this%seis_vzr%collect_group
        end if

        ! Elastic compliances
        call alloc_array(s11, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        call alloc_array(s12, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        call alloc_array(s13, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        call alloc_array(s22, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        call alloc_array(s23, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        call alloc_array(s33, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        call alloc_array(s44, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        call alloc_array(s55, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        call alloc_array(s66, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1)
        s11 = (c23**2 - c22*c33)/(c13**2*c22 &
            - 2*c12*c13*c23 + c11*c23**2 + c12**2*c33 - c11*c22*c33)
        s12 = (-c13*c23 + c12*c33)/(c13**2*c22 &
            - 2*c12*c13*c23 + c12**2*c33 + c11*(c23**2 - c22*c33))
        s13 = (c13*c22 - c12*c23)/(c13**2*c22 &
            - 2*c12*c13*c23 + c11*c23**2 + c12**2*c33 - c11*c22*c33)
        s22 = (c13**2 - c11*c33)/(c13**2*c22 &
            - 2*c12*c13*c23 + c11*c23**2 + c12**2*c33 - c11*c22*c33)
        s23 = (-c12*c13 + c11*c23)/(c13**2*c22 &
            - 2*c12*c13*c23 + c12**2*c33 + c11*(c23**2 - c22*c33))
        s33 = (c12**2 - c11*c22)/(c13**2*c22 &
            - 2*c12*c13*c23 + c11*c23**2 + c12**2*c33 - c11*c22*c33)
        where (c44 /= 0)
            s44 = 1.0/c44
        end where
        where (c55 /= 0)
            s55 = 1.0/c55
        end where
        where (c66 /= 0)
            s66 = 1.0/c66
        end where

        ! Allocate memory for arrays
        grad_c11 = zeros(nx, ny, nz)
        grad_c12 = zeros(nx, ny, nz)
        grad_c13 = zeros(nx, ny, nz)
        grad_c22 = zeros(nx, ny, nz)
        grad_c23 = zeros(nx, ny, nz)
        grad_c33 = zeros(nx, ny, nz)
        grad_c44 = zeros(nx, ny, nz)
        grad_c55 = zeros(nx, ny, nz)
        grad_c66 = zeros(nx, ny, nz)
        grad_rho = zeros(nx, ny, nz)
        grad_mt = zeros(nc_mt)

        energy_src_v = zeros(nx, ny, nz)
        energy_rec_v = zeros(nx, ny, nz)
        energy_src_a = zeros(nx, ny, nz)
        energy_rec_a = zeros(nx, ny, nz)

        call alloc_array(strainxx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainyy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainzz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainyz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainxz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainxy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainxxr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainyyr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainzzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainyzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainxzr, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(strainxyr, [nx1, nx2, ny1, ny2, nz1, nz2])

        call alloc_array(src_vx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(rec_vx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(src_vy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(rec_vy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(src_vz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(rec_vz, [nx1, nx2, ny1, ny2, nz1, nz2])

        call alloc_array(prev_stressxx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stressyy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stresszz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stressyz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stressxz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stressxy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(snapvx, [1, nx, 1, ny, 1, nz], pad=pml)
        call alloc_array(snapvy, [1, nx, 1, ny, 1, nz], pad=pml)
        call alloc_array(snapvz, [1, nx, 1, ny, 1, nz], pad=pml)

        ! Prepare boundary saving
        call prepare_boundary_saving
        call open_boundary_saving

        ! Compute gradients by cross-correlation
        l = np
        do t = nt, 1, -1

            if (yn_update_medium .and. t >= sgmtr%srcr(1)%hnt) then

                ! -------------- Wavefield reconstruction ----------------------
                ! Store previous stress wavefields for cross-correlation
                prev_stressxx = stressxx
                prev_stressyy = stressyy
                prev_stresszz = stresszz
                prev_stressyz = stressyz
                prev_stressxz = stressxz
                prev_stressxy = stressxy
                prev_vx = vx
                prev_vy = vy
                prev_vz = vz

                ! -------------- Forward wavefield reconstruction -----------------------
                if (yn_free_surface) then
                    call update_wavefield_free_surface(-dt, &
                        stressxx, stressyy, stresszz, &
                        stressyz, stressxz, stressxy, &
                        vx, vy, vz, &
                        memory_pdxxx, memory_pdyxy, memory_pdzxz, &
                        memory_pdxxy, memory_pdyyy, memory_pdzyz, &
                        memory_pdxxz, memory_pdyyz, memory_pdzzz, &
                        memory_pdxvx, memory_pdyvx, memory_pdzvx, &
                        memory_pdxvy, memory_pdyvy, memory_pdzvy, &
                        memory_pdxvz, memory_pdyvz, memory_pdzvz)
                else
                    call update_wavefield(-dt, &
                        stressxx, stressyy, stresszz, &
                        stressyz, stressxz, stressxy, &
                        vx, vy, vz, &
                        memory_pdxxx, memory_pdyxy, memory_pdzxz, &
                        memory_pdxxy, memory_pdyyy, memory_pdzyz, &
                        memory_pdxxz, memory_pdyyz, memory_pdzzz, &
                        memory_pdxvx, memory_pdyvx, memory_pdzvx, &
                        memory_pdxvy, memory_pdyvy, memory_pdzvy, &
                        memory_pdxvz, memory_pdyvz, memory_pdzvz)
                end if

                ! Read final step wavefield
                if (t == nt) then
                    call input_final_step_wavefield
                end if

                ! Read boundary wavefield
                call inject_boundary_wavefield(t)

                ! Record wavefield snapshot if necessary
                if (np /= 0 .and. l >= 1) then
                    if (t - 1 == nint(snaps(l)/dt)) then

                        snapvx = 0
                        snapvy = 0
                        snapvz = 0

                        call commute_array_group(vx, fdhalf)
                        call commute_array_group(vy, fdhalf)
                        call commute_array_group(vz, fdhalf)

                        !$omp parallel do private(i, j, k) collapse(3)
                        do k = nz1, nz2
                            do j = ny1, ny2
                                do i = nx1, nx2
                                    snapvx(i, j, k) = 0.5*sum(vx(i:i + 1, j, k))
                                    snapvy(i, j, k) = 0.5*sum(vy(i, j:j + 1, k))
                                    snapvz(i, j, k) = 0.5*sum(vz(i, j, k:k + 1))
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        call reduce_array_group(snapvx)
                        call reduce_array_group(snapvy)
                        call reduce_array_group(snapvz)

                        ! Output
                        if (rankid_group == 0) then
                            call output_array(snapvx(1:nx, 1:ny, 1:nz), tidy(dir_working)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_x_' &
                                //num2str(l)//'.bin', store=321)
                            call output_array(snapvy(1:nx, 1:ny, 1:nz), tidy(dir_working)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_y_' &
                                //num2str(l)//'.bin', store=321)
                            call output_array(snapvz(1:nx, 1:ny, 1:nz), tidy(dir_working)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_z_' &
                                //num2str(l)//'.bin', store=321)
                        end if

                        l = l - 1
                    end if
                end if

            end if

            ! -------------- Adjoint wavefield reverse-time propagation ----------------------
            if (yn_free_surface) then
                call update_wavefield_free_surface(-dt, &
                    stressxxr, stressyyr, stresszzr, &
                    stressyzr, stressxzr, stressxyr, &
                    vxr, vyr, vzr, &
                    memory_pdxxxr, memory_pdyxyr, memory_pdzxzr, &
                    memory_pdxxyr, memory_pdyyyr, memory_pdzyzr, &
                    memory_pdxxzr, memory_pdyyzr, memory_pdzzzr, &
                    memory_pdxvxr, memory_pdyvxr, memory_pdzvxr, &
                    memory_pdxvyr, memory_pdyvyr, memory_pdzvyr, &
                    memory_pdxvzr, memory_pdyvzr, memory_pdzvzr)
            else
                call update_wavefield(-dt, &
                    stressxxr, stressyyr, stresszzr, &
                    stressyzr, stressxzr, stressxyr, &
                    vxr, vyr, vzr, &
                    memory_pdxxxr, memory_pdyxyr, memory_pdzxzr, &
                    memory_pdxxyr, memory_pdyyyr, memory_pdzyzr, &
                    memory_pdxxzr, memory_pdyyzr, memory_pdzzzr, &
                    memory_pdxvxr, memory_pdyvxr, memory_pdzvxr, &
                    memory_pdxvyr, memory_pdyvyr, memory_pdzvyr, &
                    memory_pdxvzr, memory_pdyvzr, memory_pdzvzr)
            end if

            ! Read and add adjoint source
            !$omp parallel private(ir, irx, iry, irz, rgx, rgy, rgz)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    if (yn_compx) then
                        rgx = sgmtr%recr(ir)%hx
                        rgy = sgmtr%recr(ir)%gy
                        rgz = sgmtr%recr(ir)%gz
                        !$omp do collapse(3)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(rgx + irx, rgy + iry, rgz + irz) &
                                            .and. ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then
                                        vxr(rgx + irx, rgy + iry, rgz + irz) = &
                                            vxr(rgx + irx, rgy + iry, rgz + irz) &
                                            + this%seis_vxr%trace(ir)%data(t) &
                                            *sgmtr%recr(ir)%interp_hx(irx) &
                                            *sgmtr%recr(ir)%interp_iy(iry) &
                                            *sgmtr%recr(ir)%interp_iz(irz) &
                                            *sgmtr%recr(ir)%weight
                                    end if
                                end do
                            end do
                        end do
                        !$omp end do
                    end if

                    if (yn_compy) then
                        rgx = sgmtr%recr(ir)%gx
                        rgy = sgmtr%recr(ir)%hy
                        rgz = sgmtr%recr(ir)%gz
                        !$omp do collapse(3)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(rgx + irx, rgy + iry, rgz + irz) &
                                            .and. ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then
                                        vyr(rgx + irx, rgy + iry, rgz + irz) = &
                                            vyr(rgx + irx, rgy + iry, rgz + irz) &
                                            + this%seis_vyr%trace(ir)%data(t) &
                                            *sgmtr%recr(ir)%interp_ix(irx) &
                                            *sgmtr%recr(ir)%interp_hy(iry) &
                                            *sgmtr%recr(ir)%interp_iz(irz) &
                                            *sgmtr%recr(ir)%weight
                                    end if
                                end do
                            end do
                        end do
                        !$omp end do
                    end if

                    if (yn_compz) then
                        rgx = sgmtr%recr(ir)%gx
                        rgy = sgmtr%recr(ir)%gy
                        rgz = sgmtr%recr(ir)%hz
                        !$omp do collapse(3)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(rgx + irx, rgy + iry, rgz + irz) &
                                            .and. ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then
                                        vzr(rgx + irx, rgy + iry, rgz + irz) = &
                                            vzr(rgx + irx, rgy + iry, rgz + irz) &
                                            + this%seis_vzr%trace(ir)%data(t) &
                                            *sgmtr%recr(ir)%interp_ix(irx) &
                                            *sgmtr%recr(ir)%interp_iy(iry) &
                                            *sgmtr%recr(ir)%interp_hz(irz) &
                                            *sgmtr%recr(ir)%weight
                                    end if

                                end do
                            end do
                        end do
                        !$omp end do
                    end if

                end if
            end do
            !$omp end parallel

            ! Compute gradients
            if (mod(t, cc_step_interval) == 0) then
                if (yn_update_medium .and. t >= sgmtr%srcr(1)%hnt) then
                    call compute_gradient
                end if
                if (yn_update_source) then
                    call compute_gradient_source(t)
                end if
            end if

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then

                wnan = group_and(any(isnan(vxr)) .or. any(isnan(vyr)) .or. any(isnan(vzr)))

                wmin1 = group_min(minval(vxr))
                wmax1 = group_max(maxval(vxr))

                wmin2 = group_min(minval(vyr))
                wmax2 = group_max(maxval(vyr))

                wmin3 = group_min(minval(vzr))
                wmax3 = group_max(maxval(vzr))

                if (rankid_group == 0) then
                    call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                        //' FWI gradient computation step '//num2str(t)//' of '//num2str(nt))
                    if (wnan) then
                        call warn(date_time_compact()//' >> Vxr, Vyr, Vzr contain NaN!')
                        stop
                    else
                        call warn(date_time_compact()//' >> Vxr, Vyr, Vzr value range = ')
                        call warn(date_time_compact()//'      '//num2str(wmin1, '(es)')//' ~ '//num2str(wmax1, '(es)'))
                        call warn(date_time_compact()//'      '//num2str(wmin2, '(es)')//' ~ '//num2str(wmax2, '(es)'))
                        call warn(date_time_compact()//'      '//num2str(wmin3, '(es)')//' ~ '//num2str(wmax3, '(es)'))
                    end if
                end if

            end if

        end do

        ! Delete temporary files
        call close_boundary_saving(delete=.true.)

        ! Output source parameter gradient
        if (yn_update_source) then
            call allreduce_array_group(grad_mt)
            if (rankid_group == 0) then
                call grd%init(n=[nc_mt, 1, 1], d=[1.0, 1.0, 1.0], o=[0.0, 0.0, 0.0])
                grd%array = -reshape(grad_mt/maxval(abs(grad_mt)), [nc_mt, 1, 1])
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_mt.grd')
            end if
        end if
        call mpibarrier_group

        if (.not. yn_update_medium) then
            return
        end if

        ! Output medium parameter gradient
        if (yn_energy_precond) then

            call allreduce_array_group(grad_c11)
            call allreduce_array_group(grad_c12)
            call allreduce_array_group(grad_c13)
            call allreduce_array_group(grad_c22)
            call allreduce_array_group(grad_c23)
            call allreduce_array_group(grad_c33)
            call allreduce_array_group(grad_c44)
            call allreduce_array_group(grad_c55)
            call allreduce_array_group(grad_c66)
            call allreduce_array_group(grad_rho)

            call allreduce_array_group(energy_src_v)
            call allreduce_array_group(energy_rec_v)
            call allreduce_array_group(energy_src_a)
            call allreduce_array_group(energy_rec_a)

            energy_src_v = energy_src_v + 1.0e-3*maxval(energy_src_v)
            energy_rec_v = energy_rec_v + 1.0e-3*maxval(energy_rec_v)
            energy_src_v = sqrt(energy_src_v*energy_rec_v)
            grad_c11 = grad_c11/energy_src_v
            grad_c12 = grad_c12/energy_src_v
            grad_c13 = grad_c13/energy_src_v
            grad_c22 = grad_c22/energy_src_v
            grad_c23 = grad_c23/energy_src_v
            grad_c33 = grad_c33/energy_src_v
            grad_c44 = grad_c44/energy_src_v
            grad_c55 = grad_c55/energy_src_v
            grad_c66 = grad_c66/energy_src_v

            energy_src_a = energy_src_a + 1.0e-3*maxval(energy_src_a)
            energy_rec_a = energy_rec_a + 1.0e-3*maxval(energy_rec_a)
            energy_src_a = sqrt(energy_src_a*energy_rec_a)
            grad_rho = grad_rho/energy_src_a

        end if

        if (rankid_group == 0) then

            nx = this%nx
            ny = this%ny
            nz = this%nz

            ! For free-surface model, map computed gradients to regular mesh
            if (yn_free_surface) then
                call map_irregular_to_regular(grad_c11, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_c12, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_c13, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_c22, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_c23, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_c33, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_c44, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_c55, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_c66, this, [1, nx, 1, ny, 1, nz])
                call map_irregular_to_regular(grad_rho, this, [1, nx, 1, ny, 1, nz])
            end if

            call grd%init(n=[nz, ny, nx], d=[dz, dy, dx], o=[oz, oy, ox])

            if (yn_free_surface) then
                rho = rho*1.0e-3
                call map_irregular_to_regular(rho, this, [1, nx, 1, ny, 1, nz])
            else
                rho = permute(this%rho, 321)*1.0e-3
            end if

            select case (aniso_param)

                case ('iso')
                    ! Using chain rule to compute gradients of Vp and Vs from
                    ! gradients of Cij. See Vigh et al. (2014)

                    vp = vp(1:nx, 1:ny, 1:nz)*1.0e-3
                    vs = vs(1:nx, 1:ny, 1:nz)*1.0e-3

                    grad_vp = grad_c11 + grad_c12 + grad_c13 + grad_c22 + grad_c23 + grad_c33
                    grad_vp = 2*rho*vp*grad_vp

                    grad_vs = grad_c44 + grad_c55 + grad_c66 - 2*grad_c12 - 2*grad_c13 - 2*grad_c23
                    grad_vs = 2*rho*vs*grad_vs

                    grd%array = permute(return_normal(grad_vp), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                    grd%array = permute(return_normal(grad_vs), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

                case ('thomsen')

                    vp = vp(1:nx, 1:ny, 1:nz)*1.0e-3
                    vs = vs(1:nx, 1:ny, 1:nz)*1.0e-3
                    tieps = tieps(1:nx, 1:ny, 1:nz)
                    tidel = tidel(1:nx, 1:ny, 1:nz)
                    tigam = tigam(1:nx, 1:ny, 1:nz)
                    tithe = -tithe(1:nx, 1:ny, 1:nz)
                    tiphi = -tiphi(1:nx, 1:ny, 1:nz)

                    grad = &
                        +thomsen_dc11_dvp*grad_c11 &
                        + thomsen_dc12_dvp*grad_c12 &
                        + thomsen_dc13_dvp*grad_c13 &
                        + thomsen_dc22_dvp*grad_c22 &
                        + thomsen_dc23_dvp*grad_c23 &
                        + thomsen_dc33_dvp*grad_c33 &
                        + thomsen_dc44_dvp*grad_c44 &
                        + thomsen_dc55_dvp*grad_c55 &
                        + thomsen_dc66_dvp*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                    grad = &
                        +thomsen_dc11_dvs*grad_c11 &
                        + thomsen_dc12_dvs*grad_c12 &
                        + thomsen_dc13_dvs*grad_c13 &
                        + thomsen_dc22_dvs*grad_c22 &
                        + thomsen_dc23_dvs*grad_c23 &
                        + thomsen_dc33_dvs*grad_c33 &
                        + thomsen_dc44_dvs*grad_c44 &
                        + thomsen_dc55_dvs*grad_c55 &
                        + thomsen_dc66_dvs*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

                    grad = &
                        +thomsen_dc11_deps*grad_c11 &
                        + thomsen_dc12_deps*grad_c12 &
                        + thomsen_dc13_deps*grad_c13 &
                        + thomsen_dc22_deps*grad_c22 &
                        + thomsen_dc23_deps*grad_c23 &
                        + thomsen_dc33_deps*grad_c33 &
                        + thomsen_dc44_deps*grad_c44 &
                        + thomsen_dc55_deps*grad_c55 &
                        + thomsen_dc66_deps*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_epsilon.grd')

                    grad = &
                        +thomsen_dc11_ddel*grad_c11 &
                        + thomsen_dc12_ddel*grad_c12 &
                        + thomsen_dc13_ddel*grad_c13 &
                        + thomsen_dc22_ddel*grad_c22 &
                        + thomsen_dc23_ddel*grad_c23 &
                        + thomsen_dc33_ddel*grad_c33 &
                        + thomsen_dc44_ddel*grad_c44 &
                        + thomsen_dc55_ddel*grad_c55 &
                        + thomsen_dc66_ddel*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_delta.grd')

                    grad = &
                        +thomsen_dc11_dgam*grad_c11 &
                        + thomsen_dc12_dgam*grad_c12 &
                        + thomsen_dc13_dgam*grad_c13 &
                        + thomsen_dc22_dgam*grad_c22 &
                        + thomsen_dc23_dgam*grad_c23 &
                        + thomsen_dc33_dgam*grad_c33 &
                        + thomsen_dc44_dgam*grad_c44 &
                        + thomsen_dc55_dgam*grad_c55 &
                        + thomsen_dc66_dgam*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_gamma.grd')

                case ('a-t')

                    vp = vp(1:nx, 1:ny, 1:nz)*1.0e-3
                    vs = vs(1:nx, 1:ny, 1:nz)*1.0e-3
                    tieps = tieps(1:nx, 1:ny, 1:nz)
                    tieta = tieta(1:nx, 1:ny, 1:nz)
                    tigam = tigam(1:nx, 1:ny, 1:nz)
                    tithe = -tithe(1:nx, 1:ny, 1:nz)
                    tiphi = -tiphi(1:nx, 1:ny, 1:nz)

                    grad = &
                        +alkhalifah_tsvankin_dc11_dvp*grad_c11 &
                        + alkhalifah_tsvankin_dc12_dvp*grad_c12 &
                        + alkhalifah_tsvankin_dc13_dvp*grad_c13 &
                        + alkhalifah_tsvankin_dc22_dvp*grad_c22 &
                        + alkhalifah_tsvankin_dc23_dvp*grad_c23 &
                        + alkhalifah_tsvankin_dc33_dvp*grad_c33 &
                        + alkhalifah_tsvankin_dc44_dvp*grad_c44 &
                        + alkhalifah_tsvankin_dc55_dvp*grad_c55 &
                        + alkhalifah_tsvankin_dc66_dvp*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                    grad = &
                        +alkhalifah_tsvankin_dc11_dvs*grad_c11 &
                        + alkhalifah_tsvankin_dc12_dvs*grad_c12 &
                        + alkhalifah_tsvankin_dc13_dvs*grad_c13 &
                        + alkhalifah_tsvankin_dc22_dvs*grad_c22 &
                        + alkhalifah_tsvankin_dc23_dvs*grad_c23 &
                        + alkhalifah_tsvankin_dc33_dvs*grad_c33 &
                        + alkhalifah_tsvankin_dc44_dvs*grad_c44 &
                        + alkhalifah_tsvankin_dc55_dvs*grad_c55 &
                        + alkhalifah_tsvankin_dc66_dvs*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

                    grad = &
                        +alkhalifah_tsvankin_dc11_deps*grad_c11 &
                        + alkhalifah_tsvankin_dc12_deps*grad_c12 &
                        + alkhalifah_tsvankin_dc13_deps*grad_c13 &
                        + alkhalifah_tsvankin_dc22_deps*grad_c22 &
                        + alkhalifah_tsvankin_dc23_deps*grad_c23 &
                        + alkhalifah_tsvankin_dc33_deps*grad_c33 &
                        + alkhalifah_tsvankin_dc44_deps*grad_c44 &
                        + alkhalifah_tsvankin_dc55_deps*grad_c55 &
                        + alkhalifah_tsvankin_dc66_deps*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_epsilon.grd')

                    grad = &
                        +alkhalifah_tsvankin_dc11_deta*grad_c11 &
                        + alkhalifah_tsvankin_dc12_deta*grad_c12 &
                        + alkhalifah_tsvankin_dc13_deta*grad_c13 &
                        + alkhalifah_tsvankin_dc22_deta*grad_c22 &
                        + alkhalifah_tsvankin_dc23_deta*grad_c23 &
                        + alkhalifah_tsvankin_dc33_deta*grad_c33 &
                        + alkhalifah_tsvankin_dc44_deta*grad_c44 &
                        + alkhalifah_tsvankin_dc55_deta*grad_c55 &
                        + alkhalifah_tsvankin_dc66_deta*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_eta.grd')

                    grad = &
                        +alkhalifah_tsvankin_dc11_dgam*grad_c11 &
                        + alkhalifah_tsvankin_dc12_dgam*grad_c12 &
                        + alkhalifah_tsvankin_dc13_dgam*grad_c13 &
                        + alkhalifah_tsvankin_dc22_dgam*grad_c22 &
                        + alkhalifah_tsvankin_dc23_dgam*grad_c23 &
                        + alkhalifah_tsvankin_dc33_dgam*grad_c33 &
                        + alkhalifah_tsvankin_dc44_dgam*grad_c44 &
                        + alkhalifah_tsvankin_dc55_dgam*grad_c55 &
                        + alkhalifah_tsvankin_dc66_dgam*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_gamma.grd')

                case ('cij')
                    grd%array = permute(return_normal(grad_c11), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c11.grd')
                    grd%array = permute(return_normal(grad_c12), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c12.grd')
                    grd%array = permute(return_normal(grad_c13), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c13.grd')
                    grd%array = permute(return_normal(grad_c22), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c22.grd')
                    grd%array = permute(return_normal(grad_c23), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c23.grd')
                    grd%array = permute(return_normal(grad_c33), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c33.grd')
                    grd%array = permute(return_normal(grad_c44), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c44.grd')
                    grd%array = permute(return_normal(grad_c55), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c55.grd')
                    grd%array = permute(return_normal(grad_c66), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c66.grd')

            end select

            grd%array = permute(return_normal(grad_rho), 321)
            call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_rho.grd')

        end if

        call mpibarrier_group

    end subroutine

    !
    !> Compute medium parameter gradients
    !
    subroutine compute_gradient

        integer :: i, j, k
        real :: tmpxx, tmpxy, tmpxz, tmpyy, tmpyz, tmpzz
        real :: tmpxxr, tmpxyr, tmpxzr, tmpyyr, tmpyzr, tmpzzr

        !$omp parallel do private(i, j, k, tmpxx, tmpxy, tmpxz, tmpyy, tmpyz, tmpzz, &
            !$omp   tmpxxr, tmpxyr, tmpxzr, tmpyyr, tmpyzr, tmpzzr) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    tmpxx = stressxx(i, j, k) - prev_stressxx(i, j, k)
                    tmpyy = stressyy(i, j, k) - prev_stressyy(i, j, k)
                    tmpzz = stresszz(i, j, k) - prev_stresszz(i, j, k)
                    tmpyz = 0.25*sum(stressyz(i, j:j + 1, k:k + 1)) - 0.25*sum(prev_stressyz(i, j:j + 1, k:k + 1))
                    tmpxz = 0.25*sum(stressxz(i:i + 1, j, k:k + 1)) - 0.25*sum(prev_stressxz(i:i + 1, j, k:k + 1))
                    tmpxy = 0.25*sum(stressxy(i:i + 1, j:j + 1, k)) - 0.25*sum(prev_stressxy(i:i + 1, j:j + 1, k))

                    tmpxxr = stressxxr(i, j, k)
                    tmpyyr = stressyyr(i, j, k)
                    tmpzzr = stresszzr(i, j, k)
                    tmpyzr = 0.25*sum(stressyzr(i, j:j + 1, k:k + 1))
                    tmpxzr = 0.25*sum(stressxzr(i:i + 1, j, k:k + 1))
                    tmpxyr = 0.25*sum(stressxyr(i:i + 1, j:j + 1, k))

                    strainxx(i, j, k) = s11(i, j, k)*tmpxx + s12(i, j, k)*tmpyy + s13(i, j, k)*tmpzz
                    strainyy(i, j, k) = s12(i, j, k)*tmpxx + s22(i, j, k)*tmpyy + s23(i, j, k)*tmpzz
                    strainzz(i, j, k) = s13(i, j, k)*tmpxx + s23(i, j, k)*tmpyy + s33(i, j, k)*tmpzz
                    strainyz(i, j, k) = s44(i, j, k)*tmpyz
                    strainxz(i, j, k) = s55(i, j, k)*tmpxz
                    strainxy(i, j, k) = s66(i, j, k)*tmpxy

                    strainxxr(i, j, k) = s11(i, j, k)*tmpxxr + s12(i, j, k)*tmpyyr + s13(i, j, k)*tmpzzr
                    strainyyr(i, j, k) = s12(i, j, k)*tmpxxr + s22(i, j, k)*tmpyyr + s23(i, j, k)*tmpzzr
                    strainzzr(i, j, k) = s13(i, j, k)*tmpxxr + s23(i, j, k)*tmpyyr + s33(i, j, k)*tmpzzr
                    strainyzr(i, j, k) = s44(i, j, k)*tmpyzr
                    strainxzr(i, j, k) = s55(i, j, k)*tmpxzr
                    strainxyr(i, j, k) = s66(i, j, k)*tmpxyr

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
        do k = nz1_interior, nz2_interior
            do j = ny1_interior, ny2_interior
                do i = nx1_interior, nx2_interior

                    grad_c11(i, j, k) = grad_c11(i, j, k) - strainxx(i, j, k)*strainxxr(i, j, k)
                    grad_c12(i, j, k) = grad_c12(i, j, k) - strainxx(i, j, k)*strainyyr(i, j, k) - strainyy(i, j, k)*strainxxr(i, j, k)
                    grad_c13(i, j, k) = grad_c13(i, j, k) - strainxx(i, j, k)*strainzzr(i, j, k) - strainzz(i, j, k)*strainxxr(i, j, k)
                    grad_c22(i, j, k) = grad_c22(i, j, k) - strainyy(i, j, k)*strainyyr(i, j, k)
                    grad_c23(i, j, k) = grad_c23(i, j, k) - strainyy(i, j, k)*strainzzr(i, j, k) - strainzz(i, j, k)*strainyyr(i, j, k)
                    grad_c33(i, j, k) = grad_c33(i, j, k) - strainzz(i, j, k)*strainzzr(i, j, k)
                    grad_c44(i, j, k) = grad_c44(i, j, k) - strainyz(i, j, k)*strainyzr(i, j, k)
                    grad_c55(i, j, k) = grad_c55(i, j, k) - strainxz(i, j, k)*strainxzr(i, j, k)
                    grad_c66(i, j, k) = grad_c66(i, j, k) - strainxy(i, j, k)*strainxyr(i, j, k)

                end do
            end do
        end do
        !$omp end parallel do

        if (yn_energy_precond) then

            !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
            do k = nz1_interior, nz2_interior
                do j = ny1_interior, ny2_interior
                    do i = nx1_interior, nx2_interior
                        energy_src_v(i, j, k) = energy_src_v(i, j, k) &
                            + strainxx(i, j, k)**2 + strainyy(i, j, k)**2 + strainzz(i, j, k)**2 &
                            + 2*strainyz(i, j, k)**2 + 2*strainxz(i, j, k)**2 + 2*strainxy(i, j, k)**2
                        energy_rec_v(i, j, k) = energy_rec_v(i, j, k) &
                            + strainxxr(i, j, k)**2 + strainyyr(i, j, k)**2 + strainzzr(i, j, k)**2 &
                            + 2*strainyzr(i, j, k)**2 + 2*strainxzr(i, j, k)**2 + 2*strainxyr(i, j, k)**2
                    end do
                end do
            end do
            !$omp end parallel do

        end if

        !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2
                    src_vx(i, j, k) = sum(vx(i:i + 1, j, k)) - sum(prev_vx(i:i + 1, j, k))
                    rec_vx(i, j, k) = sum(vxr(i:i + 1, j, k))
                    src_vy(i, j, k) = sum(vy(i, j:j + 1, k)) - sum(prev_vy(i, j:j + 1, k))
                    rec_vy(i, j, k) = sum(vyr(i, j:j + 1, k))
                    src_vz(i, j, k) = sum(vz(i, j, k:k + 1)) - sum(prev_vz(i, j, k:k + 1))
                    rec_vz(i, j, k) = sum(vzr(i, j, k:k + 1))
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
        do k = nz1_interior, nz2_interior
            do j = ny1_interior, ny2_interior
                do i = nx1_interior, nx2_interior
                    grad_rho(i, j, k) = grad_rho(i, j, k) + src_vx(i, j, k)*rec_vx(i, j, k) &
                        + src_vy(i, j, k)*rec_vy(i, j, k) + src_vz(i, j, k)*rec_vz(i, j, k)
                end do
            end do
        end do
        !$omp end parallel do

        if (yn_energy_precond) then

            !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
            do k = nz1_interior, nz2_interior
                do j = ny1_interior, ny2_interior
                    do i = nx1_interior, nx2_interior
                        energy_src_a(i, j, k) = energy_src_a(i, j, k) + src_vx(i, j, k)**2 + src_vy(i, j, k)**2 + src_vz(i, j, k)**2
                        energy_rec_a(i, j, k) = energy_rec_a(i, j, k) + rec_vx(i, j, k)**2 + rec_vy(i, j, k)**2 + rec_vz(i, j, k)**2
                    end do
                end do
            end do
            !$omp end parallel do

        end if

    end subroutine

    !
    !> Compute source parameter gradients
    !
    subroutine compute_gradient_source(t)

        integer, intent(in) :: t

        integer :: sgx, sgy, sgz, irx, iry, irz
        integer :: i

        do i = 1, sgmtr%ns

            sgx = sgmtr%srcr(i)%gx
            sgy = sgmtr%srcr(i)%gy
            sgz = sgmtr%srcr(i)%gz
            do irz = -nkw, nkw
                do iry = -nkw, nkw
                    do irx = -nkw, nkw
                        if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                            grad_mt(1) = grad_mt(1) - &
                                stressxxr(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(2) = grad_mt(2) - &
                                stressyyr(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(3) = grad_mt(3) - &
                                stresszzr(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                        end if
                    end do
                end do
            end do

            sgx = sgmtr%srcr(i)%hx
            sgy = sgmtr%srcr(i)%hy
            sgz = sgmtr%srcr(i)%gz
            do irz = -nkw, nkw
                do iry = -nkw, nkw
                    do irx = -nkw, nkw
                        if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                            grad_mt(4) = grad_mt(4) - &
                                stressxyr(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                        end if
                    end do
                end do
            end do

            sgx = sgmtr%srcr(i)%hx
            sgy = sgmtr%srcr(i)%gy
            sgz = sgmtr%srcr(i)%hz
            do irz = -nkw, nkw
                do iry = -nkw, nkw
                    do irx = -nkw, nkw
                        if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                            grad_mt(5) = grad_mt(5) - &
                                stressxzr(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                        end if
                    end do
                end do
            end do

            sgx = sgmtr%srcr(i)%gx
            sgy = sgmtr%srcr(i)%hy
            sgz = sgmtr%srcr(i)%hz
            do irz = -nkw, nkw
                do iry = -nkw, nkw
                    do irx = -nkw, nkw
                        if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                            grad_mt(6) = grad_mt(6) - &
                                stressyzr(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                        end if
                    end do
                end do
            end do

        end do

    end subroutine

end submodule
