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


submodule(elastic_vhtiort_2d) elastic_vhtiort_2d_fwi_gradient

    use libflit
    use elastic_vhtiort_2d_vars
    use elastic_vhtiort_2d_boundary_saving
    use elastic_vhtiort_2d_wavefield

    implicit none

#include '../../lib/macro_thomsen_2d.f90'
#include '../../lib/macro_alkhalifah_tsvankin_2d.f90'

contains

    module subroutine compute_adjoint(this)

        class(wave_solver_elastic_vhtiort_2d), intent(inout) :: this

        integer :: l, ir, irx, irz, rgx, rgz
        type(grid2) :: grd
        integer :: i, j, t

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield
        call alloc_adjoint_wavefield

        yn_energy_precond = this%energy_precond
        kernel_v = this%kernel_v
        kernel_a = this%kernel_a

        ! Adjoint source
        if (yn_compx) then
            call this%seis_vxr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_x.su')
            call this%seis_vxr%resamp(nnt=nt, ddt=dt)
        end if

        if (yn_compz) then
            call this%seis_vzr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_z.su')
            call this%seis_vzr%resamp(nnt=nt, ddt=dt)
        end if

        ! Elastic compliances
        call alloc_array(s11, [1, nx, 1, nz], pad=pml + 1)
        call alloc_array(s13, [1, nx, 1, nz], pad=pml + 1)
        call alloc_array(s33, [1, nx, 1, nz], pad=pml + 1)
        call alloc_array(s55, [1, nx, 1, nz], pad=pml + 1)
        where (c13**2 - c11*c33 /= 0)
            s11 = c33/(-c13**2 + c11*c33)
            s13 = c13/(c13**2 - c11*c33)
            s33 = c11/(-c13**2 + c11*c33)
        end where
        where (c55 /= 0)
            s55 = 1.0/c55
        end where
        where (c55 == 0)
            s11 = 0
            s13 = 0
            s33 = 0
            s55 = 0
        end where

        call alloc_array(snapvx, [1, nx, 1, nz], pad=pml)
        call alloc_array(snapvz, [1, nx, 1, nz], pad=pml)

        grad_c11 = zeros(nx, nz)
        grad_c13 = zeros(nx, nz)
        grad_c33 = zeros(nx, nz)
        grad_c55 = zeros(nx, nz)
        grad_rho = zeros(nx, nz)
        grad_mt = zeros(nc_mt)

        energy_src_v = zeros(nx, nz)
        energy_rec_v = zeros(nx, nz)
        energy_src_a = zeros(nx, nz)
        energy_rec_a = zeros(nx, nz)

        ! Source wavefields
        call alloc_array(strainxx, [1, nx, 1, nz], pad=pml)
        call alloc_array(strainzz, [1, nx, 1, nz], pad=pml)
        call alloc_array(strainxz, [1, nx, 1, nz], pad=pml)
        call alloc_array(strainxxr, [1, nx, 1, nz], pad=pml)
        call alloc_array(strainzzr, [1, nx, 1, nz], pad=pml)
        call alloc_array(strainxzr, [1, nx, 1, nz], pad=pml)
        call alloc_array(src_vx, [1, nx, 1, nz], pad=pml)
        call alloc_array(rec_vx, [1, nx, 1, nz], pad=pml)
        call alloc_array(src_vz, [1, nx, 1, nz], pad=pml)
        call alloc_array(rec_vz, [1, nx, 1, nz], pad=pml)

        ! Previous step wavefields
        call alloc_array(prev_stressxx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_stresszz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_stressxz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_vx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_vz, [1, nx, 1, nz], pad=pml + fdhalf)

        call alloc_array(strainxx_lrsh, [1, nx], pad=pml)
        call alloc_array(strainxx_lrrh, [1, nx], pad=pml)
        call alloc_array(strainzz_lrsh, [1, nx], pad=pml)
        call alloc_array(strainzz_lrrh, [1, nx], pad=pml)
        call alloc_array(strainxz_lrsh, [1, nx], pad=pml)
        call alloc_array(strainxz_lrrh, [1, nx], pad=pml)

        call alloc_array(strainxx_udsh, [1, nz], pad=pml)
        call alloc_array(strainxx_udrh, [1, nz], pad=pml)
        call alloc_array(strainzz_udsh, [1, nz], pad=pml)
        call alloc_array(strainzz_udrh, [1, nz], pad=pml)
        call alloc_array(strainxz_udsh, [1, nz], pad=pml)
        call alloc_array(strainxz_udrh, [1, nz], pad=pml)

        call alloc_array(p_lrsh, [1, nx], pad=pml)
        call alloc_array(p_lrrh, [1, nx], pad=pml)
        call alloc_array(p_udsh, [1, nz], pad=pml)
        call alloc_array(p_udrh, [1, nz], pad=pml)

        ! Prepare boundary saving
        call prepare_boundary_saving
        call open_boundary_saving

        l = np
        do t = nt, 1, -1

            if (yn_update_medium .and. t >= sgmtr%srcr(1)%hnt) then

                prev_stressxx = stressxx
                prev_stresszz = stresszz
                prev_stressxz = stressxz
                prev_vx = vx
                prev_vz = vz

                ! -------------- Forward wavefield reconstruction -----------------------
                if (yn_free_surface) then
                    call update_wavefield_free_surface(-dt, &
                        stressxx, stresszz, stressxz, vx, vz, &
                        memory_pdxvx, memory_pdzvx, memory_pdxvz, memory_pdzvz, &
                        memory_pdxxx, memory_pdzxz, memory_pdxxz, memory_pdzzz)
                else
                    call update_wavefield(-dt, &
                        stressxx, stresszz, stressxz, vx, vz, &
                        memory_pdxvx, memory_pdzvx, memory_pdxvz, memory_pdzvz, &
                        memory_pdxxx, memory_pdzxz, memory_pdxxz, memory_pdzzz)
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

                        !$omp parallel do private(i, j) collapse(2)
                        do j = 1, nz
                            do i = 1, nx
                                snapvx(i, j) = 0.5*sum(vx(i:i + 1, j))
                                snapvz(i, j) = 0.5*sum(vz(i, j:j + 1))
                            end do
                        end do
                        !$omp end parallel do

                        call output_array(snapvx(1:nx, 1:nz), tidy(dir_working)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_reconstructed_wavefield_x_' &
                            //num2str(l)//'.bin', transp=.true.)
                        call output_array(snapvz(1:nx, 1:nz), tidy(dir_working)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_reconstructed_wavefield_z_' &
                            //num2str(l)//'.bin', transp=.true.)

                        l = l - 1
                    end if
                end if

            end if

            ! -------------- Adjoint wavefield backward propagation -----------------------
            if (yn_free_surface) then
                call update_wavefield_free_surface(-dt, &
                    stressxxr, stresszzr, stressxzr, vxr, vzr, &
                    memory_pdxvxr, memory_pdzvxr, memory_pdxvzr, memory_pdzvzr, &
                    memory_pdxxxr, memory_pdzxzr, memory_pdxxzr, memory_pdzzzr)
            else
                call update_wavefield(-dt, &
                    stressxxr, stresszzr, stressxzr, vxr, vzr, &
                    memory_pdxvxr, memory_pdzvxr, memory_pdxvzr, memory_pdzvzr, &
                    memory_pdxxxr, memory_pdzxzr, memory_pdxxzr, memory_pdzzzr)
            end if

            ! Read and add adjoint source
            !$omp parallel private(ir, irx, irz, rgx, rgz)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    if (yn_compx) then

                        rgx = sgmtr%recr(ir)%hx
                        rgz = sgmtr%recr(ir)%gz

                        !$omp do collapse(2)
                        do irz = -nkw, nkw
                            do irx = -nkw, nkw
                                if (ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then
                                    vxr(rgx + irx, rgz + irz) = vxr(rgx + irx, rgz + irz) &
                                        + this%seis_vxr%trace(ir)%data(t) &
                                        *sgmtr%recr(ir)%interp_hx(irx) &
                                        *sgmtr%recr(ir)%interp_iz(irz) &
                                        *sgmtr%recr(ir)%weight
                                end if
                            end do
                        end do
                        !$omp end do

                    end if

                    if (yn_compz) then

                        rgx = sgmtr%recr(ir)%gx
                        rgz = sgmtr%recr(ir)%hz

                        !$omp do collapse(2)
                        do irz = -nkw, nkw
                            do irx = -nkw, nkw
                                if (ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then
                                    vzr(rgx + irx, rgz + irz) = vzr(rgx + irx, rgz + irz) &
                                        + this%seis_vzr%trace(ir)%data(t) &
                                        *sgmtr%recr(ir)%interp_ix(irx) &
                                        *sgmtr%recr(ir)%interp_hz(irz) &
                                        *sgmtr%recr(ir)%weight
                                end if
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
                call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                    //' FWI gradient computation step '//num2str(t)//' of '//num2str(nt))
                if (any(isnan(vxr)) .or. any(isnan(vzr))) then
                    call warn(date_time_compact()//' >> Vxr, Vzr contain NaN!')
                    stop
                else
                    call warn(date_time_compact()//' >> Vxr, Vzr value range = ')
                    call warn(date_time_compact()//'      '//num2str(minval(vxr), '(es)') &
                        //' ~ '//num2str(maxval(vxr), '(es)'))
                    call warn(date_time_compact()//'      '//num2str(minval(vzr), '(es)') &
                        //' ~ '//num2str(maxval(vzr), '(es)'))
                end if
            end if

        end do

        ! Delete temporary files
        call close_boundary_saving(delete=.true.)

        ! Output source parameter gradient
        if (yn_update_source) then
            call grd%init(n=[nc_mt, 1], d=[1.0, 1.0], o=[0.0, 0.0])
            grd%array = -reshape(grad_mt, [nc_mt, 1])
            call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_mt.grd')
        end if

        if (.not. yn_update_medium) then
            return
        end if

        ! Output medium parameter gradient
        if (yn_energy_precond) then

            if (kernel_v /= '') then
                energy_src_v = energy_src_v + 1.0e-3*maxval(energy_src_v)
                energy_rec_v = energy_rec_v + 1.0e-3*maxval(energy_rec_v)
                energy_src_v = sqrt(energy_src_v*energy_rec_v)
                grad_c11 = grad_c11/energy_src_v
                grad_c13 = grad_c13/energy_src_v
                grad_c33 = grad_c33/energy_src_v
                grad_c55 = grad_c55/energy_src_v
            end if

            if (kernel_a /= '') then
                energy_src_a = energy_src_a + 1.0e-3*maxval(energy_src_a)
                energy_rec_a = energy_rec_a + 1.0e-3*maxval(energy_rec_a)
                energy_src_a = sqrt(energy_src_a*energy_rec_a)
                grad_rho = grad_rho/energy_src_a
            end if

        end if

        nx = this%nx
        nz = this%nz

        ! For free-surface model, map computed gradients to regular mesh
        if (yn_free_surface) then
            if (kernel_v /= '') then
                call map_irregular_to_regular(grad_c11, this, [1, nx, 1, nz])
                call map_irregular_to_regular(grad_c13, this, [1, nx, 1, nz])
                call map_irregular_to_regular(grad_c33, this, [1, nx, 1, nz])
                call map_irregular_to_regular(grad_c55, this, [1, nx, 1, nz])
            end if
            if (kernel_a /= '') then
                call map_irregular_to_regular(grad_rho, this, [1, nx, 1, nz])
            end if
        end if

        ! Output
        call grd%init(n=[nz, nx], d=[dz, dx], o=[oz, ox])

        if (yn_free_surface) then
            rho = rho*1.0e-3
            call map_irregular_to_regular(rho, this, [1, nx, 1, nz])
        else
            rho = rho(1:nx, 1:nz)*1.0e-3
        end if

        select case (aniso_param)

            case ('iso')

                ! Here the models are not remaped, as at the beginning we didn't map them to irregular mesh
                ! as for Cij
                vp = vp(1:nx, 1:nz)*1.0e-3
                vs = vs(1:nx, 1:nz)*1.0e-3

                ! Using chain rule to compute gradients of Vp and Vs from
                ! gradients of Cij. See Vigh et al. (2014)
                grad_vp = grad_c11 + grad_c13 + grad_c33
                grad_vp = 2*rho*vp*grad_vp

                grad_vs = grad_c55 - 2*grad_c13
                grad_vs = 2*rho*vs*grad_vs

                grd%array = transpose(return_normal(grad_vp))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                grd%array = transpose(return_normal(grad_vs))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

            case ('thomsen')

                vp = vp(1:nx, 1:nz)*1.0e-3
                vs = vs(1:nx, 1:nz)*1.0e-3
                tieps = tieps(1:nx, 1:nz)
                tidel = tidel(1:nx, 1:nz)
                tithe = -tithe(1:nx, 1:nz)

                grad_vp = thomsen_dc11_dvp*grad_c11 &
                    + thomsen_dc13_dvp*grad_c13 &
                    + thomsen_dc33_dvp*grad_c33 &
                    + thomsen_dc55_dvp*grad_c55

                grad_vs = thomsen_dc11_dvs*grad_c11 &
                    + thomsen_dc13_dvs*grad_c13 &
                    + thomsen_dc33_dvs*grad_c33 &
                    + thomsen_dc55_dvs*grad_c55

                grad_epsilon = thomsen_dc11_deps*grad_c11 &
                    + thomsen_dc13_deps*grad_c13 &
                    + thomsen_dc33_deps*grad_c33 &
                    + thomsen_dc55_deps*grad_c55

                grad_delta = thomsen_dc11_ddel*grad_c11 &
                    + thomsen_dc13_ddel*grad_c13 &
                    + thomsen_dc33_ddel*grad_c33 &
                    + thomsen_dc55_ddel*grad_c55

                grd%array = transpose(return_normal(grad_vp))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                grd%array = transpose(return_normal(grad_vs))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

                grd%array = transpose(return_normal(grad_epsilon))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_epsilon.grd')

                grd%array = transpose(return_normal(grad_delta))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_delta.grd')

            case ('a-t')

                vp = vp(1:nx, 1:nz)*1.0e-3
                vs = vs(1:nx, 1:nz)*1.0e-3
                tieps = tieps(1:nx, 1:nz)
                tieta = tieta(1:nx, 1:nz)
                tithe = -tithe(1:nx, 1:nz)

                grad_vp = alkhalifah_tsvankin_dc11_dvp*grad_c11 &
                    + alkhalifah_tsvankin_dc13_dvp*grad_c13 &
                    + alkhalifah_tsvankin_dc33_dvp*grad_c33 &
                    + alkhalifah_tsvankin_dc55_dvp*grad_c55

                grad_vs = alkhalifah_tsvankin_dc11_dvs*grad_c11 &
                    + alkhalifah_tsvankin_dc13_dvs*grad_c13 &
                    + alkhalifah_tsvankin_dc33_dvs*grad_c33 &
                    + alkhalifah_tsvankin_dc55_dvs*grad_c55

                grad_epsilon = alkhalifah_tsvankin_dc11_deps*grad_c11 &
                    + alkhalifah_tsvankin_dc13_deps*grad_c13 &
                    + alkhalifah_tsvankin_dc33_deps*grad_c33 &
                    + alkhalifah_tsvankin_dc55_deps*grad_c55

                grad_eta = alkhalifah_tsvankin_dc11_deta*grad_c11 &
                    + alkhalifah_tsvankin_dc13_deta*grad_c13 &
                    + alkhalifah_tsvankin_dc33_deta*grad_c33 &
                    + alkhalifah_tsvankin_dc55_deta*grad_c55

                grd%array = transpose(return_normal(grad_vp))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                grd%array = transpose(return_normal(grad_vs))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

                grd%array = transpose(return_normal(grad_epsilon))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_epsilon.grd')

                grd%array = transpose(return_normal(grad_eta))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_eta.grd')

            case ('cij')

                grd%array = transpose(return_normal(grad_c11))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c11.grd')

                grd%array = transpose(return_normal(grad_c13))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c13.grd')

                grd%array = transpose(return_normal(grad_c33))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c33.grd')

                grd%array = transpose(return_normal(grad_c55))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c55.grd')

        end select

        grd%array = transpose(return_normal(grad_rho))
        call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_rho.grd')

        call mpibarrier_group

    end subroutine

    !
    !> Compute medium parameter gradients
    !
    subroutine compute_gradient

        integer :: i, j
        integer :: sgnh
        real :: tmpxx, tmpzz, tmpxz, tmpxxr, tmpzzr, tmpxzr

        !$omp parallel do private(i, j, tmpxx, tmpzz, tmpxz, tmpxxr, tmpzzr, tmpxzr) collapse(2) schedule(auto)
        do j = -pml + 1, nz + pml
            do i = -pml + 1, nx + pml

                tmpxx = stressxx(i, j) - prev_stressxx(i, j)
                tmpzz = stresszz(i, j) - prev_stresszz(i, j)
                tmpxz = 0.25*sum(stressxz(i:i + 1, j:j + 1)) - 0.25*sum(prev_stressxz(i:i + 1, j:j + 1))

                tmpxxr = stressxxr(i, j)
                tmpzzr = stresszzr(i, j)
                tmpxzr = 0.25*sum(stressxzr(i:i + 1, j:j + 1))

                strainxx(i, j) = s11(i, j)*tmpxx + s13(i, j)*tmpzz
                strainzz(i, j) = s13(i, j)*tmpxx + s33(i, j)*tmpzz
                strainxz(i, j) = s55(i, j)*tmpxz
                strainxxr(i, j) = s11(i, j)*tmpxxr + s13(i, j)*tmpzzr
                strainzzr(i, j) = s13(i, j)*tmpxxr + s33(i, j)*tmpzzr
                strainxzr(i, j) = s55(i, j)*tmpxzr

            end do
        end do
        !$omp end parallel do

        if (kernel_v /= '') then

            if (kernel_v == 'full') then

                !$omp parallel do private(i, j) collapse(2) schedule(auto)
                do j = 1, nz
                    do i = 1, nx

                        grad_c11(i, j) = grad_c11(i, j) - (strainxxr(i, j)*strainxx(i, j))
                        grad_c13(i, j) = grad_c13(i, j) - (strainxxr(i, j)*strainzz(i, j) + strainzzr(i, j)*strainxx(i, j))
                        grad_c33(i, j) = grad_c33(i, j) - (strainzzr(i, j)*strainzz(i, j))
                        grad_c55(i, j) = grad_c55(i, j) - (strainxzr(i, j)*strainxz(i, j))

                    end do
                end do
                !$omp end parallel do

            else

                ! Along x
                if (index(kernel_v, 'lowx') /= 0) then
                    sgnh = 1
                else if (index(kernel_v, 'highx') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    !$omp parallel do private(i, j, &
                        !$omp strainxx_lrsh, strainxx_lrrh, &
                        !$omp strainzz_lrsh, strainzz_lrrh, &
                        !$omp strainxz_lrsh, strainxz_lrrh) schedule(auto)
                    do j = 1, nz

                        strainxx_lrsh = strainxx(:, j)
                        strainzz_lrsh = strainzz(:, j)
                        strainxz_lrsh = strainxz(:, j)
                        strainxx_lrrh = strainxxr(:, j)
                        strainzz_lrrh = strainzzr(:, j)
                        strainxz_lrrh = strainxzr(:, j)

                        call hilbert_transform(strainxx_lrsh)
                        call hilbert_transform(strainzz_lrsh)
                        call hilbert_transform(strainxz_lrsh)
                        call hilbert_transform(strainxx_lrrh)
                        call hilbert_transform(strainzz_lrrh)
                        call hilbert_transform(strainxz_lrrh)

                        grad_c11(1:nx, j) = grad_c11(1:nx, j) &
                            - (strainxx(1:nx, j)*strainxxr(1:nx, j) + sgnh*strainxx_lrsh(1:nx)*strainxx_lrrh(1:nx))
                        grad_c13(1:nx, j) = grad_c13(1:nx, j) &
                            - (strainxx(1:nx, j)*strainzzr(1:nx, j) + sgnh*strainxx_lrsh(1:nx)*strainzz_lrrh(1:nx)) &
                            - (strainzz(1:nx, j)*strainxxr(1:nx, j) + sgnh*strainzz_lrsh(1:nx)*strainxx_lrrh(1:nx))
                        grad_c33(1:nx, j) = grad_c33(1:nx, j) &
                            - (strainzz(1:nx, j)*strainzzr(1:nx, j) + sgnh*strainzz_lrsh(1:nx)*strainzz_lrrh(1:nx))
                        grad_c55(1:nx, j) = grad_c55(1:nx, j) &
                            - (strainxz(1:nx, j)*strainxzr(1:nx, j) + sgnh*strainxz_lrsh(1:nx)*strainxz_lrrh(1:nx))

                    end do
                    !$omp end parallel do

                end if

                ! Along z
                if (index(kernel_v, 'lowz') /= 0) then
                    sgnh = 1
                else if (index(kernel_v, 'highz') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    !$omp parallel do private(i, j, &
                        !$omp strainxx_udsh, strainxx_udrh, &
                        !$omp strainzz_udsh, strainzz_udrh, &
                        !$omp strainxz_udsh, strainxz_udrh) schedule(auto)
                    do i = 1, nx

                        strainxx_udsh = strainxx(i, :)
                        strainzz_udsh = strainzz(i, :)
                        strainxz_udsh = strainxz(i, :)
                        strainxx_udrh = strainxxr(i, :)
                        strainzz_udrh = strainzzr(i, :)
                        strainxz_udrh = strainxzr(i, :)

                        call hilbert_transform(strainxx_udsh)
                        call hilbert_transform(strainzz_udsh)
                        call hilbert_transform(strainxz_udsh)
                        call hilbert_transform(strainxx_udrh)
                        call hilbert_transform(strainzz_udrh)
                        call hilbert_transform(strainxz_udrh)

                        grad_c11(i, 1:nz) = grad_c11(i, 1:nz) &
                            - (strainxx(i, 1:nz)*strainxxr(i, 1:nz) + sgnh*strainxx_udsh(1:nz)*strainxx_udrh(1:nz))
                        grad_c13(i, 1:nz) = grad_c13(i, 1:nz) &
                            - (strainxx(i, 1:nz)*strainzzr(i, 1:nz) + sgnh*strainxx_udsh(1:nz)*strainzz_udrh(1:nz)) &
                            - (strainzz(i, 1:nz)*strainxxr(i, 1:nz) + sgnh*strainzz_udsh(1:nz)*strainxx_udrh(1:nz))
                        grad_c33(i, 1:nz) = grad_c33(i, 1:nz) &
                            - (strainzz(i, 1:nz)*strainzzr(i, 1:nz) + sgnh*strainzz_udsh(1:nz)*strainzz_udrh(1:nz))
                        grad_c55(i, 1:nz) = grad_c55(i, 1:nz) &
                            - (strainxz(i, 1:nz)*strainxzr(i, 1:nz) + sgnh*strainxz_udsh(1:nz)*strainxz_udrh(1:nz))

                    end do
                    !$omp end parallel do

                end if

            end if

            if (yn_energy_precond) then

                !$omp parallel do private(i, j) collapse(2) schedule(auto)
                do j = 1, nz
                    do i = 1, nx
                        energy_src_v(i, j) = energy_src_v(i, j) + strainxx(i, j)**2 + strainzz(i, j)**2 + 2*strainxz(i, j)**2
                        energy_rec_v(i, j) = energy_rec_v(i, j) + strainxxr(i, j)**2 + strainzzr(i, j)**2 + 2*strainxzr(i, j)**2
                    end do
                end do
                !$omp end parallel do

            end if

        end if

        if (kernel_a /= '') then

            !$omp parallel do private(i, j) collapse(2) schedule(auto)
            do j = -pml + 1, nz + pml
                do i = -pml + 1, nx + pml
                    src_vx(i, j) = sum(vx(i:i + 1, j)) - sum(prev_vx(i:i + 1, j))
                    rec_vx(i, j) = sum(vxr(i:i + 1, j))
                    src_vz(i, j) = sum(vz(i, j:j + 1)) - sum(prev_vz(i, j:j + 1))
                    rec_vz(i, j) = sum(vzr(i, j:j + 1))
                end do
            end do
            !$omp end parallel do

            if (kernel_a == 'full') then

                !$omp parallel do private(i, j) collapse(2) schedule(auto)
                do j = 1, nz
                    do i = 1, nx
                        grad_rho(i, j) = grad_rho(i, j) + src_vx(i, j)*rec_vx(i, j) + src_vz(i, j)*rec_vz(i, j)
                    end do
                end do
                !$omp end parallel do

            else

                ! Along x
                if (index(kernel_a, 'lowx') /= 0) then
                    sgnh = 1
                else if (index(kernel_a, 'highx') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    !$omp parallel do private(j, p_lrsh, p_lrrh) schedule(auto)
                    do j = 1, nz

                        p_lrsh = src_vx(-pml + 1:nx + pml, j)
                        p_lrrh = rec_vx(-pml + 1:nx + pml, j)
                        call hilbert_transform(p_lrsh)
                        call hilbert_transform(p_lrrh)

                        grad_rho(1:nx, j) = grad_rho(1:nx, j) &
                            + (src_vx(1:nx, j)*rec_vx(1:nx, j) + sgnh*p_lrsh(1:nx)*p_lrrh(1:nx))

                        p_lrsh = src_vz(-pml + 1:nx + pml, j)
                        p_lrrh = rec_vz(-pml + 1:nx + pml, j)
                        call hilbert_transform(p_lrsh)
                        call hilbert_transform(p_lrrh)

                        grad_rho(1:nx, j) = grad_rho(1:nx, j) &
                            + (src_vz(1:nx, j)*rec_vz(1:nx, j) + sgnh*p_lrsh(1:nx)*p_lrrh(1:nx))

                    end do
                    !$omp end parallel do

                end if

                ! Along z
                if (index(kernel_a, 'lowz') /= 0) then
                    sgnh = 1
                else if (index(kernel_a, 'highz') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    !$omp parallel do private(i, p_udsh, p_udrh) schedule(auto)
                    do i = 1, nx

                        p_udsh = src_vx(i, -pml + 1:nz + pml)
                        p_udrh = rec_vx(i, -pml + 1:nz + pml)
                        call hilbert_transform(p_udsh)
                        call hilbert_transform(p_udrh)
                        grad_rho(i, 1:nz) = grad_rho(i, 1:nz) &
                            + (src_vx(i, 1:nz)*rec_vx(i, 1:nz) + sgnh*p_udsh(1:nz)*p_udrh(1:nz))

                        p_udsh = src_vz(i, -pml + 1:nz + pml)
                        p_udrh = rec_vz(i, -pml + 1:nz + pml)
                        call hilbert_transform(p_udsh)
                        call hilbert_transform(p_udrh)
                        grad_rho(i, 1:nz) = grad_rho(i, 1:nz) &
                            + (src_vz(i, 1:nz)*rec_vz(i, 1:nz) + sgnh*p_udsh(1:nz)*p_udrh(1:nz))

                    end do
                    !$omp end parallel do

                end if

            end if

            if (yn_energy_precond) then

                !$omp parallel do private(i, j) collapse(2) schedule(auto)
                do j = 1, nz
                    do i = 1, nx
                        energy_src_a(i, j) = energy_src_a(i, j) + src_vx(i, j)**2 + src_vz(i, j)**2
                        energy_rec_a(i, j) = energy_rec_a(i, j) + rec_vx(i, j)**2 + rec_vz(i, j)**2
                    end do
                end do
                !$omp end parallel do

            end if

        end if

    end subroutine

    !
    !> Compute source parameter gradients
    !
    subroutine compute_gradient_source(t)

        integer, intent(in) :: t

        integer :: sgx, sgz, i, irx, irz

        do i = 1, sgmtr%ns

            sgx = sgmtr%srcr(i)%gx
            sgz = sgmtr%srcr(i)%gz
            do irz = -nkw, nkw
                do irx = -nkw, nkw
                    if (ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                        grad_mt(1) = grad_mt(1) - &
                            stressxxr(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_ix(irx) &
                            *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                        grad_mt(3) = grad_mt(3) - &
                            stresszzr(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_ix(irx) &
                            *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                    end if
                end do
            end do

            sgx = sgmtr%srcr(i)%hx
            sgz = sgmtr%srcr(i)%hz
            do irz = -nkw, nkw
                do irx = -nkw, nkw
                    if (ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                        grad_mt(5) = grad_mt(5) - &
                            stressxzr(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_hx(irx) &
                            *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                    end if
                end do
            end do

        end do

    end subroutine

end submodule
