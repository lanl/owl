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


submodule(elastic_tti_2d) elastic_tti_2d_fwi_gradient

    use libflit
    use elastic_tti_2d_vars
    use elastic_tti_2d_boundary_saving
    use elastic_tti_2d_wavefield
    use mod_anisotropy

    implicit none

#include '../../lib/macro_thomsen_2d.f90'
#include '../../lib/macro_alkhalifah_tsvankin_2d.f90'

contains

    module subroutine compute_adjoint(this)

        class(wave_solver_elastic_tti_2d), intent(inout) :: this

        integer :: l, ir, irx, irz, rgx, rgz
        type(grid2) :: grd
        integer :: i, j, t
        real :: amp1, amp2

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
        else
            call this%seis_vxr%init(nt=nt, dt=dt, nr=sgmtr%nr)
        end if

        if (yn_compz) then
            call this%seis_vzr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_z.su')
            call this%seis_vzr%resamp(nnt=nt, ddt=dt)
        else
            call this%seis_vzr%init(nt=nt, dt=dt, nr=sgmtr%nr)
        end if

        ! Elastic compliances
        call alloc_array(s11, [1, nx, 1, nz], pad=pml + 1)
        call alloc_array(s13, [1, nx, 1, nz], pad=pml + 1)
        call alloc_array(s15, [1, nx, 1, nz], pad=pml + 1)
        call alloc_array(s33, [1, nx, 1, nz], pad=pml + 1)
        call alloc_array(s35, [1, nx, 1, nz], pad=pml + 1)
        call alloc_array(s55, [1, nx, 1, nz], pad=pml + 1)
        !$omp parallel do private(i, j)
        do j = -pml, nz + pml + 1
            do i = -pml, nx + pml + 1
                call cij_to_sij(c11(i, j), c13(i, j), c15(i, j), c33(i, j), c35(i, j), c55(i, j), &
                    s11(i, j), s13(i, j), s15(i, j), s33(i, j), s35(i, j), s55(i, j))
            end do
        end do
        !$omp end parallel do

        call alloc_array(snapvx, [1, nx, 1, nz], pad=pml)
        call alloc_array(snapvz, [1, nx, 1, nz], pad=pml)

        grad_c11 = zeros(nx, nz)
        grad_c13 = zeros(nx, nz)
        grad_c15 = zeros(nx, nz)
        grad_c33 = zeros(nx, nz)
        grad_c35 = zeros(nx, nz)
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
        call alloc_array(prev_stressxx_ixiz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_stressxx_hxhz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_stresszz_ixiz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_stresszz_hxhz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_stressxz_ixiz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_stressxz_hxhz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_vx_hxiz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_vx_ixhz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_vz_hxiz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_vz_ixhz, [1, nx, 1, nz], pad=pml + fdhalf)

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

        ! prepare boundary saving
        call prepare_boundary_saving
        call open_boundary_saving

        l = np
        do t = nt, 1, -1

            if (yn_update_medium .and. t >= sgmtr%srcr(1)%hnt) then

                prev_stressxx_ixiz = stressxx_ixiz
                prev_stressxx_hxhz = stressxx_hxhz
                prev_stresszz_ixiz = stresszz_ixiz
                prev_stresszz_hxhz = stresszz_hxhz
                prev_stressxz_ixiz = stressxz_ixiz
                prev_stressxz_hxhz = stressxz_hxhz
                prev_vx_hxiz = vx_hxiz
                prev_vx_ixhz = vx_ixhz
                prev_vz_hxiz = vz_hxiz
                prev_vz_ixhz = vz_ixhz

                ! -------------- Forward wavefield reconstruction -----------------------
                if (yn_free_surface) then
                    call update_wavefield_free_surface(-dt, &
                        stressxx_ixiz, stresszz_ixiz, stressxz_ixiz, &
                        stressxx_hxhz, stresszz_hxhz, stressxz_hxhz, &
                        vx_ixhz, vz_ixhz, &
                        vx_hxiz, vz_hxiz, &
                        memory_pdxvx_hxhz, &
                        memory_pdxvz_hxhz, &
                        memory_pdxvx_ixiz, &
                        memory_pdxvz_ixiz, &
                        memory_pdxxx_hxiz, &
                        memory_pdxxz_hxiz, &
                        memory_pdxxx_ixhz, &
                        memory_pdxxz_ixhz, &
                        memory_pdzvx_hxhz, &
                        memory_pdzvz_hxhz, &
                        memory_pdzvx_ixiz, &
                        memory_pdzvz_ixiz, &
                        memory_pdzxx_hxiz, &
                        memory_pdzxz_hxiz, &
                        memory_pdzzz_hxiz, &
                        memory_pdzxx_ixhz, &
                        memory_pdzxz_ixhz, &
                        memory_pdzzz_ixhz)
                else
                    call update_wavefield(-dt, &
                        stressxx_ixiz, stresszz_ixiz, stressxz_ixiz, &
                        stressxx_hxhz, stresszz_hxhz, stressxz_hxhz, &
                        vx_ixhz, vz_ixhz, &
                        vx_hxiz, vz_hxiz, &
                        memory_pdxvx_hxhz, &
                        memory_pdxvz_hxhz, &
                        memory_pdxvx_ixiz, &
                        memory_pdxvz_ixiz, &
                        memory_pdxxx_hxiz, &
                        memory_pdxxz_hxiz, &
                        memory_pdxxx_ixhz, &
                        memory_pdxxz_ixhz, &
                        memory_pdzvx_hxhz, &
                        memory_pdzvz_hxhz, &
                        memory_pdzvx_ixiz, &
                        memory_pdzvz_ixiz, &
                        memory_pdzxz_hxiz, &
                        memory_pdzzz_hxiz, &
                        memory_pdzxz_ixhz, &
                        memory_pdzzz_ixhz)
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

                        if (yn_free_surface) then

                            call alloc_array(snapvx, [1, nx, 1, nz], pad=pml)
                            call alloc_array(snapvz, [1, nx, 1, nz], pad=pml)

                            !$omp parallel do private(i, j) collapse(2)
                            do j = -pml + 1, nz + pml
                                do i = -pml + 1, nx + pml
                                    snapvx(i, j) = 0.5*(0.5*sum(vx_hxiz(i:i + 1, j)) + 0.5*sum(vx_ixhz(i, j:j + 1)))
                                    snapvz(i, j) = 0.5*(0.5*sum(vz_hxiz(i:i + 1, j)) + 0.5*sum(vz_ixhz(i, j:j + 1)))
                                end do
                            end do
                            !$omp end parallel do

                            open(3, file=tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_x_' &
                                //num2str(l)//'.txt')
                            do i = -pml + 1, nx + pml
                                do j = 1, nz + pml
                                    write(3, *) (i - 1)*dx, zz_i(i, j), snapvx(i, j)
                                end do
                            end do
                            close(3)

                            open(3, file=tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_z_' &
                                //num2str(l)//'.txt')
                            do i = -pml + 1, nx + pml
                                do j = 1, nz + pml
                                    write(3, *) (i - 1)*dx, zz_i(i, j), snapvz(i, j)
                                end do
                            end do
                            close(3)

                            call map_irregular_to_regular(snapvx, this, [1, this%nx, 1, this%nz])
                            call map_irregular_to_regular(snapvz, this, [1, this%nx, 1, this%nz])

                            call output_array(snapvx, tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_x_' &
                                //num2str(l)//'.bin', transp=.true.)
                            call output_array(snapvz, tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_z_' &
                                //num2str(l)//'.bin', transp=.true.)

                        else

                            !$omp parallel do private(i, j) collapse(2)
                            do j = -pml + 1, nz + pml - 1
                                do i = -pml + 1, nx + pml - 1
                                    snapvx(i, j) = 0.5*(0.5*sum(vx_hxiz(i:i + 1, j)) + 0.5*sum(vx_ixhz(i, j:j + 1)))
                                    snapvz(i, j) = 0.5*(0.5*sum(vz_hxiz(i:i + 1, j)) + 0.5*sum(vz_ixhz(i, j:j + 1)))
                                end do
                            end do

                            call output_array(snapvx(1:nx, 1:nz), tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_x_' &
                                //num2str(l)//'.bin', transp=.true.)
                            call output_array(snapvz(1:nx, 1:nz), tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_reconstructed_wavefield_z_' &
                                //num2str(l)//'.bin', transp=.true.)

                        end if


                        l = l - 1
                    end if
                end if

            end if

            ! -------------- Adjoint wavefield reverse-time propagation -----------------------
            if (yn_free_surface) then
                call update_wavefield_free_surface(-dt, &
                    stressxxr_ixiz, stresszzr_ixiz, stressxzr_ixiz, &
                    stressxxr_hxhz, stresszzr_hxhz, stressxzr_hxhz, &
                    vxr_ixhz, vzr_ixhz, &
                    vxr_hxiz, vzr_hxiz, &
                    memory_pdxvxr_hxhz, &
                    memory_pdxvzr_hxhz, &
                    memory_pdxvxr_ixiz, &
                    memory_pdxvzr_ixiz, &
                    memory_pdxxxr_hxiz, &
                    memory_pdxxzr_hxiz, &
                    memory_pdxxxr_ixhz, &
                    memory_pdxxzr_ixhz, &
                    memory_pdzvxr_hxhz, &
                    memory_pdzvzr_hxhz, &
                    memory_pdzvxr_ixiz, &
                    memory_pdzvzr_ixiz, &
                    memory_pdzxxr_hxiz, &
                    memory_pdzxzr_hxiz, &
                    memory_pdzzzr_hxiz, &
                    memory_pdzxxr_ixhz, &
                    memory_pdzxzr_ixhz, &
                    memory_pdzzzr_ixhz)
            else
                call update_wavefield(-dt, &
                    stressxxr_ixiz, stresszzr_ixiz, stressxzr_ixiz, &
                    stressxxr_hxhz, stresszzr_hxhz, stressxzr_hxhz, &
                    vxr_ixhz, vzr_ixhz, &
                    vxr_hxiz, vzr_hxiz, &
                    memory_pdxvxr_hxhz, &
                    memory_pdxvzr_hxhz, &
                    memory_pdxvxr_ixiz, &
                    memory_pdxvzr_ixiz, &
                    memory_pdxxxr_hxiz, &
                    memory_pdxxzr_hxiz, &
                    memory_pdxxxr_ixhz, &
                    memory_pdxxzr_ixhz, &
                    memory_pdzvxr_hxhz, &
                    memory_pdzvzr_hxhz, &
                    memory_pdzvxr_ixiz, &
                    memory_pdzvzr_ixiz, &
                    memory_pdzxzr_hxiz, &
                    memory_pdzzzr_hxiz, &
                    memory_pdzxzr_ixhz, &
                    memory_pdzzzr_ixhz)
            end if

            ! Read and add adjoint source
            !$omp parallel private(ir, irx, irz, rgx, rgz, amp1, amp2)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    rgx = sgmtr%recr(ir)%hx
                    rgz = sgmtr%recr(ir)%gz
                    !$omp do collapse(2)
                    do irz = -nkw, nkw
                        do irx = -nkw, nkw
                            if (ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then
                                amp1 = sgmtr%recr(ir)%interp_hx(irx) &
                                    *sgmtr%recr(ir)%interp_iz(irz) &
                                    *sgmtr%recr(ir)%weight*0.5
                                vxr_hxiz(rgx + irx, rgz + irz) = vxr_hxiz(rgx + irx, rgz + irz) &
                                    + this%seis_vxr%trace(ir)%data(t)*amp1
                                vzr_hxiz(rgx + irx, rgz + irz) = vzr_hxiz(rgx + irx, rgz + irz) &
                                    + this%seis_vzr%trace(ir)%data(t)*amp1
                            end if
                        end do
                    end do
                    !$omp end do

                    rgx = sgmtr%recr(ir)%gx
                    rgz = sgmtr%recr(ir)%hz
                    !$omp do collapse(2)
                    do irz = -nkw, nkw
                        do irx = -nkw, nkw
                            if (ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then
                                amp2 = sgmtr%recr(ir)%interp_ix(irx) &
                                    *sgmtr%recr(ir)%interp_hz(irz) &
                                    *sgmtr%recr(ir)%weight*0.5
                                vxr_ixhz(rgx + irx, rgz + irz) = vxr_ixhz(rgx + irx, rgz + irz) &
                                    + this%seis_vxr%trace(ir)%data(t)*amp2
                                vzr_ixhz(rgx + irx, rgz + irz) = vzr_ixhz(rgx + irx, rgz + irz) &
                                    + this%seis_vzr%trace(ir)%data(t)*amp2
                            end if
                        end do
                    end do
                    !$omp end do

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
                if (any(isnan(vxr_ixhz)) .or. any(isnan(vzr_ixhz))) then
                    call warn(date_time_compact()//' >> Vxr, Vzr contain NaN!')
                    stop
                else
                    call warn(date_time_compact()//' >> Vxr, Vzr value range = ')
                    call warn(date_time_compact()//'      '//num2str(min(minval(vxr_hxiz), minval(vxr_ixhz)), '(es)') &
                        //' ~ '//num2str(max(maxval(vxr_hxiz), maxval(vxr_ixhz)), '(es)'))
                    call warn(date_time_compact()//'      '//num2str(min(minval(vzr_hxiz), minval(vzr_ixhz)), '(es)') &
                        //' ~ '//num2str(max(maxval(vzr_hxiz), maxval(vzr_ixhz)), '(es)'))
                end if
            end if

        end do

        ! Delete temporary files
        call close_boundary_saving(delete=.true.)

        ! Output source parameter gradient
        if (yn_update_source) then
            call grd%init(n=[nc_mt, 1], d=[1.0, 1.0], o=[0.0, 0.0])
            grd%array = -reshape(grad_mt/maxval(abs(grad_mt)), [nc_mt, 1])
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
                grad_c15 = grad_c15/energy_src_v
                grad_c33 = grad_c33/energy_src_v
                grad_c35 = grad_c35/energy_src_v
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
                call map_irregular_to_regular(grad_c15, this, [1, nx, 1, nz])
                call map_irregular_to_regular(grad_c33, this, [1, nx, 1, nz])
                call map_irregular_to_regular(grad_c35, this, [1, nx, 1, nz])
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

        ! Compute for medium parameter gradients using first-order optimization
        select case (aniso_param)

            case ('iso')

                vp = vp(1:nx, 1:nz)*1.0e-3
                vs = vs(1:nx, 1:nz)*1.0e-3

                ! Using chain rule to compute the gradients of Vp and Vs from the gradients of Cij
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
                    + thomsen_dc15_dvp*grad_c15 &
                    + thomsen_dc33_dvp*grad_c33 &
                    + thomsen_dc35_dvp*grad_c35 &
                    + thomsen_dc55_dvp*grad_c55

                grad_vs = thomsen_dc11_dvs*grad_c11 &
                    + thomsen_dc13_dvs*grad_c13 &
                    + thomsen_dc15_dvs*grad_c15 &
                    + thomsen_dc33_dvs*grad_c33 &
                    + thomsen_dc35_dvs*grad_c35 &
                    + thomsen_dc55_dvs*grad_c55

                grad_epsilon = thomsen_dc11_deps*grad_c11 &
                    + thomsen_dc13_deps*grad_c13 &
                    + thomsen_dc15_deps*grad_c15 &
                    + thomsen_dc33_deps*grad_c33 &
                    + thomsen_dc35_deps*grad_c35 &
                    + thomsen_dc55_deps*grad_c55

                grad_delta = thomsen_dc11_ddel*grad_c11 &
                    + thomsen_dc13_ddel*grad_c13 &
                    + thomsen_dc15_ddel*grad_c15 &
                    + thomsen_dc33_ddel*grad_c33 &
                    + thomsen_dc35_ddel*grad_c35 &
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
                    + alkhalifah_tsvankin_dc15_dvp*grad_c15 &
                    + alkhalifah_tsvankin_dc33_dvp*grad_c33 &
                    + alkhalifah_tsvankin_dc35_dvp*grad_c35 &
                    + alkhalifah_tsvankin_dc55_dvp*grad_c55

                grad_vs = alkhalifah_tsvankin_dc11_dvs*grad_c11 &
                    + alkhalifah_tsvankin_dc13_dvs*grad_c13 &
                    + alkhalifah_tsvankin_dc15_dvs*grad_c15 &
                    + alkhalifah_tsvankin_dc33_dvs*grad_c33 &
                    + alkhalifah_tsvankin_dc35_dvs*grad_c35 &
                    + alkhalifah_tsvankin_dc55_dvs*grad_c55

                grad_epsilon = alkhalifah_tsvankin_dc11_deps*grad_c11 &
                    + alkhalifah_tsvankin_dc13_deps*grad_c13 &
                    + alkhalifah_tsvankin_dc15_deps*grad_c15 &
                    + alkhalifah_tsvankin_dc33_deps*grad_c33 &
                    + alkhalifah_tsvankin_dc35_deps*grad_c35 &
                    + alkhalifah_tsvankin_dc55_deps*grad_c55

                grad_eta = alkhalifah_tsvankin_dc11_deta*grad_c11 &
                    + alkhalifah_tsvankin_dc13_deta*grad_c13 &
                    + alkhalifah_tsvankin_dc15_deta*grad_c15 &
                    + alkhalifah_tsvankin_dc33_deta*grad_c33 &
                    + alkhalifah_tsvankin_dc35_deta*grad_c35 &
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

                grd%array = transpose(return_normal(grad_c15))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c15.grd')

                grd%array = transpose(return_normal(grad_c33))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c33.grd')

                grd%array = transpose(return_normal(grad_c35))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c35.grd')

                grd%array = transpose(return_normal(grad_c55))
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c55.grd')

        end select

        grd%array = transpose(return_normal(grad_rho))
        call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_rho.grd')

        call mpibarrier_group

    end subroutine

    subroutine compute_gradient

        integer :: i, j
        integer :: sgnh
        real :: tmpxx, tmpzz, tmpxz, tmpxxr, tmpzzr, tmpxzr

        !$omp parallel do private(i, j, tmpxx, tmpzz, tmpxz, tmpxxr, tmpzzr, tmpxzr) collapse(2) schedule(auto)
        do j = -pml + 1, nz + pml
            do i = -pml + 1, nx + pml

                tmpxx = (stressxx_ixiz(i, j) + 0.25*sum(stressxx_hxhz(i:i + 1, j:j + 1))) &
                    - (prev_stressxx_ixiz(i, j) + 0.25*sum(prev_stressxx_hxhz(i:i + 1, j:j + 1)))
                tmpzz = (stresszz_ixiz(i, j) + 0.25*sum(stresszz_hxhz(i:i + 1, j:j + 1))) &
                    - (prev_stresszz_ixiz(i, j) + 0.25*sum(prev_stresszz_hxhz(i:i + 1, j:j + 1)))
                tmpxz = (stressxz_ixiz(i, j) + 0.25*sum(stressxz_hxhz(i:i + 1, j:j + 1))) &
                    - (prev_stressxz_ixiz(i, j) + 0.25*sum(prev_stressxz_hxhz(i:i + 1, j:j + 1)))

                tmpxxr = stressxxr_ixiz(i, j) + 0.25*sum(stressxxr_hxhz(i:i + 1, j:j + 1))
                tmpzzr = stresszzr_ixiz(i, j) + 0.25*sum(stresszzr_hxhz(i:i + 1, j:j + 1))
                tmpxzr = stressxzr_ixiz(i, j) + 0.25*sum(stressxzr_hxhz(i:i + 1, j:j + 1))

                strainxx(i, j) = s11(i, j)*tmpxx + s13(i, j)*tmpzz + s15(i, j)*tmpxz
                strainzz(i, j) = s13(i, j)*tmpxx + s33(i, j)*tmpzz + s35(i, j)*tmpxz
                strainxz(i, j) = s15(i, j)*tmpxx + s35(i, j)*tmpzz + s55(i, j)*tmpxz
                strainxxr(i, j) = s11(i, j)*tmpxxr + s13(i, j)*tmpzzr + s15(i, j)*tmpxzr
                strainzzr(i, j) = s13(i, j)*tmpxxr + s33(i, j)*tmpzzr + s35(i, j)*tmpxzr
                strainxzr(i, j) = s15(i, j)*tmpxxr + s35(i, j)*tmpzzr + s55(i, j)*tmpxzr

            end do
        end do
        !$omp end parallel do

        if (kernel_v /= '') then

            if (kernel_v == 'full') then

                !$omp parallel do private(i, j) collapse(2) schedule(auto)
                do j = 1, nz
                    do i = 1, nx

                        grad_c11(i, j) = grad_c11(i, j) - (strainxxr(i, j)*strainxx(i, j))
                        grad_c13(i, j) = grad_c13(i, j) - (strainxxr(i, j)*strainzz(i, j) &
                            + strainzzr(i, j)*strainxx(i, j))
                        grad_c15(i, j) = grad_c15(i, j) - (strainxxr(i, j)*strainxz(i, j) &
                            + strainxx(i, j)*strainxzr(i, j))
                        grad_c33(i, j) = grad_c33(i, j) - (strainzzr(i, j)*strainzz(i, j))
                        grad_c35(i, j) = grad_c35(i, j) - (strainzzr(i, j)*strainxz(i, j) &
                            + strainzz(i, j)*strainxzr(i, j))
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
                        grad_c15(1:nx, j) = grad_c15(1:nx, j) &
                            - (strainxx(1:nx, j)*strainxzr(1:nx, j) + sgnh*strainxx_lrsh(1:nx)*strainxz_lrrh(1:nx))
                        grad_c33(1:nx, j) = grad_c33(1:nx, j) &
                            - (strainzz(1:nx, j)*strainzzr(1:nx, j) + sgnh*strainzz_lrsh(1:nx)*strainzz_lrrh(1:nx))
                        grad_c35(1:nx, j) = grad_c35(1:nx, j) &
                            - (strainzz(1:nx, j)*strainxzr(1:nx, j) + sgnh*strainzz_lrsh(1:nx)*strainxz_lrrh(1:nx))
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
                        grad_c15(i, 1:nz) = grad_c15(i, 1:nz) &
                            - (strainxx(i, 1:nz)*strainxzr(i, 1:nz) + sgnh*strainxx_udsh(1:nz)*strainxz_udrh(1:nz))
                        grad_c33(i, 1:nz) = grad_c33(i, 1:nz) &
                            - (strainzz(i, 1:nz)*strainzzr(i, 1:nz) + sgnh*strainzz_udsh(1:nz)*strainzz_udrh(1:nz))
                        grad_c35(i, 1:nz) = grad_c35(i, 1:nz) &
                            - (strainzz(i, 1:nz)*strainxzr(i, 1:nz) + sgnh*strainzz_udsh(1:nz)*strainxz_udrh(1:nz))
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

                    src_vx(i, j) = (0.5*sum(vx_hxiz(i:i + 1, j)) + 0.5*sum(vx_ixhz(i, j:j + 1))) &
                        - (0.5*sum(prev_vx_hxiz(i:i + 1, j)) + 0.5*sum(prev_vx_ixhz(i, j:j + 1)))
                    rec_vx(i, j) = 0.5*sum(vxr_hxiz(i:i + 1, j)) + 0.5*sum(vxr_ixhz(i, j:j + 1))

                    src_vz(i, j) = (0.5*sum(vz_hxiz(i:i + 1, j)) + 0.5*sum(vz_ixhz(i, j:j + 1))) &
                        - (0.5*sum(prev_vz_hxiz(i:i + 1, j)) + 0.5*sum(prev_vz_ixhz(i, j:j + 1)))
                    rec_vz(i, j) = 0.5*sum(vzr_hxiz(i:i + 1, j)) + 0.5*sum(vzr_ixhz(i, j:j + 1))

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

            ! ix-iz
            sgx = sgmtr%srcr(i)%gx
            sgz = sgmtr%srcr(i)%gz
            do irz = -nkw, nkw
                do irx = -nkw, nkw
                    if (ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                        grad_mt(1) = grad_mt(1) - &
                            stressxxr_ixiz(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_ix(irx) &
                            *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                        grad_mt(3) = grad_mt(3) - &
                            stresszzr_ixiz(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_ix(irx) &
                            *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                        grad_mt(5) = grad_mt(5) - &
                            stressxzr_ixiz(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_ix(irx) &
                            *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                    end if
                end do
            end do

            ! hx-hz
            sgx = sgmtr%srcr(i)%hx
            sgz = sgmtr%srcr(i)%hz
            do irz = -nkw, nkw
                do irx = -nkw, nkw
                    if (ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                        grad_mt(1) = grad_mt(1) - &
                            stressxxr_hxhz(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_hx(irx) &
                            *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                        grad_mt(3) = grad_mt(3) - &
                            stresszzr_hxhz(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_hx(irx) &
                            *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                        grad_mt(5) = grad_mt(5) - &
                            stressxzr_hxhz(sgx + irx, sgz + irz) &
                            *sgmtr%srcr(i)%interp_hx(irx) &
                            *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                    end if
                end do
            end do

        end do

    end subroutine

end submodule
