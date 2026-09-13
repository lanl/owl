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


module regularization

    use libflit
    use mod_parameters
    use mod_utility
    use mod_model

    implicit none

    real :: tv_mu
    real :: tv_lambda1
    real :: tv_lambda2
    real :: tv_norm
    integer :: tv_niter
    real :: tikhonov_lambda
    real :: reg_smoothx
    real :: reg_smoothy
    real :: reg_smoothz
    type(andf_param) :: param
    character(len=1024) :: file_soraux, file_sorcoh

contains

    subroutine tgpv_denoise(u, m)

        real, dimension(:, :, :), intent(in) :: u
        real, dimension(:, :, :), intent(inout) :: m

        if (rankid == 0) then
            call warn(date_time_compact()//' tv_mu = '//num2str(tv_mu))
            call warn(date_time_compact()//' tv_lambda1 = '//num2str(tv_lambda1))
            call warn(date_time_compact()//' tv_lambda2 = '//num2str(tv_lambda2))
            call warn(date_time_compact()//' tv_norm = '//num2str(tv_norm))
            call warn(date_time_compact()//' tv_niter = '//num2str(tv_niter))
        end if

        call mpibarrier

        m = tgpv_filt_mpi(u, tv_mu, tv_lambda1, tv_lambda2, tv_niter, tv_norm)

        call mpibarrier

    end subroutine tgpv_denoise

    !
    !> Solve the Tikhonov regularization in a hard way:
    !> Decompose the inversion using ADMM, and then solve
    !>      ||\nabla u||_2^2 = ||u-m^{k+1}||_2^2.
    !> The first-order optimality condition of the above equation is
    !>      \nabla^T \nabla u = u - m^{k+1},
    !> which can be solved using the Gauss-Seidel method.
    !
    subroutine tikhonov_denoise(m, u)

        real, dimension(:, :, :), intent(in) :: m
        real, dimension(:, :, :), intent(out) :: u

        real, allocatable, dimension(:, :, :) :: pu
        integer :: i, j, k, iter
        real :: sumu
        integer :: n1, n2, n3, niter
        integer :: i1, i2, j1, j2, k1, k2
        integer :: ii1, ii2, jj1, jj2, kk1, kk2
        real :: mu

        u = m
        mu = 1.0
        niter = 200
        n1 = size(u, 1)
        n2 = size(u, 2)
        n3 = size(u, 3)

        ! TGpV iteration
        do iter = 1, niter

            pu = u

            ! minimization u
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        if (i == 1) then
                            i1 = 0
                            ii1 = i
                        else
                            i1 = 1
                            ii1 = i - 1
                        end if

                        if (i == n1) then
                            i2 = 0
                            ii2 = i
                        else
                            i2 = 1
                            ii2 = i + 1
                        end if

                        if (j == 1) then
                            j1 = 0
                            jj1 = j
                        else
                            j1 = 1
                            jj1 = j - 1
                        end if
                        if (j == n2) then
                            j2 = 0
                            jj2 = j
                        else
                            j2 = 1
                            jj2 = j + 1
                        end if

                        if (k == 1) then
                            k1 = 0
                            kk1 = k
                        else
                            k1 = 1
                            kk1 = k - 1
                        end if
                        if (k == n3) then
                            k2 = 0
                            kk2 = k
                        else
                            k2 = 1
                            kk2 = k + 1
                        end if

                        sumu = tikhonov_lambda*( &
                            +i2*u(ii2, j, k) + i1*u(ii1, j, k) &
                            + j2*u(i, jj2, k) + j1*u(i, jj1, k) &
                            + k2*u(i, j, kk2) + k1*u(i, j, kk1)) &
                            + mu*m(i, j, k)

                        u(i, j, k) = sumu/(mu + (6.0 - (1 - i1) - (1 - i2) - (1 - j1) - (1 - j2) - (1 - k1) - (1 - k2))*tikhonov_lambda)

                    end do
                end do
            end do

            ! progress
            if (rankid == 0 .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
                call warn(date_time_compact()//'>> Tikhonov iteration '//tidy(num2str(iter, '(i)')) &
                    //' of '//tidy(num2str(niter, '(i)')// &
                    ' relative norm2 diff = '//tidy(num2str(norm2(u - pu), '(es)'))))
            end if

        end do

    end subroutine tikhonov_denoise

    !
    !> Apply a series of regularization to a single model parameter
    !
    subroutine regularize_single_parameter(m, mr, name, default_value, smooth_inverse)

        real, dimension(:, :, :), intent(in) :: m
        real, dimension(:, :, :), intent(inout) :: mr
        character(len=*), intent(in) :: name
        real, intent(in) :: default_value
        logical, intent(in), optional :: smooth_inverse

        integer :: i, n1, n2, n3
        real, allocatable, dimension(:, :, :) :: mt, maux, mcoh

        n1 = size(m, 1)
        n2 = size(m, 2)
        n3 = size(m, 3)
        mt = m

        call readpar_int(file_parameter, 'rankx', rank3, 1)
        call readpar_int(file_parameter, 'ranky', rank2, 1)
        call readpar_int(file_parameter, 'rankz', rank1, 1)

        do i = 1, size(model_regularization_method)

            select case (model_regularization_method(i))

                case ('Tikhonov', 'tikhonov')
                    call readpar_xfloat(file_parameter, 'reg_tikhonov_lambda', tikhonov_lambda, 10.0, iter*1.0)
                    call tikhonov_denoise(mt, mr)
                    mt = mr
                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Tikhonov regularization finished. ')
                    end if

                case ('smooth')
                    call readpar_xfloat(file_parameter, 'reg_smooth_x', reg_smoothx, -1.0, iter*1.0)
                    if (reg_smoothx < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_x_'//tidy(name), reg_smoothx, 1.0*dx, iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_smooth_y', reg_smoothy, -1.0, iter*1.0)
                    if (reg_smoothy < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_y_'//tidy(name), reg_smoothy, 1.0*dy, iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_smooth_z', reg_smoothz, -1.0, iter*1.0)
                    if (reg_smoothz < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_z_'//tidy(name), reg_smoothz, 1.0*dz, iter*1.0)
                    end if
                    if (present(smooth_inverse) .and. smooth_inverse) then
                        mt = 1.0/gauss_filt(1.0/mt, [reg_smoothz/dz, reg_smoothy/dy, reg_smoothx/dx])
                    else
                        mt = gauss_filt(mt, [reg_smoothz/dz, reg_smoothy/dy, reg_smoothx/dx])
                    end if
                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Smoothing regularization finished. ')
                    end if

                case ('TGpV', 'tgpv')
                    call readpar_xfloat(file_parameter, 'reg_tv_mu_'//tidy(name), tv_mu, default_value, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_tv_lambda1', tv_lambda1, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_tv_lambda2', tv_lambda2, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_tv_norm', tv_norm, 0.5, iter*1.0)
                    call readpar_xint(file_parameter, 'reg_tv_niter', tv_niter, 50, iter*1.0)
                    call tgpv_denoise(mt, mr)
                    mt = mr
                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> TGpV regularization finished. ')
                    end if

                case ('structure')
                    call readpar_xfloat(file_parameter, 'reg_andf_alpha', param%lambda1, 0.001, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_beta', param%lambda2, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_gamma', param%lambda3, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_smooth_x', param%smooth3, 2.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_smooth_y', param%smooth2, 2.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_smooth_z', param%smooth1, 8.0, iter*1.0)
                    call readpar_xint(file_parameter, 'reg_andf_t', param%niter, 5, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_sigma', param%sigma, 10.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_andf_powerm', param%powerm, 4.0, iter*1.0)
                    call readpar_xstring(file_parameter, 'reg_andf_aux', file_soraux, '', iter*1.0)
                    mr = mt
                    if (file_soraux /= '' .and. file_sorcoh == '') then
                        maux = load(file_soraux, n1, n2, n3)
                        mr = andf_filt(mr, param, aux=maux)
                    else if (file_soraux == '' .and. file_sorcoh /= '') then
                        mcoh = load(file_sorcoh, n1, n2, n3)
                        mr = andf_filt(mr, param, acoh=mcoh)
                    else if (file_soraux /= '' .and. file_sorcoh /= '') then
                        maux = load(file_soraux, n1, n2, n3)
                        mcoh = load(file_sorcoh, n1, n2, n3)
                        mr = andf_filt(mr, param, aux=maux, acoh=mcoh)
                    else
                        mr = andf_filt(mr, param)
                    end if
                    mt = mr
                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Structure-oriented regularization finished. ')
                    end if

            end select

        end do

        call mpibarrier
        mr = mt

        ! Output regularized model
        if (rankid == 0) then
            call output_array(mr, dir_iter_model(iter)//'/reg_'//tidy(name)//'.bin')
        end if

    end subroutine

    !
    !> Apply regularization to model parameters
    !
    subroutine model_regularization

        integer :: i

        do i = 1, nmodel
            if (model_m(i)%name == 'vp' .or. model_m(i)%name == 'vs') then
                call regularize_single_parameter(model_m(i)%array, &
                    model_reg(i)%array, model_name(i), 100.0/(maxval(abs(model_m(i)%array)) + float_tiny), &
                    smooth_inverse=.true.)
            else
                call regularize_single_parameter(model_m(i)%array, &
                    model_reg(i)%array, model_name(i), 100.0/(maxval(abs(model_m(i)%array)) + float_tiny))
            end if
        end do

    end subroutine

end module
