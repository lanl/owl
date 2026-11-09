!
! © 2025-2026. Triad National Security, LLC. All rights reserved.
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


module vars

    use libflit
    use mod_parameters
    use mod_model
    use mod_utility

    implicit none

contains

    !
    !> Initialize regularization auxiliary variables
    !
    subroutine init_reg

        integer :: i

        if (model_regularization_method(1) /= '') then

            allocate (model_reg(1:nmodel))

            do i = 1, nmodel
                call alloc_array(model_reg(i)%array, [1, nz, 1, nx], source=model_m(i)%array)
            end do

        end if

    end subroutine init_reg

    !
    !> Input medium parameter models for inversion
    !
    subroutine prepare_model

        integer :: i
        character(len=1024) :: fname

        ! Model static
        if (nmodel_aux >= 1) then

            allocate (model_aux(1:nmodel_aux))
            do i = 1, nmodel_aux

                ! Models that will not be updated
                call readpar_string(file_parameter, 'file_'//tidy(model_name_aux(i)), fname, '', .true.)
                call prepare_model_single_parameter(model_aux(i)%array, model_name_aux(i), fname, update=.false.)

                if (model_aux(i)%name == 'rho' .and. any(model_aux(i)%array == 0)) then
                    call warn(' <prepare_model> Error: One or multiple rho = 0')
                    call mpibarrier
                    call mpistop
                end if

                ! Model name
                model_aux(i)%name = model_name_aux(i)

            end do

        end if

        ! Model to be updated
        allocate (model_m(1:nmodel))
        allocate (model_grad(1:nmodel))
        allocate (model_srch(1:nmodel))
        allocate (model_m_backup(1:nmodel))

        do i = 1, nmodel

            call readpar_string(file_parameter, 'file_'//tidy(model_name(i)), fname, '')
            call prepare_model_single_parameter(model_m(i)%array, model_name(i), fname, update=.true.)

            if (model_m(i)%name == 'rho' .and. any(model_m(i)%array == 0)) then
                call warn(' <prepare_model> Error: One or multiple rho = 0')
                call mpibarrier
                call mpistop
            end if

            model_m(i)%name = model_name(i)
            model_grad(i)%name = model_name(i)

            model_grad(i)%array = zeros_like(model_m(i)%array)
            model_srch(i)%array = zeros_like(model_m(i)%array)
            model_m_backup(i)%array = zeros_like(model_m(i)%array)

        end do

    end subroutine

    !
    !> Clip Vp or Vs to ensure that Vp/Vs ratio is within a reasonable range
    !
    subroutine clip_vpvsratio

        integer :: i
        real, allocatable, dimension(:, :) :: mvp, mvs, r
        logical :: vp_from_aux, vs_from_aux

        if (index(which_medium, 'acoustic') /= 0) then
            return
        end if

        vp_from_aux = .false.
        vs_from_aux = .false.

        ! Vp or Vs might come from static models
        do i = 1, nmodel_aux
            if (model_aux(i)%name == 'vp') then
                mvp = model_aux(i)%array
                vp_from_aux = .true.
            end if
            if (model_aux(i)%name == 'vs') then
                mvs = model_aux(i)%array
                vs_from_aux = .true.
            end if
        end do

        ! Vp or/and Vs might come from update models
        do i = 1, nmodel
            if (model_m(i)%name == 'vp') then
                mvp = model_m(i)%array
                vp_from_aux = .false.
            end if
            if (model_m(i)%name == 'vs') then
                mvs = model_m(i)%array
                vs_from_aux = .false.
            end if
        end do

        ! If Vp is static model but Vs is update model, or both Vp and Vs are update model
        ! Then use Vp/Vs ratio to constrain Vs
        if ((vp_from_aux .and. .not. vs_from_aux) .or. (.not. vp_from_aux .and. .not. vs_from_aux)) then

            do i = 1, nmodel
                if (model_m(i)%name == 'vs') then
                    where (mvs /= 0)
                        mvs = clip(mvs, model_min(i), model_max(i))
                    end where
                end if
            end do

            r = clip(mvp/mvs, min_vpvsratio, max_vpvsratio)
            where (mvs == 0)
                r = 0
            end where

            if (vpvsratio_smoothx /= 0 .or. vpvsratio_smoothz /= 0) then
                r = smooth_vpvsratio(r, mvp /= 0 .and. mvs /= 0)
            end if
            r = clip(r, min_vpvsratio, max_vpvsratio)

            do i = 1, nmodel
                if (model_m(i)%name == 'vs') then
                    model_m(i)%array = mvp/r
                    where (mvs == 0)
                        model_m(i)%array = 0
                    end where
                end if
            end do

        end if

        ! If Vs is static but Vp is update,
        ! Then use Vp/Vs ratio to contrain Vp
        if (vs_from_aux .and. .not. vp_from_aux) then

            do i = 1, nmodel
                if (model_m(i)%name == 'vp') then
                    where (mvp /= 0)
                        mvp = clip(mvp, model_min(i), model_max(i))
                    end where
                end if
            end do

            r = clip(mvp/mvs, min_vpvsratio, max_vpvsratio)
            where (mvp == 0)
                r = 0
            end where

            if (vpvsratio_smoothx /= 0 .or. vpvsratio_smoothz /= 0) then
                r = smooth_vpvsratio(r, mvp /= 0 .and. mvs /= 0)
            end if
            r = clip(r, min_vpvsratio, max_vpvsratio)

            do i = 1, nmodel
                if (model_m(i)%name == 'vp') then
                    model_m(i)%array = mvs*r
                    where (mvp == 0)
                        model_m(i)%array = 0
                    end where
                    ! Where Vs is zero, e.g., in water, there is no ratio
                    ! to constrain Vp, which keeps its value
                    where (mvs == 0)
                        model_m(i)%array = mvp
                    end where
                end if
            end do

        end if

    end subroutine

    !
    !> Gaussian smoothing of the Vp/Vs ratio r in clip_vpvsratio. Points where Vp or
    !> Vs is zero, e.g., in air or water, have no ratio (valid is false). When such
    !> points exist, the smoothing uses the ratio of the valid points alone,
    !> normalized by the smoothed indicator of the valid points, and leaves r
    !> unchanged at the other points
    !
    function smooth_vpvsratio(r, valid) result(rs)

        real, dimension(:, :), intent(in) :: r
        logical, dimension(:, :), intent(in) :: valid
        real, allocatable, dimension(:, :) :: rs

        real, allocatable, dimension(:, :) :: w

        if (all(valid)) then
            rs = gauss_filt(r, [vpvsratio_smoothz/dz, vpvsratio_smoothx/dx])
        else
            rs = gauss_filt(merge(r, 0.0, valid), [vpvsratio_smoothz/dz, vpvsratio_smoothx/dx])
            w = gauss_filt(merge(1.0, 0.0, valid), [vpvsratio_smoothz/dz, vpvsratio_smoothx/dx])
            where (valid)
                rs = rs/w
            elsewhere
                rs = r
            end where
        end if

    end function smooth_vpvsratio

end module

