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


module mod_data_processing

    use libflit
    use mod_parameters
    use mod_su
    use mod_utility
    use mod_model
    use mod_source_receiver
    use, intrinsic :: ieee_arithmetic

    implicit none

    private :: top_mute, bottom_mute

contains

    !
    !> Top mute
    !
    subroutine top_mute(seis, srcindex)

        type(su), intent(inout) :: seis
        integer, intent(in) :: srcindex

        integer :: i, nr, nc
        real :: ds
        real, allocatable, dimension(:) :: tt
        integer :: nsample
        integer, allocatable, dimension(:) :: nca

        ! Dimension parameters
        nr = seis%nr
        ds = seis%dt
        nsample = seis%nt
        nca = zeros(nr)

        if (topmute_vel > 0) then
            ! if mute by near surface vp, which should be homogeneous

            !$omp parallel do private(i)
            do i = 1, nr
                if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                    nca(i) = nint((gmtr(set_gmtrid(srcindex))%recr(i)%aoff/topmute_vel &
                        + topmute_width + topmute_shift + gmtr(set_gmtrid(srcindex))%srcr(1)%t0)/ds)
                end if
            end do
            !$omp end parallel do

        end if

        if (dir_topmute_time /= '') then
            ! if mute by pre-computed traveltime, which should be fixed during iterations

            tt = load(tidy(dir_topmute_time)//'/shot_'//num2str(set_srcid(srcindex))//'_traveltime_p.bin', nr)

            !$omp parallel do private(i)
            do i = 1, nr
                if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                    ! Here the time delay should not be added any more, because
                    ! giving a FAT from file is like giving a speicific, "hard" traveltime value
                    ! for the mute to start
                    nca(i) = max(nint((tt(i) + topmute_width + topmute_shift)/ds), nca(i))
                end if
            end do
            !$omp end parallel do

        end if

        !$omp parallel do private(i, nc)
        do i = 1, nr
            if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then

                ! Mute
                if (topmute_const > 0) then
                    nc = min(nint(topmute_const/ds) + 1, nsample)
                else
                    nc = clip(nca(i), 1, nsample)
                end if
                seis%trace(i)%data(1:nc) = 0.0

                ! Taper
                if (topmute_taper > 0) then
                    seis%trace(i)%data(nc:) = taper(seis%trace(i)%data(nc:), len=[nint(topmute_taper/ds), 0], method=['hann', ''])
                end if

            else

                seis%trace(i)%data = 0.0

            end if
        end do
        !$omp end parallel do

    end subroutine

    !
    !> Bottom mute
    !
    subroutine bottom_mute(seis, srcindex)

        type(su), intent(inout) :: seis
        integer, intent(in) :: srcindex

        integer :: i, nr, nc
        real :: ds
        real, allocatable, dimension(:) :: tt
        integer :: nsample
        integer, allocatable, dimension(:) :: nca, ncb

        ! Dimension parameters
        nr = seis%nr
        ds = seis%dt
        nsample = seis%nt
        nca = zeros(nr)
        ncb = zeros(nr)

        if (btmmute_fromtop > 0) then
            ! If the bottom mute is modified from the top mute by shifting a constant value,
            ! then first compute the topmute start value

            if (topmute_vel > 0) then

                !$omp parallel do private(i)
                do i = 1, nr
                    if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                        nca(i) = nint((gmtr(set_gmtrid(srcindex))%recr(i)%aoff/topmute_vel &
                            + gmtr(set_gmtrid(srcindex))%srcr(1)%t0)/ds)
                    end if
                end do
                !$omp end parallel do

            end if

            if (dir_topmute_time /= '') then

                tt = load(tidy(dir_topmute_time)//'/shot_'//num2str(set_srcid(srcindex))//'_traveltime_p.bin', nr)

                !$omp parallel do private(i)
                do i = 1, nr
                    if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                        nca(i) = max(nint(tt(i)/ds), nca(i))
                    end if
                end do
                !$omp end parallel do

            end if

            nca = nca + nint(btmmute_fromtop/ds)

        end if

        if (btmmute_vel > 0) then
            ! If mute by a given constant surface vp value

            !$omp parallel do private(i) schedule(dynamic)
            do i = 1, nr
                if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                    ! Here the mute_width should be not used, because the bottom mute is
                    ! supposed to assign zeros to the samples after the ncb value
                    ncb(i) = nint((gmtr(set_gmtrid(srcindex))%recr(i)%aoff/btmmute_vel + &
                        gmtr(set_gmtrid(srcindex))%srcr(1)%t0)/ds)
                end if
            end do
            !$omp end parallel do

        end if

        if (dir_btmmute_time /= '') then
            ! if mute by pre-computed traveltime, which should be fixed during iterations

            tt = load(tidy(dir_btmmute_time)//'/shot_'//num2str(set_srcid(srcindex))//'_traveltime_p.bin', nr)

            !$omp parallel do private(i) schedule(dynamic)
            do i = 1, nr
                if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                    ! Same here, mute_width should not be used anymore
                    ! And choose the larger value
                    ncb(i) = max(nint(tt(i)/ds), ncb(i))
                end if
            end do
            !$omp end parallel do

        end if

        !$omp parallel do private(i, nc) schedule(dynamic)
        do i = 1, nr
            if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                ! Muting
                if (btmmute_const > 0) then
                    nc = min(nint(btmmute_const/ds) + 1, nsample)
                else
                    nc = clip(max(nca(i) + 1, ncb(i) + 1), 1, nsample)
                end if
                seis%trace(i)%data(nc + 1:) = 0.0
                ! Tapering
                if (btmmute_taper > 0) then
                    seis%trace(i)%data(:nc) = taper(seis%trace(i)%data(:nc), len=[0, nint(btmmute_taper/ds)], method=['', 'hann'])
                end if
            else
                seis%trace(i)%data = 0.0
            end if
        end do
        !$omp end parallel do

    end subroutine

    !
    !> Surgical mute
    !
    subroutine surgical_mute(seis, srcindex)

        type(su), intent(inout) :: seis
        integer, intent(in) :: srcindex

        integer :: i, nr
        real :: ds
        real, allocatable, dimension(:) :: tt, tp
        integer, allocatable, dimension(:) :: nca, ncb

        ! Dimension parameters
        nr = seis%nr
        ds = seis%dt
        nca = zeros(nr)
        ncb = zeros(nr)

        if (dir_surgicalmute_time /= '') then
            ! if mute by pre-computed traveltime, which should be fixed during iterations

            tt = load(tidy(dir_surgicalmute_time)//'/shot_'//num2str(set_srcid(srcindex))//'_traveltime_p.bin', 2*nr)

            !$omp parallel do private(i) schedule(dynamic)
            do i = 1, nr
                if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                    nca(i) = nint(tt(i)/ds)
                    ncb(i) = nint(tt(i + nr)/ds)
                end if
            end do
            !$omp end parallel do

        else

            !$omp parallel do private(i) schedule(dynamic)
            do i = 1, nr
                if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                    nca(i) = nint((gmtr(set_gmtrid(srcindex))%recr(i)%aoff/surgicalmute_vel(1) + &
                        gmtr(set_gmtrid(srcindex))%srcr(1)%t0)/ds)
                    ncb(i) = nint((gmtr(set_gmtrid(srcindex))%recr(i)%aoff/surgicalmute_vel(2) + &
                        gmtr(set_gmtrid(srcindex))%srcr(1)%t0 + surgicalmute_width)/ds)
                end if
            end do
            !$omp end parallel do

        end if

        ! Taper the trace
        tp = ones(seis%nt)
        !$omp parallel do private(i, tp) schedule(dynamic)
        do i = 1, nr
            if (gmtr(set_gmtrid(srcindex))%recr(i)%weight /= 0) then
                tp = 1.0
                tp = taper(tp, len=nint(surgicalmute_taper/ds), protect=[nca(i), ncb(i)], method=['hann', 'hann'])
                if (surgicalmute_inverse) then
                    seis%trace(i)%data = tp*seis%trace(i)%data
                else
                    seis%trace(i)%data = (1 - tp)*seis%trace(i)%data
                end if
            else
                seis%trace(i)%data = 0.0
            end if
        end do
        !$omp end parallel do

    end subroutine

    !
    !> Correct synthetic data source time function to match record
    !
    subroutine stf_correction(seis_syn_, seis_obs_)

        real, dimension(:, :), intent(inout) :: seis_syn_
        real, dimension(:, :), intent(in) :: seis_obs_

        integer :: i, n, nv, nr
        real, allocatable, dimension(:, :) :: seis_obs, seis_syn
        real, allocatable, dimension(:) :: mf
        integer :: ief
        complex, allocatable, dimension(:) :: vs, vd
        complex, allocatable, dimension(:, :) :: ms, md
        real :: eps

        n = size(seis_syn_, 1)
        nr = size(seis_syn_, 2)

        nv = next_power_235(2*n)
        ms = zeros(nv, nr)
        md = zeros(nv, nr)

        seis_syn = pad(taper(seis_syn_, [0, nint(0.05*n), 0, 0]), [0, nv - n, 0, 0])
        seis_obs = pad(taper(seis_obs_, [0, nint(0.05*n), 0, 0]), [0, nv - n, 0, 0])

        !$omp parallel do private(i)
        do i = 1, nr
            ms(:, i) = fft(seis_syn(:, i))
            md(:, i) = fft(seis_obs(:, i))
        end do
        !$omp end parallel do

        vs = zeros(nv)
        vd = zeros(nv)
        !$omp parallel do private(i) reduction(+: vs) reduction(+: vd)
        do i = 1, nr
            vs = vs + md(:, i)*conjg(ms(:, i))
            vd = vd + ms(:, i)*conjg(ms(:, i))
        end do
        !$omp end parallel do

        eps = 0.001*maxval(abs(vd))

        mf = ifft(vs/(vd + eps), real=.true.)
        mf = fftshift(mf)

        ief = maxloc(envelope(mf), dim=1)
        nv = nint(0.1*n)
        mf = mf(ief - nv:ief + nv)

        !$omp parallel do private(i)
        do i = 1, nr
            seis_syn(:, i) = conv(seis_syn(:, i), mf, 'same')
        end do
        !$omp end parallel do

        seis_syn_ = seis_syn(1:n, :)

    end subroutine

    !
    !> Processing single component
    !
    subroutine process_single_component(dirin, dirout, srcindex, component, processings)

        character(len=*), intent(in) :: dirin, dirout, component
        character(len=*), dimension(:), intent(in) :: processings
        integer, intent(in) :: srcindex

        integer :: i, j, t
        real :: maxoff
        real, allocatable, dimension(:, :) :: w
        integer :: wn1, wn2
        type(andf_param) :: param
        type(su) :: seis, dseis, seisx, seisy, seisz
        real, allocatable, dimension(:) :: p
        real, allocatable, dimension(:, :) :: u, d
        character(len=1024) :: dir_processed

        ! Original data
        seis%nr = gmtr(set_gmtrid(srcindex))%nr
        call seis%load(tidy(dirin)//'/shot_'//num2str(set_srcid(srcindex))//'_seismogram_'//tidy(component)//'.su')
        call seis%resamp(nint(tmax/seis%dt + 1), seis%dt)

        ! Data processing
        do i = 1, size(processings)

            select case (processings(i))

                case ('stf_correction')
                    ! Synthetic data
                    u = seis%to_array()

                    ! Observed data
                    if (record_processing(1) /= '') then
                        dir_processed = tidy(dir_working)//'/record_processed'
                    else
                        dir_processed = tidy(dir_record)
                    end if

                    call dseis%load(tidy(dir_processed)//'/shot_'//num2str(set_srcid(srcindex)) &
                        //'_seismogram_'//tidy(component)//'.su')
                    call dseis%resamp(nint(tmax/seis%dt + 1), seis%dt)
                    d = dseis%to_array()

                    ! STF correction with deconvolution and convolution
                    call stf_correction(u, d)

                    call seis%from_array(u)

                case ('time_shift')
                    call seis%shift(tshift)

                case ('time_deriv')
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        seis%trace(j)%data = deriv(seis%trace(j)%data)/dt
                    end do
                    !$omp end parallel do

                case ('time_integ')
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        seis%trace(j)%data = integ(seis%trace(j)%data)*dt
                    end do
                    !$omp end parallel do

                case ('top_mute')
                    call top_mute(seis, srcindex)

                case ('bottom_mute')
                    call bottom_mute(seis, srcindex)

                case ('surgical_mute')
                    call surgical_mute(seis, srcindex)

                case ('subtract_base')
                    call dseis%load(tidy(dir_base)//'/shot_'//num2str(set_srcid(srcindex)) &
                        //'_seismogram_'//tidy(component)//'.su')
                    seis = seis - dseis

                case ('t_power')
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        seis%trace(j)%data = seis%trace(j)%data*regspace(0.0, seis%dt, (seis%nt - 1)*seis%dt)**tpow
                    end do
                    !$omp end parallel do

                case ('t_balance')
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        if (maxval(abs(seis%trace(j)%data)) /= 0) then
                            seis%trace(j)%data = balance_filt(seis%trace(j)%data, &
                                nint(movingbalwin/seis%trace(j)%header%d1), 0.001)
                        end if
                    end do
                    !$omp end parallel do

                case ('freq_filt')
                    call seis%freqfilt(freqs, amps)

                case ('rms_balance')
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        seis%trace(j)%data = seis%trace(j)%data/(norm2(seis%trace(j)%data) + float_tiny)
                    end do
                    !$omp end parallel do

                case ('offset_balance')
                    maxoff = 0.0
                    do j = 1, seis%nr
                        if (gmtr(set_gmtrid(srcindex))%recr(j)%aoff > maxoff) then
                            maxoff = gmtr(set_gmtrid(srcindex))%recr(j)%aoff
                        end if
                    end do
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        seis%trace(j)%data = seis%trace(j)%data*(1.0 - &
                            gmtr(set_gmtrid(srcindex))%recr(j)%aoff/maxoff)
                    end do
                    !$omp end parallel do

                case ('max_balance')
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        if (maxval(abs(seis%trace(j)%data)) /= 0) then
                            seis%trace(j)%data = seis%trace(j)%data/maxval(abs(seis%trace(j)%data))
                        end if
                    end do
                    !$omp end parallel do

                case ('dip_filt')
                    w = seis%to_array()
                    w = dip_filt(w, [seis%dt, 1.0], dips, dipamps)
                    call seis%from_array(w)

                case ('andf_filt')
                    wn1 = seis%nt
                    wn2 = seis%nr
                    w = seis%to_array()
                    param = proc_andfparam
                    param%smooth1 = param%smooth1/seis%dt
                    w = andf_filt(w, param)
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        seis%trace(j)%data(:) = w(:, j)
                    end do
                    !$omp end parallel do

                case ('remove_nan')
                    !$omp parallel do private(j)
                    do j = 1, seis%nr
                        seis%trace(j)%data = return_normal(seis%trace(j)%data)
                    end do
                    !$omp end parallel do

                case ('rotate')
                    ! Rotate multi-component data
                    !   +: ccw
                    !   -: cw

                    seisx%nr = gmtr(set_gmtrid(srcindex))%nr
                    seisy%nr = gmtr(set_gmtrid(srcindex))%nr
                    seisz%nr = gmtr(set_gmtrid(srcindex))%nr

                    ! Read x component
                    call seisx%load(tidy(dirin)//'/shot_'//num2str(set_srcid(srcindex)) &
                        //'_seismogram_x.su')
                    call seisx%resamp(nnt=nint(tmax/seisx%dt + 1), ddt=seisx%dt)

                    ! Read y component
                    call seisy%load(tidy(dirin)//'/shot_'//num2str(set_srcid(srcindex)) &
                        //'_seismogram_y.su')
                    call seisy%resamp(nnt=nint(tmax/seisy%dt + 1), ddt=seisy%dt)

                    ! Read z component
                    call seisz%load(tidy(dirin)//'/shot_'//num2str(set_srcid(srcindex)) &
                        //'_seismogram_z.su')
                    call seisz%resamp(nnt=nint(tmax/seisz%dt + 1), ddt=seisz%dt)

                    !$omp parallel do private(j, t, p)
                    do j = 1, seis%nr
                        do t = 1, seis%nt
                            p = rotate_point([seisx%trace(j)%data(t), seisy%trace(j)%data(t), seisz%trace(j)%data(t)], &
                                [drotx, droty, drotz]*real(const_deg2rad), [0.0, 0.0, 0.0], 'xyz')
                            select case(component)
                                case('x')
                                    seis%trace(j)%data(t) = p(1)
                                case('y')
                                    seis%trace(j)%data(t) = p(2)
                                case('z')
                                    seis%trace(j)%data(t) = p(3)
                            end select
                        end do
                    end do
                    !$omp end parallel do

            end select

            if (rankid_group == 0) then
                call warn(date_time_compact()//' >> Shot '//num2str(set_srcid(srcindex)) &
                    //' data processing (' // tidy(processings(i)) // ') completed. ')
            end if

        end do

        ! Zero unrelevant traces and zero NaN and Inf
        !$omp parallel do private(i)
        do i = 1, seis%nr
            if (gmtr(set_gmtrid(srcindex))%recr(i)%weight == 0) then
                seis%trace(i)%data = 0.0
            end if
            where (.not. ieee_is_finite(seis%trace(i)%data) .or. ieee_is_nan(seis%trace(i)%data))
                seis%trace(i)%data = 0.0
            end where
        end do
        !$omp end parallel do

        ! Output processed data
        call seis%output(tidy(dirout)//'/shot_'//num2str(set_srcid(srcindex)) &
            //'_seismogram_'//tidy(component)//'.su')

    end subroutine

    !
    !> Process record data if necessary
    !>        To avoid duplicate work, the record data is processed
    !>        only once before the iterations start if no source encoding
    !>        For source-encoded FWI, the encoded record data is processed
    !>        at each iteration
    !
    subroutine process_record

        character(len=1024) :: dir_raw, dir_processed
        integer :: i, idata
        !        logical :: enc

        if (record_processing(1) /= '') then

            if (iter == 0) then

                dir_raw = tidy(dir_record)
                dir_processed = tidy(dir_working)//'/record_processed'

                call make_directory(dir_processed)

                do i = shot_in_group(groupid, 1), shot_in_group(groupid, 2)
                    do idata = 1, ndata
                        call process_single_component(dir_raw, dir_processed, &
                            i, data_name(idata), record_processing)
                    end do
                end do


                dir_record = dir_processed

            end if

        end if

        call mpibarrier

    end subroutine

    !
    !> Process synthetic data
    !
    subroutine process_synthetic(ishot)

        integer, intent(in) :: ishot
        integer :: idata

        synthetic_processed = .false.

        if (synthetic_processing(1) /= '') then

            call make_directory(dir_synthetic_processed)

            do idata = 1, ndata
                call process_single_component(dir_synthetic, dir_synthetic_processed, &
                    ishot, data_name(idata), synthetic_processing)
            end do

            ! Turn on synthetic processed flag
            synthetic_processed = .true.

        end if

    end subroutine

    !
    !> Process synthetic data
    !
    subroutine process_adjoint_source(ishot)

        integer, intent(in) :: ishot
        integer :: idata

        if (adjoint_source_processing(1) /= '') then

            do idata = 1, ndata
                call process_single_component(dir_adjoint, dir_adjoint, &
                    ishot, data_name(idata), adjoint_source_processing)
            end do

        end if

    end subroutine

end module
