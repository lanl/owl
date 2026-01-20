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


module mod_dtw

    use libflit

    implicit none

    ! For simplicity, r_inst and r_cuml are both defined to be L2
#define r_cuml(x) ((x)**2)
#define r_inst(x) ((x)**2)

    type gdtw

        ! search space size
        integer :: n
        integer :: n_defualt = 300
        integer :: m = 300
        integer :: m_max = 300
        double precision :: eta = 0.15d0

        ! params and loss, regularizer functionals
        double precision :: cost
        integer, allocatable, dimension(:) :: path
        double precision, allocatable, dimension(:) :: tau
        double precision :: lambda_cuml = 1.0d0
        double precision :: lambda_inst = 0.1d0
        character(len=12) :: loss = 'l2'

        ! slope constraints and boundary conditions
        double precision :: s_min = 1.0d-8
        double precision :: s_max = 1.0d2
        double precision :: s_beta = 0.0d0
        logical :: bc_start_stop = .true.

        ! termination conditions
        integer :: max_iters = 10
        double precision :: epsilon_abs = 1.0d-1
        double precision :: epsilon_rel = 1.0d-2

        logical :: verbose = .true.

    contains
        procedure :: solve => solve_gdtw

    end type

    public :: gdtw
    public :: dtwn

contains

    !
    !> Generalized dynamic time warping
    !
    subroutine solve_gdtw(this, x, y)

        class(gdtw), intent(inout) :: this
        double precision, dimension(:, :), intent(in) :: x, y

        integer :: n, m, iter, i, j, k
        double precision, allocatable, dimension(:) :: t, u, l, u_orig, l_orig, tau_range, yy
        double precision, allocatable, dimension(:, :) :: a, b, tau_graph
        double precision :: dt, slope
        double precision, allocatable, dimension(:, :) :: dist, nu, f
        integer :: j_center, j_opt
        integer, allocatable, dimension(:, :) :: p
        double precision :: path_cost, e_ijk, min_cost
        double precision :: cost_prev
        double precision :: delta = 1.0d-10
        real :: power

        n = this%n
        if (mod(this%m, 2) == 0) then
            m = this%m + 1
        else
            m = this%m
        end if

        this%tau = zeros(n)
        this%path = zeros(n)

        t = linspace(0.0d0, 1.0d0, n)
        dt = 1.0d0/(n - 1.0d0)

        cost_prev = float_tiny

        do iter = 1, this%max_iters

            dist = zeros(n, m)
            nu = zeros(n, m)
            f = zeros(n, m)
            p = ones(n, m)

            ! =================================================================
            ! Compute taus

            if (iter == 1) then

                u = minval(reshape([this%s_beta + this%s_max*t, &
                    this%s_beta + 1 - this%s_min*(1 - t), &
                    ones(n)*1.0d0], [n, 3]), dim=2)
                l = maxval(reshape([this%s_min*t, &
                    -this%s_beta + 1 - this%s_max*(1 - t), &
                    zeros(n)*0.0d0], [n, 3]), dim=2)

                u_orig = u
                l_orig = l

            else

                tau_range = this%eta*(u - l)/2.0

                u = minval(reshape([this%tau + tau_range, u_orig], [n, 2]), dim=2)
                l = maxval(reshape([this%tau - tau_range, l_orig], [n, 2]), dim=2)

            end if

            a = reshape([l, u - l], [n, 2])
            b = transpose(reshape([ones(m)*1.0d0, linspace(0.0d0, 1.0d0, m)], [m, 2]))
            tau_graph = matmul(a, b)

            ! =================================================================
            ! Compute distance matrix

            power = extract_float(this%loss(2:))
            do k = 1, size(x, 2)
                do j = 1, m
                    yy = ginterp(t, x(:, k), tau_graph(:, j), 'linear') - y(:, k)
                    dist(:, j) = dist(:, j) + abs(yy)**power
                end do
            end do

            ! =================================================================
            ! Solve the GDTW

            ! Initialize node costs and f matrix
            do i = 1, n
                do j = 1, m
                    nu(i, j) = dist(i, j) + this%lambda_cuml*r_cuml(tau_graph(i, j) - t(i))
                    f(i, j) = huge(1.0d0)
                end do
            end do

            ! Initialize first row of f matrix
            j_center = (m + 1)/2
            if (this%bc_start_stop) then
                f(1, :) = huge(1.0d0)
                f(1, j_center) = nu(1, j_center)
            else
                f(1, :) = nu(1, :)
            end if

            ! Dynamic programming to fill f, n, and p matrices
            do i = 1, n - 1
                do j = 1, m
                    do k = 1, m

                        slope = (tau_graph(i + 1, k) - tau_graph(i, j))/dt
                        if (slope < this%s_min - delta .or. slope > this%s_max + delta) then
                            cycle
                        end if

                        e_ijk = this%lambda_inst*r_inst(slope)
                        path_cost = f(i, j) + dt*(e_ijk + nu(i + 1, k))
                        if (path_cost < f(i + 1, k)) then
                            f(i + 1, k) = path_cost
                            p(i + 1, k) = j
                        end if

                    end do
                end do
            end do

            ! Determine optimal path
            if (this%bc_start_stop) then
                j_opt = j_center
            else
                min_cost = huge(1.0d0)
                do j = 1, m
                    if (f(n, j) < min_cost) then
                        min_cost = f(n, j)
                        j_opt = j
                    end if
                end do
            end if

            ! Net cost
            this%cost = f(n, j_opt)

            ! Traceback to find optimal tau and path
            do i = n, 1, -1
                this%tau(i) = tau_graph(i, j_opt)
                this%path(i) = j_opt
                j_opt = p(i, j_opt)
            end do

            ! =================================================================
            ! Check stop criterion
            if (this%verbose) then
                call warn(date_time_compact()//' >> GDTW iteration = '//num2str(iter)// &
                    ', cost = '//num2str(this%cost, '(es12.5)'))
            end if

            if (abs(this%cost - cost_prev) <= this%epsilon_abs + this%epsilon_rel*cost_prev) then
                exit
            end if

            cost_prev = this%cost

        end do

        this%tau = (this%tau - linspace(0.0, 1.0, n))*(n - 1)

    end subroutine

    !
    !> Conventional dynamic time warping
    !
    function dtwn(source, target, maxlag, order) result(u)

        real, dimension(:, :) :: source, target
        integer :: maxlag
        integer, optional :: order
        integer, allocatable, dimension(:) :: u

        integer :: n, i, j, k, m
        real, allocatable, dimension(:, :) :: e, d
        integer, allocatable, dimension(:) :: js

        if (present(order)) then
            m = order
        else
            m = 2
        end if

        n = size(target, 1)

        e = zeros(n, n) + float_huge/2.0
        do j = 1, n
            do i = 1, n

                if (abs(j - i) <= maxlag) then
                    ! Shift target to match source
                    e(i, j) = sum((target(i, :) - source(j, :))**2)
                    ! Shift source to match target
                    ! e(i, j) = sum((source(i, :) - target(j, :))**2)
                end if

            end do
        end do

        d = zeros(n, n)
        do j = 1, n
            d(1, j) = e(1, j)
        end do
        do k = 1, m
            do i = 2, n
                d(i, k) = e(i, k) + minval(d(i - 1, 1:k))
            end do
        end do
        do j = m + 1, n
            do i = 3, n
                d(i, j) = e(i, j) + minval(d(i - 1, j - m:j))
            end do
        end do

        u = zeros(n)
        u(n) = n
        do i = n - 1, 1, -1
            js = [u(i + 1) - m:u(i + 1)]
            js = pack(js, mask=(js >= 1))
            u(i) = js(minloc(d(i, js), dim=1))
        end do

        u = u - regspace(0.0, 1.0, n - 1.0)

    end function

end module

module inversion_adjoint_source

    use libflit
    use mod_parameters
    use mod_model
    use mod_source_receiver
    use mod_su
    use mod_data_processing
    use mod_dtw

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    integer, allocatable, dimension(:, :) :: trace_in_group_rank
    integer :: adaptive_hw = 3

contains

    !
    !> Time-shift-based penalty function for AWI and LAWI
    !
    function xcorr_penalty(dt, maxlag) result(w)

        real, intent(in) :: dt
        integer, intent(in) :: maxlag
        real, allocatable, dimension(:) :: w

        integer :: i
        real :: scaling

        scaling = 1.0/sqrt(2.0*log(100.0))

        call alloc_array(w, [-maxlag, maxlag])

        select case (penalty_method)
            case ('linear')
                do i = -maxlag, maxlag
                    w(i) = abs(i*dt/tlag_max)
                end do
            case ('power')
                do i = -maxlag, maxlag
                    w(i) = abs(i*dt/tlag_max)**abs(penalty_power)
                end do
            case ('gaussian')
                do i = -maxlag, maxlag
                    w(i) = exp(-(i*dt)**2/(scaling*tlag_max**2))
                end do
                w = rescale(w, [0.0, 1.0])
                w = w**abs(penalty_power)
                w = 1.0 - w
            case ('exp')
                do i = -maxlag, maxlag
                    w(i) = exp(-abs(i*dt)/(scaling*tlag_max))
                end do
                w = rescale(w, [0.0, 1.0])
                w = w**abs(penalty_power)
                w = 1.0 - w
            case default
                do i = -maxlag, maxlag
                    w(i) = abs(i*dt/tlag_max)
                end do
        end select

    end function xcorr_penalty

    !
    !> L2-norm waveform misfit; the original form of FWI, developed by
    !>    Tarantola, 1984:
    !>    Inversion of seismic reflection data in the acoustic approximation
    !>    Geophysics, doi: 10.1190/1.1441754
    !
    subroutine compute_adjsrc_waveform(seis_obs, seis_syn, seis_adj, val_misfit)

        real, dimension(:), intent(in) :: seis_obs, seis_syn
        real, dimension(:), intent(out) :: seis_adj
        real, intent(out) :: val_misfit

        val_misfit = sum((seis_obs - seis_syn)**2)

        if (yn_save_adjsrc) then
            seis_adj = seis_syn - seis_obs
        end if

    end subroutine compute_adjsrc_waveform

    !
    !> Zero-lag correlation misfit; the method is modified from
    !>    Zhang et al., 2015:
    !>    A stable and practical implementation of least-squares reverse time migration
    !>    Geophysics, doi: 10.1190/GEO2013-0461.1
    !
    subroutine compute_adjsrc_corr(seis_obs, seis_syn, seis_adj, val_misfit)

        real, dimension(:), intent(in) :: seis_obs, seis_syn
        real, dimension(:), intent(inout) :: seis_adj
        real, intent(out) :: val_misfit

        if (maxval(abs(seis_syn)) == 0 .and. maxval(abs(seis_obs)) /= 0) then
            val_misfit = sum((seis_obs/norm2(seis_obs))**2)
            if (yn_save_adjsrc) then
                seis_adj = -seis_obs/norm2(seis_obs)
            end if
            return
        end if

        val_misfit = sum((seis_syn/norm2(seis_syn) - seis_obs/norm2(seis_obs))**2)

        if (yn_save_adjsrc) then
            seis_adj = (dot_product(seis_syn, seis_obs)/norm2(seis_syn)**2*seis_syn - seis_obs) &
                /(norm2(seis_obs)*norm2(seis_syn))
        end if

    end subroutine compute_adjsrc_corr

    !
    !> Envelope misfit; the method is developed by
    !>    Chi et al., 2014:
    !>    Full waveform inversion method using envelope objective function without low frequency data
    !>    Journal of Applied Geophysics, doi: 10.1016/j.jappgeo.2014.07.010
    !>
    !>    Wu et al., 2014:
    !>    Seismic envelope inversion and modulation signal model
    !>    Geophysics, doi: 10.1190/GEO2013-0294.1
    !
    subroutine compute_adjsrc_envelope(seis_obs, seis_syn, seis_adj, val_misfit)

        real, dimension(:), intent(in) :: seis_obs, seis_syn
        real, dimension(:), intent(inout) :: seis_adj
        real, intent(inout) :: val_misfit

        real, allocatable, dimension(:) :: evlp1, evlp2, tempseis, tempcoef

        evlp1 = envelope(seis_obs)
        evlp2 = envelope(seis_syn)
        tempseis = zeros(nt)
        tempcoef = zeros(nt)

        val_misfit = sum((evlp1 - evlp2)**2)

        if (yn_save_adjsrc) then

            if (maxval(evlp2) == 0) then
                tempcoef = 0.0
            else
                tempcoef = (evlp2 - evlp1)/(evlp2 + 1.0e-2*maxval(evlp2))
            end if

            seis_adj = seis_syn*tempcoef - hilbert(hilbert(seis_syn)*tempcoef)

        end if

    end subroutine

    !
    !> Phase misfit; the method is developed by
    !>   Ebru Bozdağ, Jeannot Trampert, Jeroen Tromp, 2011:
    !>   Misfit functions for full waveform inversion based on instantaneous phase and envelope measurements
    !>   Geophysical Journal International, doi: 10.1111/j.1365-246X.2011.04970.x
    !
    subroutine compute_adjsrc_phase(seis_obs, seis_syn, seis_adj, val_misfit)

        real, dimension(:), intent(in) :: seis_obs, seis_syn
        real, dimension(:), intent(inout) :: seis_adj
        real, intent(inout) :: val_misfit

        real, allocatable, dimension(:) :: phiu, phid, env

        phiu = instant_phase(seis_syn)
        phid = instant_phase(seis_obs)

        val_misfit = sum((phid - phiu)**2)

        if (yn_save_adjsrc) then

            env = envelope(seis_syn)
            env = env**2
            env = env + 5.0e-2*maxval(env)

            seis_adj = -(phiu - phid)*hilbert(seis_syn)/env + hilbert((phiu - phid)*seis_syn/env)

        end if

    end subroutine compute_adjsrc_phase

    !
    !> Weighted normalized deconvolution misfit; the method is developed by
    !>   Warner and Guasch, 2016:
    !>   Adaptive waveform inversion: Theory
    !>   Geophysics, doi: 10.1190/geo2015-0387.1
    !
    subroutine compute_adjsrc_adaptive(seis_obs, seis_syn, dt, seis_adj, val_misfit)

        real, dimension(:), intent(in) :: seis_obs, seis_syn
        real, intent(in) :: dt
        real, dimension(:), intent(out) :: seis_adj
        real, intent(out) :: val_misfit

        real, allocatable, dimension(:) :: u, d, xc, p, w, a
        complex, allocatable, dimension(:) :: uu, dd
        integer :: maxlag, n, nnt, tl
        real :: eps

        n = size(seis_obs)
        maxlag = ceiling(tlag_max/dt)
        nnt = next_power_235(2*n)
        tl = nint(0.025*n)

        u = pad(taper(seis_syn, len=[tl, tl]), [0, nnt - n])
        d = pad(taper(seis_obs, len=[tl, tl]), [0, nnt - n])

        uu = fft(u)
        dd = fft(d)

        eps = deconv_eps*maxval(abs(dd*conjg(dd)))

        ! The mathcing filter
        xc = ifft(uu*conjg(dd)/(dd*conjg(dd) + eps), real=.true.)

        ! Penalty function
        call alloc_array(p, [-maxlag, maxlag])
        p = xcorr_penalty(dt, maxlag)
        w = zeros(nnt)
        w(1:maxlag + 1) = p(0:maxlag)
        w(nnt - maxlag + 1:nnt) = p(-maxlag:-1)

        ! Misfit
        val_misfit = 0.5*sum((w*xc)**2)/sum(xc**2)

        ! Adjoint source
        if (yn_save_adjsrc) then

            a = ifft(fft((w**2 - 2.0*val_misfit)*xc/sum(xc**2))*dd/(dd*conjg(dd)+ eps), real=.true.)
            seis_adj = a(1:n)
            seis_adj = taper(seis_adj, [0, tl])

        end if

    end subroutine

    !
    !> Localized weighted normalized deconvolution misfit; the method is developed by
    !>   Yong et al., 2023:
    !>   Localized adaptive waveform inversion: theory and numerical verification
    !>   Geophysical Journal International, doi: 10.1093/gji/ggac496
    !
    subroutine compute_adjsrc_local_adaptive(seis_obs, seis_syn, dt, seis_adj, val_misfit)

        real, dimension(:), intent(in) :: seis_obs, seis_syn
        real, intent(in) :: dt
        real, dimension(:), intent(out) :: seis_adj
        real, intent(out) :: val_misfit

        real, allocatable, dimension(:) :: u, d, xc, p, w, a, g, s
        complex, allocatable, dimension(:) :: uu, dd, u0, d0, af
        integer :: maxlag, n, nnt, tl
        real :: eps, eta, vt
        integer :: i

        n = size(seis_obs)
        maxlag = ceiling(tlag_max/dt)
        nnt = next_power_235(2*n)
        tl = nint(0.025*n)

        u = pad(taper(seis_syn, len=[tl, tl]), [0, nnt - n])
        d = pad(taper(seis_obs, len=[tl, tl]), [0, nnt - n])

        a = regspace(1.0, 1.0, nnt*1.0)
        af = zeros(nnt)
        seis_adj = zeros(nnt)
        s = zeros(nnt)

        val_misfit = 0.0

        ! Compute a common small value for samples to avoid division by zero
        ! According to Yong et al. (2023), less noisy signals need smaller eps
        u0 = fft(u)
        d0 = fft(d)
        eps = deconv_eps*maxval(abs(d0*conjg(d0)))
        if (eps == 0) then
            eps = float_small
        end if

        xc = ifft(u0*conjg(d0)/(d0*conjg(d0) + eps), real=.true.)
        eta = deconv_eps*maxval(xc**2)
        if (eta == 0) then
            eta = float_small
        end if

        ! Weights for the matching filter
        call alloc_array(p, [-maxlag, maxlag])
        p = xcorr_penalty(dt, maxlag)
        w = zeros(nnt)
        w(1:maxlag + 1) = p(0:maxlag)
        w(nnt - maxlag + 1:nnt) = p(-maxlag:-1)

        ! Loop through all time samples
        do i = 1, n

            ! Gaussian kernel of the Gabor transform
            g = gaussian(a, i*1.0, lawi_sigma/dt)

            uu = fft(u*g)
            dd = fft(d*g)

            ! The local mathcing filter
            xc = ifft(uu*conjg(dd)/(dd*conjg(dd) + eps), real=.true.)

            ! Misfit
            vt = sum((w*xc)**2)/(sum(xc**2) + eta)
            val_misfit = val_misfit + vt**2

            ! Adjoint source
            if (yn_save_adjsrc) then
                af = fft(vt*(w**2 - vt)*xc/(sum(xc**2) + eta))*dd/(dd*conjg(dd)+ eps)
                seis_adj = seis_adj + ifft(af, real=.true.)*g
                s = s + g**2
            end if

        end do

        if (yn_save_adjsrc) then
            seis_adj = seis_adj(1:n)/s(1:n)
            seis_adj = taper(seis_adj, [0, tl])
        end if

    end subroutine

    !
    !> GDTW-based time-variant phase shift misfit for vector data
    !
    subroutine compute_adjsrc_adaptive_spacetime(seis_obs, seis_syn, dt, seis_adj, val_misfit)

        real, dimension(:, :), intent(in) :: seis_obs, seis_syn
        real, intent(in) :: dt
        real, dimension(:), intent(out) :: seis_adj
        real, intent(out) :: val_misfit

        real, allocatable, dimension(:, :) :: u, d
        complex, allocatable, dimension(:, :) :: uu, dd
        real, allocatable, dimension(:) :: xc, p, w, a
        complex, allocatable, dimension(:) :: uf
        real, allocatable, dimension(:) :: df
        integer :: maxlag, n, nnt, tl, nc, i, l
        real :: eps

        n = size(seis_obs, 1)
        nc = size(seis_obs, 2)
        maxlag = ceiling(tlag_max/dt)
        nnt = next_power_235(4*n)
        tl = nint(0.025*n)

        u = pad(taper(seis_syn, len=[tl, tl, 0, 0]), [0, nnt - n, 0, 0])
        d = pad(taper(seis_obs, len=[tl, tl, 0, 0]), [0, nnt - n, 0, 0])

        uu = fft(u, along=1)
        dd = fft(d, along=1)

        eps = 0.0
        uf = zeros(nnt)
        df = zeros(nnt)
        xc = zeros(nnt)
        l = 0
        do i = 1, nc
            uf = uu(:, i)*conjg(dd(:, i))
            df = dd(:, i)*conjg(dd(:, i))
            eps = deconv_eps*maxval(df)
            if (eps > 0) then
                xc = xc + ifft(uf/(df + eps), real=.true.)
                l = l + 1
            end if
        end do
        if (l > 0) then
            xc = xc/l
        end if

        ! Penalty function
        call alloc_array(p, [-maxlag, maxlag])
        p = xcorr_penalty(dt, maxlag)
        w = zeros(nnt)
        w(1:maxlag + 1) = p(0:maxlag)
        w(nnt - maxlag + 1:nnt) = p(-maxlag:-1)

        ! Misfit
        if (l > 0) then
            val_misfit = 0.5*sum((w*xc)**2)/sum(xc**2)
        else
            val_misfit = 0
        end if

        ! Adjoint source
        if (yn_save_adjsrc) then

            i = (nc + 1)/2
            df = dd(:, i)*conjg(dd(:, i))
            eps = deconv_eps*maxval(df)
            a = ifft(fft((w**2 - 2*val_misfit)*xc/sum(xc**2))*dd(:, i)/(df + eps), real=.true.)
            seis_adj = a(1:n)
            seis_adj = taper(seis_adj, [0, tl])

        end if

    end subroutine

    !
    !> Generalized GDTW-based time-variant phase shift misfit for scalar data; in reference to
    !>
    !>   Ma and Hale, 2013:
    !>   Wave-equation reflection traveltime inversion with dynamic warping and full-waveform inversion
    !>   Geophysics, doi: 10.1190/GEO2013-0004.1
    !>
    !>   DTW can be viewed as a case of time-variant phase shift based methods, including
    !>   DTW, soft DTW, optimal transport, graph-space optimal transport, etc.
    !>   The general idea of these method is estimating a time shift field to move
    !>   d_obs to d_syn, that is, d_obs(t + tau) ~ d_syn(t).
    !
    subroutine compute_adjsrc_dtw(seis_obs, seis_syn, d, weight, seis_adj, val_misfit, tau_map)

        real, dimension(:, :, :), intent(in) :: seis_obs, seis_syn
        real, intent(in) :: d
        real, dimension(:), intent(in) :: weight
        real, dimension(:, :, :), intent(out) :: seis_adj
        real, dimension(:), intent(out) :: val_misfit
        real, dimension(:, :), intent(out), optional :: tau_map

        integer :: n, nr, nc, i, nn, ic
        real :: dd
        double precision, allocatable, dimension(:, :, :) :: dobs, dsyn
        type(gdtw), allocatable, dimension(:):: g
        real, allocatable, dimension(:, :) :: tau
        real, allocatable, dimension(:) :: t, warped

        integer :: j, nw, nv, l
        real, allocatable, dimension(:, :) :: q
        real, allocatable, dimension(:) :: u, v
        real :: eta
        real :: m_phase, m_amp

        seis_adj = 0.0
        val_misfit = 0.0

        n = size(seis_obs, 1)
        nr = size(seis_obs, 2)
        nc = size(seis_obs, 3)

        if (adj_nt > 0) then
            nn = adj_nt
            dd = (n - 1.0)*d/(nn - 1.0)
        else
            nn = n
            dd = d
        end if

        ! Resample data
        dsyn = zeros(nn, nr, nc)
        dobs = zeros(nn, nr, nc)

        !$omp parallel do private(i, j) schedule(auto)
        do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

            if (weight(i) == 0) then
                cycle
            end if

            ! Normalize data and resample
            ! The normalize is based on all components, if the input is a vector data
            do j = 1, nc
                if (nn /= n) then
                    dsyn(:, i, j) = interp_to(seis_syn(:, i, j)/norm2(seis_syn(:, i, :)), nn, method='cubic')
                    dobs(:, i, j) = interp_to(seis_obs(:, i, j)/norm2(seis_obs(:, i, :)), nn, method='cubic')
                else
                    dsyn(:, i, j) = seis_syn(:, i, j)/norm2(seis_syn(:, i, :))
                    dobs(:, i, j) = seis_obs(:, i, j)/norm2(seis_obs(:, i, :))
                end if
            end do

        end do
        !$omp end parallel do

        call allreduce_array_group(dsyn)
        call allreduce_array_group(dobs)
        call mpibarrier_group

        ! GDTW
        allocate(g(1:nr))
        tau = zeros(nn, nr)

        nw = dtw_trc
        call pad_array(dsyn, [0, 0, nw, nw, 0, 0])
        call pad_array(dobs, [0, 0, nw, nw, 0, 0])

        !$omp parallel do private(i) schedule(auto)
        do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

            if (weight(i) == 0) then
                cycle
            end if

            ! Compute time shifts
            g(i)%n = nn
            g(i)%m = nint(2*tlag_max/dd) + 1
            g(i)%max_iters = dtw_niter
            g(i)%lambda_inst = dtw_rinst
            g(i)%lambda_cuml = dtw_rcuml
            g(i)%epsilon_abs = dtw_epsabs
            g(i)%epsilon_rel = dtw_epsrel
            g(i)%loss = dtw_loss
            g(i)%verbose = .false.

            call g(i)%solve(reshape(dobs(1:nn, i - nw:i + nw, :), [nn, nc*(2*nw + 1)]), &
                reshape(dsyn(1:nn, i - nw:i + nw, :), [nn, nc*(2*nw + 1)]))
            tau(:, i) = g(i)%tau*dd

            deallocate(g(i)%tau, g(i)%path)

        end do
        !$omp end parallel do

        call allreduce_array_group(tau)
        call mpibarrier_group

        ! Smoothing the phase shifts
        if (dtw_smooth_median > 0 .and. nr > 3) then
            tau = median_filt(tau, [0, nint(dtw_smooth_median)])
        end if

        if (dtw_smooth_gaussian > 0 .and. nr > max(3.0, dtw_smooth_gaussian)) then
            tau = gauss_filt(tau, [0.0, dtw_smooth_gaussian])
        end if

        !$omp parallel do private(i) schedule(auto)
        do i = 1, nr
            if (weight(i) == 0) then
                tau(:, i) = 0
            end if
        end do
        !$omp end parallel do

        ! Warping d_obs
        select case (dtw_form)

            case default

                q = zeros(nn, nr)
                nw = nint(tlag_max/dd)
                nv = 2*nw + 1
                u = zeros(nv)
                v = zeros(nv)

                t = regspace(0.0, 1.0, nn - 1.0)*dd
                !$omp parallel do private(i, j, l, ic, u, v, warped, eta) schedule(auto)
                do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

                    if (weight(i) == 0) then
                        cycle
                    end if

                    do ic = 1, nc

                        dobs(:, i, ic) = return_normal(ginterp(t, real(dobs(:, i, ic)), tau(:, i) + t, 'linear'))
                        eta = max(norm2(dsyn(:, i, ic)), norm2(dobs(:, i, ic)))*1.0e-2
                        if (eta == 0) then
                            eta = 1.0e-9
                        end if

                        do j = 1, nn
                            u = 0
                            v = 0
                            l = min(nn, j + nw) - max(1, j - nw) + 1
                            u(1:l) = dsyn(max(1, j - nw):min(nn, j + nw), i, ic)
                            v(1:l) = dobs(max(1, j - nw):min(nn, j + nw), i, ic)
                            q(j, i) = q(j, i) + sum(u*v)/(norm2(u)*norm2(v) + eta)
                        end do

                    end do

                    q(:, i) = max(0.0, abs(q(:, i))**0.2)
                    q(:, i) = median_filt(q(:, i), 3)
                    q(:, i) = gauss_filt(q(:, i), 3.0)
                    q(:, i) = q(:, i)/norm2(q(:, i))
                    q(:, i) = return_normal(rescale(q(:, i), [0.0, 1.0]))

                end do
                !$omp end parallel do

                call allreduce_array_group(q)

                if (dtw_smooth_median > 0) then
                    q = median_filt(q, [0, nint(dtw_smooth_median)])
                end if

                if (dtw_smooth_gaussian > 0) then
                    q = gauss_filt(q, [0.0, dtw_smooth_gaussian])
                end if

                !$omp parallel do private(i) schedule(auto)
                do i = 1, nr
                    if (maxval(abs(q(:, i))) == 0) then
                        q(:, i) = 1
                    end if
                end do
                !$omp end parallel do

            case ('amp')
                q = ones(nn, nr)

        end select

        ! Interpolate computed time shifts to output dimensions
        if (nn /= n) then
            tau = interp_to(tau, [n, nr], method=['pchip', 'nearest'])
            q = interp_to(q, [n, nr], method=['pchip', 'nearest'])
        end if

        ! Compute adjoint source
        warped = zeros(n)
        val_misfit = 0
        t = regspace(0.0, 1.0, n - 1.0)*d

        if (yn_save_adjsrc) then
            if (present(tau_map)) then
                tau_map = tau
            end if
        end if

        !$omp parallel do private(i, warped, ic, m_phase, m_amp) schedule(auto)
        do i = 1, nr

            if (weight(i) == 0) then
                cycle
            end if

            select case (dtw_form)

                case ('phase')

                    val_misfit(i) = sum(q(:, i)*tau(:, i)**2)

                    if (yn_save_adjsrc) then
                        do ic = 1, nc
                            warped = return_normal(ginterp(t, seis_obs(:, i, ic), tau(:, i) + t, 'cubic'))
                            seis_adj(:, i, ic) = q(:, i)*tau(:, i)*deriv(warped)/d
                        end do
                    end if

                case ('amp')

                    do ic = 1, nc

                        warped = return_normal(ginterp(t, seis_obs(:, i, ic), tau(:, i) + t, 'cubic'))
                        warped = seis_syn(:, i, ic) - warped
                        val_misfit(i) = val_misfit(i) + sum(warped**2)

                        if (yn_save_adjsrc) then
                            seis_adj(:, i, ic) = warped
                        end if

                    end do

                case ('phase+amp')

                    m_phase = sum(q(:, i)*tau(:, i)**2)

                    do ic = 1, nc

                        val_misfit(i) = val_misfit(i) + m_phase*(1.0 + dtw_amp_weight)

                        if (yn_save_adjsrc) then
                            warped = return_normal(ginterp(t, seis_obs(:, i, ic), tau(:, i) + t, 'cubic'))
                            m_amp = sum((seis_syn(:, i, ic) - warped)**2)
                            seis_adj(:, i, ic) = q(:, i)*tau(:, i)*deriv(warped)/d &
                                + dtw_amp_weight*m_phase/(m_amp + float_tiny)*(seis_syn(:, i, ic) - warped)
                        end if

                    end do

            end select

        end do
        !$omp end parallel do

    end subroutine

    !
    !> Compute adjoint source
    !
    subroutine compute_adjoint_source(srcindex, misfit_sum, dir_syn, dir_obs, dir_adj)

        integer, intent(in) :: srcindex
        real, intent(inout) :: misfit_sum
        character(len=*), intent(in) :: dir_syn, dir_obs, dir_adj

        integer :: i, nr, nn, nc, ic
        integer :: nt_base
        real :: dt_base
        character(len=1024) :: filename, file_syn, file_adj
        type(su) :: seismo_obs, seismo_syn
        real, allocatable, dimension(:) :: err_trc, dobs, dsyn, dadj
        real, allocatable, dimension(:, :) :: weight
        real :: dd, amp_obs, amp_syn, m_obs, m_syn
        real, allocatable, dimension(:) :: norm_obs, norm_syn
        real, allocatable, dimension(:, :, :) :: data_obs, data_syn, data_adj
        real, allocatable, dimension(:, :) :: dsyn2, dobs2
        real, allocatable, dimension(:, :, :) :: tau

        misfit_sum = 0.0

        ! Get dimensions
        nc = size(data_name)
        nr = gmtr(set_gmtrid(srcindex))%nr

        ! Divide the traces into different group workers
        call alloc_array(trace_in_group_rank, [0, nrank_group - 1, 1, 2])
        call cut(1, nr, nrank_group, trace_in_group_rank)

        ! Set base nt, dt
        filename = '/shot_'//num2str(set_srcid(srcindex))//'_'//'seismogram_'//tidy(data_name(1))//'.su'
        file_syn = tidy(dir_syn)//tidy(filename)
        call seismo_syn%load(file_syn, nr=1)
        nt_base = seismo_syn%nt
        dt_base = seismo_syn%dt

        ! Get data
        data_obs = zeros(nt_base, nr, nc)
        data_syn = zeros(nt_base, nr, nc)

        weight = zeros(nr, nc)
        do ic = 1, nc

            filename = '/shot_'//num2str(set_srcid(srcindex))//'_'//'seismogram_'//tidy(data_name(ic))//'.su'

            call seismo_syn%load(tidy(dir_syn)//tidy(filename), nr=nr)
            call seismo_obs%load(tidy(dir_obs)//tidy(filename), nr=nr)

            norm_obs = zeros(nr)
            norm_syn = zeros(nr)
            !$omp parallel do private(i)
            do i = 1, nr
                norm_obs(i) = norm2(seismo_obs%trace(i)%data)
                norm_syn(i) = norm2(seismo_syn%trace(i)%data)
            end do
            !$omp end parallel do

            m_obs = maxval(norm_obs)
            m_syn = maxval(norm_syn)

            !$omp parallel do private(i, amp_obs, amp_syn) schedule(auto)
            do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

                amp_obs = norm2(seismo_obs%trace(i)%data)
                amp_syn = norm2(seismo_syn%trace(i)%data)

                if (.not.(gmtr(set_gmtrid(srcindex))%recr(i)%weight == 0 &
                        .or. amp_obs < trace_discard_threshold*m_obs &
                        .or. amp_syn < trace_discard_threshold*m_syn &
                        .or. amp_obs*amp_syn == 0)) then

                    if (seismo_syn%nt /= nt_base .or. seismo_syn%dt /= dt_base) then
                        data_syn(:, i, ic) = interp(seismo_syn%trace(i)%data, seismo_syn%nt, seismo_syn%dt, 0.0, &
                            nt_base, dt_base, 0.0, 'cubic')
                    else
                        data_syn(:, i, ic) = seismo_syn%trace(i)%data
                    end if
                    if (seismo_obs%nt /= nt_base .or. seismo_obs%dt /= dt_base) then
                        data_obs(:, i, ic) = interp(seismo_obs%trace(i)%data, seismo_obs%nt, seismo_obs%dt, 0.0, &
                            nt_base, dt_base, 0.0, 'cubic')
                    else
                        data_obs(:, i, ic) = seismo_obs%trace(i)%data
                    end if

                    weight(i, ic) = 1.0

                end if

            end do
            !$omp end parallel do

        end do

        call mpibarrier_group
        call allreduce_array_group(weight)
        call allreduce_array_group(data_syn)
        call allreduce_array_group(data_obs)

        if (adj_nt > 0) then
            nn = adj_nt
            dd = (nt_base - 1.0)*dt_base/(nn - 1.0)
        else
            nn = nt_base
            dd = dt_base
        end if

        err_trc = zeros(nr)
        data_adj = zeros(nt_base, nr, nc)
        dobs = zeros(nn)
        dsyn = zeros(nn)
        dadj = zeros(nn)

        tau = zeros_like(data_adj)

        select case (misfit_type)

            case default
                ! Component by component, trace by trace

                do ic = 1, nc

                    err_trc = 0.0

                    !$omp parallel do private(i, dobs, dsyn, dadj) schedule(auto)
                    do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

                        if (weight(i, ic) /= 0) then

                            ! Resample data to the desired length (adj_nt)
                            if (nt_base /= nn) then
                                dsyn = interp_to(data_syn(:, i, ic), nn, method='cubic')
                                dobs = interp_to(data_obs(:, i, ic), nn, method='cubic')
                            else
                                dsyn = data_syn(:, i, ic)
                                dobs = data_obs(:, i, ic)
                            end if

                            ! Compute misfit and adjoint source
                            select case (misfit_type)

                                case ('waveform')
                                    call compute_adjsrc_waveform(dobs, dsyn, dadj, err_trc(i))

                                case ('corr')
                                    call compute_adjsrc_corr(dobs, dsyn, dadj, err_trc(i))

                                case ('adaptive')
                                    call compute_adjsrc_adaptive(dobs, dsyn, dd, dadj, err_trc(i))

                                case ('local-adaptive')
                                    call compute_adjsrc_local_adaptive(dobs, dsyn, dd, dadj, err_trc(i))

                                case ('envelope')
                                    call compute_adjsrc_envelope(dobs, dsyn, dadj, err_trc(i))

                                case ('phase')
                                    call compute_adjsrc_phase(dobs, dsyn, dadj, err_trc(i))

                                case default
                                    call compute_adjsrc_waveform(dobs, dsyn, dadj, err_trc(i))

                            end select

                            ! Resample adjoint to
                            if (nt_base /= nn) then
                                data_adj(:, i, ic) = interp_to(dadj, nt_base, method='cubic')
                            else
                                data_adj(:, i, ic) = dadj
                            end if

                        end if

                    end do
                    !$omp end parallel do

                    call mpibarrier_group

                    call allreduce_array_group(err_trc)
                    if (yn_save_adjsrc) then
                        call allreduce_array_group(data_adj(:, :, ic))
                    end if

                    misfit_sum = misfit_sum + sum(err_trc)

                end do

            case ('adaptive-spacetime')
                ! Group nearby traces for computing the matching filter

                call readpar_xint(file_parameter, 'adaptive_half_window', adaptive_hw, 3, iter*1.0)

                do ic = 1, nc

                    err_trc = 0.0

                    if (nt_base /= nn) then
                        dsyn2 = interp_to(data_syn(:, :, ic), [nn, nr], method=['cubic', 'nearest'])
                        dobs2 = interp_to(data_obs(:, :, ic), [nn, nr], method=['cubic', 'nearest'])
                    else
                        dsyn2 = data_syn(:, :, ic)
                        dobs2 = data_obs(:, :, ic)
                    end if

                    call alloc_array(dsyn2, [1, nn, -adaptive_hw + 1, nr + adaptive_hw], &
                        source=pad(dsyn2, [0, 0, adaptive_hw, adaptive_hw]))
                    call alloc_array(dobs2, [1, nn, -adaptive_hw + 1, nr + adaptive_hw], &
                        source=pad(dobs2, [0, 0, adaptive_hw, adaptive_hw]))

                    !$omp parallel do private(i, dadj) schedule(auto)
                    do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

                        if (weight(i, ic) /= 0) then

                            call compute_adjsrc_adaptive_spacetime(dobs2(:, i - adaptive_hw:i + adaptive_hw), &
                                dsyn2(:, i - adaptive_hw:i + adaptive_hw), dd, dadj, err_trc(i))

                            ! Resample adjoint to
                            if (nt_base /= nn) then
                                data_adj(:, i, ic) = interp_to(dadj, nt_base, method='cubic')
                            else
                                data_adj(:, i, ic) = dadj
                            end if

                        end if

                    end do
                    !$omp end parallel do

                    call mpibarrier_group

                    call allreduce_array_group(err_trc)
                    if (yn_save_adjsrc) then
                        call allreduce_array_group(data_adj(:, :, ic))
                    end if

                    misfit_sum = misfit_sum + sum(err_trc)

                end do

            case ('dtw')
                ! For DTW, the misfit is computed shot gather by shot gathter

                do ic = 1, nc
                    call compute_adjsrc_dtw(data_obs(:, :, ic:ic), data_syn(:, :, ic:ic), dt_base, &
                        weight(:, ic), data_adj(:, :, ic:ic), err_trc, tau(:, :, ic))
                    misfit_sum = misfit_sum + sum(err_trc)
                end do

            case ('dtw-vector')
                ! For vector DTW, the misfit is computed shot gather by shot gather, but also all components together

                call compute_adjsrc_dtw(data_obs, data_syn, dt_base, sum(weight, dim=2), data_adj, err_trc, tau(:, :, 1))
                misfit_sum = misfit_sum + sum(err_trc)

        end select

        ! Save adjoint source if necessary
        if (yn_save_adjsrc .and. rankid_group == 0) then

            call make_directory(dir_adj)

            do ic = 1, nc

                filename = '/shot_'//num2str(set_srcid(srcindex))//'_seismogram_'//tidy(data_name(ic))//'.su'
                file_adj = tidy(dir_adj)//tidy(filename)
                call seismo_syn%from_array(data_adj(:, :, ic))
                call seismo_syn%clean()
                call seismo_syn%output(file_adj)

                select case (misfit_type)
                    case ('dtw')
                        filename = '/shot_'//num2str(set_srcid(srcindex))//'_time_shift_'//tidy(data_name(ic))//'.su'
                        file_adj = tidy(dir_adj)//tidy(filename)
                        call seismo_syn%from_array(tau(:, :, ic))
                        call seismo_syn%clean()
                        call seismo_syn%output(file_adj)
                    case ('dtw-vector')
                        if (ic == 1) then
                            filename = '/shot_'//num2str(set_srcid(srcindex))//'_time_shift.su'
                            file_adj = tidy(dir_adj)//tidy(filename)
                            call seismo_syn%from_array(tau(:, :, 1))
                            call seismo_syn%clean()
                            call seismo_syn%output(file_adj)
                        end if
                end select

            end do

        end if

        call mpibarrier_group

    end subroutine

end module
