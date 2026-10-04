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


module mod_source_receiver

    use libflit
    use mod_parameters
    use mod_model

    implicit none

    ! a receiver
    type, public :: receiver

        ! receiver location in x, y and z, unit in m
        real :: x = 0.0, y = 0.0, z = 0.0
        ! absolute offset
        real :: aoff
        ! gain
        real :: weight = 1.0

        ! Nearest integer grid point on the mesh
        integer :: gx, gy, gz
        ! Nearest half-integer grid point on the mesh
        integer :: hx, hy, hz
        ! Kaiser windowed sinc function for interpolation
        real, allocatable, dimension(:) :: interp_ix, interp_iy, interp_iz
        real, allocatable, dimension(:) :: interp_hx, interp_hy, interp_hz

        ! valid
        logical :: valid = .true.

        !    contains
        !        procedure, private, pass(this) :: copy_receiver
        !        generic :: assignment(=) => copy_receiver

    end type receiver

    ! a source
    type, public :: source

        !        ! source ID
        !        integer :: id = 1

        ! source location in x, y and z, unit in meter
        real :: x = 0.0, y = 0.0, z = 0.0

        ! central frequency, unit in Hz
        real :: f0 = 20.0
        ! amplitude
        real :: amp = 1.0
        ! phase delay, unit in sec.
        real :: t0 = 0.0
        ! polar and azimuth angle if force vector, unit deg.
        real :: polar = 0.0
        real :: azimuth = 0.0

        integer :: nt, hnt

        ! Soruce mechanism
        character(len=24) :: mechanism = 'explosion'

        ! Source time function: gaussian, gaussian_deriv, ricker, ricker_deriv, ormsby
        character(len=24) :: wavelet = 'ricker'

        integer :: morlet_cycle
        real, allocatable, dimension(:) :: ormsby_freqs

        character(len=1024) :: file_custom_stf

        real, allocatable, dimension(:) :: stf_t, stf_amp, stf

        ! moment tensor, 3x3 array
        real, dimension(3, 3) :: moment_tensor

        ! Nearest integer position on the mesh
        integer :: gx, gy, gz
        ! Nearest half-integer position on the mesh
        integer :: hx, hy, hz

        ! Kaiser windowed sinc function for interpolation
        real, allocatable, dimension(:) :: interp_ix, interp_iy, interp_iz
        real, allocatable, dimension(:) :: interp_hx, interp_hy, interp_hz

        ! valid
        logical :: valid = .true.

        !    contains
        !        procedure, private, pass(this) :: copy_source
        !        generic :: assignment(=) => copy_source
        integer :: time_integration = 0

    end type source

    type source_receiver_geometry

        integer :: id = 1
        integer :: ns = 1
        type(source), allocatable, dimension(:) :: srcr
        integer :: nr = 1
        type(receiver), allocatable, dimension(:) :: recr

        real :: ox = 0, oy = 0, oz = 0
        real :: dx = 1, dy = 1, dz = 1
        integer :: nx = 1, ny = 1, nz = 1
        real :: xmin = -float_huge, xmax = +float_huge
        real :: ymin = -float_huge, ymax = +float_huge
        real :: zmin = -float_huge, zmax = +float_huge
        real :: sxmin = -float_huge, sxmax = +float_huge
        real :: symin = -float_huge, symax = +float_huge
        real :: szmin = -float_huge, szmax = +float_huge
        real :: rxmin = -float_huge, rxmax = +float_huge
        real :: rymin = -float_huge, rymax = +float_huge
        real :: rzmin = -float_huge, rzmax = +float_huge

        integer :: nt
        real :: dt
        real :: f0factor = 1.0
        real, allocatable, dimension(:) :: filtering_freqs, filtering_coefs
        real, allocatable, dimension(:, :, :) :: z_i, z_h, dz_i, dz_h

    contains
        procedure, public :: prepare_geometry
        procedure, public :: prepare_stf

    end type source_receiver_geometry

    type(source_receiver_geometry), allocatable, dimension(:) :: gmtr

    integer, parameter :: nkw = 4
    double precision, parameter, private :: kaiser_b0 = 6.31

contains

    !
    ! Set shot id name
    !
    function set_srcid(srcindex) result(srcid)

        integer :: srcindex
        integer :: srcid

        !        if (yn_source_encoding) then
        !            srcid = srcindex
        !        else
        srcid = gmtr(srcindex)%id
        !        end if

    end function

    function set_gmtrid(srcindex) result(gmtrid)

        integer :: srcindex
        integer :: gmtrid

        !        if (yn_source_encoding) then
        !            gmtrid = shot_in_super(srcindex, 1)
        !        else
        gmtrid = srcindex
        !        end if

    end function

    !
    !> Load geometry
    !
    subroutine load_geometry

        integer :: i, j, l
        character(len=1024) :: dir_geometry, fs, file_source
        logical, allocatable, dimension(:) :: qs
        integer, allocatable, dimension(:) :: usid
        type(source_receiver_geometry), allocatable, dimension(:) :: g
        character(len=1024) :: ff
        real, allocatable, dimension(:) :: fq
        character(len=:), allocatable, dimension(:) :: sq
        integer :: nq

        logical, allocatable, dimension(:) :: src_below_free_surface
        logical, allocatable, dimension(:) :: rec_below_free_surface
        integer :: n
        real, allocatable, dimension(:, :) :: topo
#ifdef _dim2_
        real, allocatable, dimension(:) :: topox, topo_i
#endif
#ifdef _dim3_
        real, allocatable, dimension(:) :: topox, topoy
        real, allocatable, dimension(:, :) :: topoz, topo_ixiy
#endif
        real, allocatable, dimension(:) :: stp, rtp
        real :: topo_max

        if (rankid == 0) then
            call warn(date_time_compact()//' Loading geometry... ')
        end if

        if (allocated(gmtr)) then
            deallocate (gmtr)
        end if
        allocate (gmtr(1:ns))
        if (.not. file_exists(file_geometry)) then
            call warn(date_time_compact()//' @'//get_hostname()//' Error: Geometry file not found. ')
            stop
        else
            dir_geometry = get_file_directory(file_geometry)
        end if

        ! Read in names of separate geometry files
        open (11, file=tidy(file_geometry), action='read', status='old')
        do i = 1, ns

            read (11, *) file_source

            ! Read source parameters
            fs = tidy(dir_geometry)//tidy(file_source)

            open (3, file=tidy(fs), status='old', action='read')

            read (3, *) gmtr(i)%id
            read (3, *) gmtr(i)%ns
            allocate (gmtr(i)%srcr(1:gmtr(i)%ns))

            do j = 1, gmtr(i)%ns

                ! Source location: x, y, z
                read (3, *) gmtr(i)%srcr(j)%x, gmtr(i)%srcr(j)%y, gmtr(i)%srcr(j)%z

                ! Source mechanism: type of mechanism, polar, azimuth, moment tensor
                read (3, '(a)') ff
                sq = split_string(ff, ' ')
                gmtr(i)%srcr(j)%mechanism = tidy(sq(1))
                call assert(any(gmtr(i)%srcr(j)%mechanism == ['explosion', 'force', 'mt']), &
                    ' <load_geometry> Error: mechanism must be one of explosion, force, mt. ')
                select case (gmtr(i)%srcr(j)%mechanism)
                    case default
                        ! Default is explosion
                    case ('explosion')
                        ! Explosion source
                    case ('force')
                        ! Force vector
                        gmtr(i)%srcr(j)%polar = extract_float(tidy(sq(2)))*const_deg2rad
                        gmtr(i)%srcr(j)%azimuth = extract_float(tidy(sq(3)))*const_deg2rad
                    case ('mt')
                        ! Moment tensor
                        gmtr(i)%srcr(j)%moment_tensor(1, 1) = extract_float(tidy(sq(2)))
                        gmtr(i)%srcr(j)%moment_tensor(2, 2) = extract_float(tidy(sq(3)))
                        gmtr(i)%srcr(j)%moment_tensor(3, 3) = extract_float(tidy(sq(4)))
                        gmtr(i)%srcr(j)%moment_tensor(1, 2) = extract_float(tidy(sq(5)))
                        gmtr(i)%srcr(j)%moment_tensor(1, 3) = extract_float(tidy(sq(6)))
                        gmtr(i)%srcr(j)%moment_tensor(2, 3) = extract_float(tidy(sq(7)))
                end select

                ! Source time function: type of wavelet, f0, max amplitudue, delay time
                read (3, '(a)') ff
                sq = split_string(ff, ' ')
                gmtr(i)%srcr(j)%wavelet = tidy(sq(1))
                gmtr(i)%srcr(j)%f0 = extract_float(sq(2))
                gmtr(i)%srcr(j)%amp = extract_float(sq(3))
                gmtr(i)%srcr(j)%t0 = extract_float(sq(4))

                ! Source time function: if custom wavelet then read the name of the ascii file containing the stf
                if (gmtr(i)%srcr(j)%wavelet == 'custom') then
                    read (3, '(a)') gmtr(i)%srcr(j)%file_custom_stf
                end if

                ! Source time function: frequency filtering of the source time function
                read (3, '(a)') ff
                call extract_nfloat(tidy(ff), ' ', fq)
                nq = size(fq)
                nq = floor(nq/2.0)
                if (nq <= 1) then
                    if (gmtr(i)%srcr(j)%wavelet == 'ormsby') then
                        gmtr(i)%srcr(j)%ormsby_freqs = [0.01, 0.1, 2.9, 3.0]*gmtr(i)%srcr(j)%f0
                    end if
                else
                    if (gmtr(i)%srcr(j)%wavelet == 'ormsby') then
                        gmtr(i)%srcr(j)%ormsby_freqs = fq(1:nq)
                    else
                        if (any(fq(1:nq) /= 0)) then
                            gmtr(i)%filtering_freqs = fq(1:nq)
                            gmtr(i)%filtering_coefs = fq(nq + 1:2*nq)
                        end if
                    end if
                end if

                ! If frequency filtering is specificed through parameters, then do filtering
                ! This is mainly for the convenience of multiscale FWI
                if (src_filt_freqs(1) >= 0) then
                    gmtr(i)%filtering_freqs = src_filt_freqs
                    gmtr(i)%filtering_coefs = src_filt_coefs
                end if

            end do

            ! Read receiver parameters
            read (3, *) gmtr(i)%nr
            allocate (gmtr(i)%recr(1:gmtr(i)%nr))
            do j = 1, gmtr(i)%nr
                read (3, *) gmtr(i)%recr(j)%x, gmtr(i)%recr(j)%y, gmtr(i)%recr(j)%z, gmtr(i)%recr(j)%weight
            end do

            close (3)

            ! For 2D model, set source and receiver y to a proper value
            if (ny == 1) then
                gmtr(i)%srcr(:)%y = 0.5*(ymin + ymax)
                gmtr(i)%recr(:)%y = 0.5*(ymin + ymax)
            end if

            ! f0 factor
            gmtr(i)%f0factor = f0_factor

            ! Check validity of geometry
            if (yn_free_surface .and. file_topo /= '') then

                n = count_nonempty_lines(file_topo)

#ifdef _dim2_

                topo = load(file_topo, n, 2, ascii=.true.)
                topox = regspace(0.0, dx, (nx + 2*pml + 1)*dx) - (pml + 1)*dx
                call alloc_array(topo_i, [-pml, nx + pml + 1], source=ginterp(topo(:, 1), topo(:, 2), topox, method=topo_interp))

                if (measure_source_depth_from_surface) then
                    src_below_free_surface = trues(gmtr(i)%ns)
                else
                    stp = ginterp(topox, -topo_i, gmtr(i)%srcr(:)%x, topo_interp)
                    src_below_free_surface = gmtr(i)%srcr(:)%z >= stp
                end if

                if (measure_receiver_depth_from_surface) then
                    rec_below_free_surface = trues(gmtr(i)%nr)
                else
                    rtp = ginterp(topox, -topo_i, gmtr(i)%recr(:)%x, topo_interp)
                    rec_below_free_surface = gmtr(i)%recr(:)%z >= rtp
                end if

                topo_max = maxval(topo_i)

#endif

#ifdef _dim3_

                topo = load(file_topo, n, 3, ascii=.true.)
                topox = meshgrid([nx + 2*pml + 2, ny + 2*pml + 2]*2 + [1, 1], [dx, dy]/2.0, [-(pml + 1)*dx + ox, -(pml + 1)*dy + oy], dim=1)
                topoy = meshgrid([nx + 2*pml + 2, ny + 2*pml + 2]*2 + [1, 1], [dx, dy]/2.0, [-(pml + 1)*dx + ox, -(pml + 1)*dy + oy], dim=2)
                topoz = reshape(ginterp(topo(:, 1), topo(:, 2), topo(:, 3), topox, topoy), [nx + 2*pml + 2, ny + 2*pml + 2]*2 + [1, 1])
                topo_ixiy = topoz(2:size(topoz, 1):2, 2:size(topoz, 2):2)

                topox = meshgrid([nx + 2*pml + 2, ny + 2*pml + 2], [dx, dy], [-(pml + 1)*dx + ox, -(pml + 1)*dy + oy], dim=1)
                topoy = meshgrid([nx + 2*pml + 2, ny + 2*pml + 2], [dx, dy], [-(pml + 1)*dx + ox, -(pml + 1)*dy + oy], dim=2)

                if (measure_source_depth_from_surface) then
                    src_below_free_surface = trues(gmtr(i)%ns)
                else
                    stp = ginterp(topox, topoy, flatten(-topo_ixiy), gmtr(i)%srcr(:)%x, gmtr(i)%srcr(:)%y)
                    src_below_free_surface = gmtr(i)%srcr(:)%z >= stp
                end if

                if (measure_receiver_depth_from_surface) then
                    rec_below_free_surface = trues(gmtr(i)%nr)
                else
                    rtp = ginterp(topox, topoy, flatten(-topo_ixiy), gmtr(i)%recr(:)%x, gmtr(i)%recr(:)%y)
                    rec_below_free_surface = gmtr(i)%recr(:)%z >= rtp
                end if

                topo_max = maxval(topo_ixiy)

#endif

            else

                src_below_free_surface = trues(gmtr(i)%ns)
                rec_below_free_surface = trues(gmtr(i)%nr)
                topo_max = 0

            end if

            ! Check source validity
            l = 0
            do j = 1, gmtr(i)%ns
                if ( &
                        gmtr(i)%srcr(j)%x >= xmin .and. gmtr(i)%srcr(j)%x <= xmax .and. &
                        gmtr(i)%srcr(j)%y >= ymin .and. gmtr(i)%srcr(j)%y <= ymax .and. &
                        gmtr(i)%srcr(j)%z >= zmin - topo_max .and. gmtr(i)%srcr(j)%z <= zmax - topo_max .and. &
                        gmtr(i)%srcr(j)%x >= sxmin .and. gmtr(i)%srcr(j)%x <= sxmax .and. &
                        gmtr(i)%srcr(j)%y >= symin .and. gmtr(i)%srcr(j)%y <= symax .and. &
                        gmtr(i)%srcr(j)%z - topo_max >= szmin .and. gmtr(i)%srcr(j)%z <= szmax - topo_max .and. &
                        gmtr(i)%id >= sid_min .and. gmtr(i)%id <= sid_max .and. &
                        src_below_free_surface(j)) then
                    l = l + 1
                else
                    gmtr(i)%srcr(j)%amp = 0
                end if
            end do

            if (l == 0) then
                gmtr(i)%ns = 0
            end if

            ! Check receiver validity
            l = 0
            do j = 1, gmtr(i)%nr

                ! Compute absolute offset
                gmtr(i)%recr(j)%aoff = sqrt( &
                    (gmtr(i)%recr(j)%x - mean(gmtr(i)%srcr(:)%x))**2 &
                    + (gmtr(i)%recr(j)%y - mean(gmtr(i)%srcr(:)%y))**2 &
                    + (gmtr(i)%recr(j)%z - mean(gmtr(i)%srcr(:)%z))**2)

                ! Check if receiver in a proper spatial range
                if ( &
                        gmtr(i)%recr(j)%x >= xmin .and. gmtr(i)%recr(j)%x <= xmax .and. &
                        gmtr(i)%recr(j)%y >= ymin .and. gmtr(i)%recr(j)%y <= ymax .and. &
                        gmtr(i)%recr(j)%z >= zmin - topo_max .and. gmtr(i)%recr(j)%z <= zmax - topo_max .and. &
                        gmtr(i)%recr(j)%x >= rxmin .and. gmtr(i)%recr(j)%x <= rxmax .and. &
                        gmtr(i)%recr(j)%y >= rymin .and. gmtr(i)%recr(j)%y <= rymax .and. &
                        gmtr(i)%recr(j)%z >= rzmin - topo_max .and. gmtr(i)%recr(j)%z <= rzmax - topo_max .and. &
                        gmtr(i)%recr(j)%aoff >= offset_min .and. gmtr(i)%recr(j)%aoff <= offset_max .and. &
                        j >= rec_min .and. j <= rec_max .and. mod(j - rec_min, rec_every) == 0 .and. &
                        .not. any(rec_exclude == j) .and. &
                        gmtr(i)%recr(j)%weight /= 0 .and. rec_below_free_surface(j)) then
                    l = l + 1
                else
                    gmtr(i)%recr(j)%weight = 0.0
                end if

            end do

            ! Check if any receiver in the model
            if (l == 0) then
                gmtr(i)%nr = 0
            end if

            ! Display progress
            if ((mod(i, max(nint(ns/10.0), 1)) == 0 .or. i == ns .or. i == 1) .and. rankid == 0) then
                call warn(date_time_compact()//' Loading geometry '//num2str(i)//' of '//num2str(ns))
            end if

        end do
        close (11)

        ! Remove common-shot gathers that are not qualified
        qs = falses(ns)
        ! Select sources based on offset, gain and representational id restrictions
        do i = shot_min, shot_max, shot_every
            if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0) then
                qs(i) = .true.
            end if
        end do
        if (allocated(src_select)) then
            if (any(src_select /= 0)) then
                ! Select sources based on their representional id
                qs = .false.
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. any(src_select == i)) then
                        qs(i) = .true.
                    end if
                end do
            end if
        end if
        if (allocated(sid_select)) then
            if (any(sid_select /= 0)) then
                ! Select source based on their field id
                qs = .false.
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. any(sid_select == gmtr(i)%id)) then
                        qs(i) = .true.
                    end if
                end do
            end if
        end if
        if (allocated(src_exclude)) then
            if (any(src_exclude /= 0)) then
                ! Exclude sources based on their representional id
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. qs(i) .and. any(src_exclude == i)) then
                        qs(i) = .false.
                    end if
                end do
            end if
        end if
        if (allocated(sid_exclude)) then
            if (any(sid_exclude /= 0)) then
                ! Exclude source based on their field id
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. qs(i) .and. any(sid_exclude == gmtr(i)%id)) then
                        qs(i) = .false.
                    end if
                end do
            end if
        end if

        ! If no qualified shots, then stop
        if (.not. any(qs)) then
            call warn(date_time_compact()//' @'//get_hostname()//' Error: No shot or receiver in the model. ')
            stop
        end if

        ! Select all qualified shots
        ! Simpler way is gmtr = pack(gmtr, mask=qs), but may not work depending on compiler's version
        allocate (g(1:ns))
        l = 1
        do i = 1, ns
            if (qs(i)) then
                g(l) = gmtr(i)
                l = l + 1
            end if
        end do

        gmtr = g(1:l - 1)
        deallocate (g)

        ns = size(gmtr)

        ! If any duplicate source id, then stop
        if (ns >= 2) then
            call alloc_array(sid_select, [1, ns])
            do i = 1, ns
                sid_select(i) = gmtr(i)%id
            end do
            usid = unique(sid_select)
            if (size(usid) < ns) then
                call warn(date_time_compact()//' <load_geometry> Erorr: Duplicate source ID found. ')
                stop
            end if
        end if

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' Geometry is loaded with '//num2str(ns)//' sources. ')
        end if

        call mpibarrier

    end subroutine

    !
    !> Prepare geometry
    !
    subroutine prepare_geometry(this)

        class(source_receiver_geometry), intent(inout) :: this
        logical :: depth_varying_dz

        real :: dis
        integer :: i, igx, igy

        if (this%ns == 0 .or. this%nr == 0) then
            call warn(' Error: No shot or receiver in the model. ')
            stop
        end if

        depth_varying_dz = allocated(this%z_i) .and. allocated(this%z_h) .and. rov(this%z_i) > 0 .and. rov(this%z_h) > 0

        ! Source
        !        !$omp parallel do private(i, j, dis, igx, igy)
        do i = 1, this%ns

            this%srcr(i)%valid = this%srcr(i)%x >= this%xmin &
                .and. this%srcr(i)%x <= this%xmax &
                .and. this%srcr(i)%y >= this%ymin &
                .and. this%srcr(i)%y <= this%ymax &
                .and. this%srcr(i)%z >= this%zmin &
                .and. this%srcr(i)%z <= this%zmax &
                .and. this%srcr(i)%x >= this%sxmin &
                .and. this%srcr(i)%x <= this%sxmax &
                .and. this%srcr(i)%y >= this%symin &
                .and. this%srcr(i)%y <= this%symax &
                .and. this%srcr(i)%z >= this%szmin &
                .and. this%srcr(i)%z <= this%szmax

            if (this%srcr(i)%valid) then

                ! Find the nearest integer grid point for the source: [i - 1/2, i + 1/2) E i
                this%srcr(i)%gx = nint((this%srcr(i)%x - this%ox)/this%dx) + 1
                if (this%dy > 0) then
                    this%srcr(i)%gy = nint((this%srcr(i)%y - this%oy)/this%dy) + 1
                end if

                if (depth_varying_dz) then
                    ! Note that here minloc will always treat lower bound of array as 1, not true lower bound
                    igx = ifelse(size(this%z_i, 1) == 1, 1, this%srcr(i)%gx)
                    igy = ifelse(size(this%z_i, 2) == 1, 1, this%srcr(i)%gy)
                    this%srcr(i)%gz = minloc(abs(this%srcr(i)%z - this%oz - this%z_i(igx, igy, :)), dim=1) + lbound(this%z_i, dim=3) - 1
                else
                    this%srcr(i)%gz = nint((this%srcr(i)%z - this%oz)/this%dz) + 1
                end if

                ! Find the nearest half integer grid point for the source: [i, i + 1) E i + 1/2
                this%srcr(i)%hx = nint((this%srcr(i)%x - this%ox + 0.5*this%dx)/this%dx) + 1
                if (this%dy > 0) then
                    this%srcr(i)%hy = nint((this%srcr(i)%y - this%oy + 0.5*this%dy)/this%dy) + 1
                end if

                if (depth_varying_dz) then
                    igx = ifelse(size(this%z_h, 1) == 1, 1, this%srcr(i)%gx)
                    igy = ifelse(size(this%z_h, 2) == 1, 1, this%srcr(i)%gy)
                    ! For z = 0, we want the relevant source/receiver to be below the free surface
                    this%srcr(i)%hz = minloc(abs(this%srcr(i)%z + ifelse(this%srcr(i)%z == 0, float_small, 0.0) &
                        - this%oz - this%z_h(igx, igy, :)), dim=1) + lbound(this%z_h, dim=3) - 1
                else
                    this%srcr(i)%hz = nint((this%srcr(i)%z - this%oz + 0.5*this%dz)/this%dz) + 1
                end if

                ! Allocate memory for Kaiser-windowed sinc interpolators
                call alloc_array(this%srcr(i)%interp_ix, [-nkw, nkw])
                call alloc_array(this%srcr(i)%interp_iy, [-nkw, nkw])
                call alloc_array(this%srcr(i)%interp_iz, [-nkw, nkw])
                call alloc_array(this%srcr(i)%interp_hx, [-nkw, nkw])
                call alloc_array(this%srcr(i)%interp_hy, [-nkw, nkw])
                call alloc_array(this%srcr(i)%interp_hz, [-nkw, nkw])

                ! Integer finite-diference grids
                ! Check: x = 0.1, then gx = 1, dis = 0.1 - (1 - 1) = 0.1
                ! x = 0.5, then gx = 2, dis = 0.5 - (2 - 1) = -0.5
                ! x = 0.7, then gx = 2, dis = 0.7 - (2 - 1) = -0.3
                dis = (this%srcr(i)%x - this%ox)/this%dx - (this%srcr(i)%gx - 1)
                this%srcr(i)%interp_ix = kaiser_sinc_kernel(dis)

                if (this%dy > 0) then
                    dis = (this%srcr(i)%y - this%oy)/this%dy - (this%srcr(i)%gy - 1)
                    this%srcr(i)%interp_iy = kaiser_sinc_kernel(dis)
                end if

                if (depth_varying_dz) then
                    igx = ifelse(size(this%z_i, 1) == 1, 1, this%srcr(i)%gx)
                    igy = ifelse(size(this%z_i, 2) == 1, 1, this%srcr(i)%gy)
                    dis = this%srcr(i)%z - this%oz - this%z_i(igx, igy, this%srcr(i)%gz)
                    dis = dis/this%dz_i(igx, igy, this%srcr(i)%gz)
                else
                    dis = (this%srcr(i)%z - this%oz)/this%dz - (this%srcr(i)%gz - 1)
                end if
                this%srcr(i)%interp_iz = kaiser_sinc_kernel(dis)

                ! Half-integer finite-diference grids
                ! In all the modeling codes, the half-integer to integer conversion convention is
                ! -0.5 --> 0, 0.5 --> 1, 1.5 --> 2, ...
                ! That is, using nint(float)
                ! Here, the integer grid points have been computed to be nint(float) + 1
                ! Therefore, for example, assume we want to get the Kaiser-windowed sinc interpolator
                ! for a point at x = 0.5, then it is +0.5 to 0, but is +0 to the variable assigned at
                ! the half-integer points (vx in acoustic-wave modeling, for instance)
                ! Therefore, the interpolator is 0, 0, 0, 0, 1, 0, 0, 0, 0
                ! x = 0.1, then hgx = 2 (1/2), then dis = 0.1 + 0.5 - (2 - 1) = -0.4
                ! x = 0.5, then hgx = 2 (1/2), then dis = 0.5 + 0.5 - (2 - 1) = 0,
                ! x = 0.9, then hgx = 2 (1/2), then dis = 0.9 + 0.5 - (2 - 1) = 0.4
                !
                dis = (this%srcr(i)%x - this%ox + 0.5*this%dx)/this%dx - (this%srcr(i)%hx - 1)
                this%srcr(i)%interp_hx = kaiser_sinc_kernel(dis)

                if (this%dy > 0) then
                    dis = (this%srcr(i)%y - this%oy + 0.5*this%dy)/this%dy - (this%srcr(i)%hy - 1)
                    this%srcr(i)%interp_hy = kaiser_sinc_kernel(dis)
                end if

                if (depth_varying_dz) then
                    igx = ifelse(size(this%z_h, 1) == 1, 1, this%srcr(i)%gx)
                    igy = ifelse(size(this%z_h, 2) == 1, 1, this%srcr(i)%gy)
                    dis = this%srcr(i)%z - this%oz - this%z_h(igx, igy, this%srcr(i)%hz)
                    dis = dis/this%dz_h(igx, igy, this%srcr(i)%hz)
                else
                    dis = (this%srcr(i)%z - this%oz + 0.5*this%dz)/this%dz - (this%srcr(i)%hz - 1)
                end if
                this%srcr(i)%interp_hz = kaiser_sinc_kernel(dis)

            end if

        end do
        !        !$omp end parallel do

        ! Receiver
        !        !$omp parallel do private(i, dis, igx, igy)
        do i = 1, this%nr

            this%recr(i)%valid = this%recr(i)%x >= this%xmin &
                .and. this%recr(i)%x <= this%xmax &
                .and. this%recr(i)%y >= this%ymin &
                .and. this%recr(i)%y <= this%ymax &
                .and. this%recr(i)%z >= this%zmin &
                .and. this%recr(i)%z <= this%zmax &
                .and. this%recr(i)%x >= this%rxmin &
                .and. this%recr(i)%x <= this%rxmax &
                .and. this%recr(i)%y >= this%rymin &
                .and. this%recr(i)%y <= this%rymax &
                .and. this%recr(i)%z >= this%rzmin &
                .and. this%recr(i)%z <= this%rzmax

            if (this%recr(i)%weight /= 0 .and. this%recr(i)%valid) then

                ! Find the nearest integer and half-integer grid points
                this%recr(i)%gx = nint((this%recr(i)%x - this%ox)/this%dx) + 1
                if (this%dy > 0) then
                    this%recr(i)%gy = nint((this%recr(i)%y - this%oy)/this%dy) + 1
                end if

                if (depth_varying_dz) then
                    ! Note that here minloc will always treat lower bound of array as 1, not true lower bound
                    igx = ifelse(size(this%z_i, 1) == 1, 1, this%recr(i)%gx)
                    igy = ifelse(size(this%z_i, 2) == 1, 1, this%recr(i)%gy)
                    this%recr(i)%gz = minloc(abs(this%recr(i)%z - this%oz - this%z_i(igx, igy, :)), dim=1) + lbound(this%z_i, dim=3) - 1
                else
                    this%recr(i)%gz = nint((this%recr(i)%z - this%oz)/this%dz) + 1
                end if

                this%recr(i)%hx = nint((this%recr(i)%x - this%ox + 0.5*this%dx)/this%dx) + 1
                if (this%dy > 0) then
                    this%recr(i)%hy = nint((this%recr(i)%y - this%oy + 0.5*this%dy)/this%dy) + 1
                end if

                if (depth_varying_dz) then
                    igx = ifelse(size(this%z_h, 1) == 1, 1, this%recr(i)%gx)
                    igy = ifelse(size(this%z_h, 2) == 1, 1, this%recr(i)%gz)
                    ! Note that here minloc will always treat lower bound of array as 1, not true lower bound
                    this%recr(i)%hz = minloc(abs(this%recr(i)%z  + ifelse(this%recr(i)%z == 0, float_small, 0.0) &
                        - this%oz - this%z_h(igx, igy, :)), dim=1) + lbound(this%z_h, dim=3) - 1
                else
                    this%recr(i)%hz = nint((this%recr(i)%z - this%oz + 0.5*this%dz)/this%dz) + 1
                end if

                ! For interpolation
                call alloc_array(this%recr(i)%interp_ix, [-nkw, nkw])
                call alloc_array(this%recr(i)%interp_iy, [-nkw, nkw])
                call alloc_array(this%recr(i)%interp_iz, [-nkw, nkw])
                call alloc_array(this%recr(i)%interp_hx, [-nkw, nkw])
                call alloc_array(this%recr(i)%interp_hy, [-nkw, nkw])
                call alloc_array(this%recr(i)%interp_hz, [-nkw, nkw])

                ! Integer finite-difference grids
                dis = (this%recr(i)%x - this%ox)/this%dx - (this%recr(i)%gx - 1)
                this%recr(i)%interp_ix = kaiser_sinc_kernel(dis)

                if (this%dy > 0) then
                    dis = (this%recr(i)%y - this%oy)/this%dy - (this%recr(i)%gy - 1)
                    this%recr(i)%interp_iy = kaiser_sinc_kernel(dis)
                end if

                if (depth_varying_dz) then
                    igx = ifelse(size(this%z_i, 1) == 1, 1, this%recr(i)%gx)
                    igy = ifelse(size(this%z_i, 2) == 1, 1, this%recr(i)%gy)
                    dis = this%recr(i)%z - this%oz - this%z_i(igx, igy, this%recr(i)%gz)
                    dis = dis/this%dz_i(igx, igy, this%recr(i)%gz)
                else
                    dis = (this%recr(i)%z - this%oz)/this%dz - (this%recr(i)%gz - 1)
                end if
                this%recr(i)%interp_iz = kaiser_sinc_kernel(dis)

                ! Half-integer finite-difference grids
                dis = (this%recr(i)%x - this%ox + 0.5*this%dx)/this%dx - (this%recr(i)%hx - 1)
                this%recr(i)%interp_hx = kaiser_sinc_kernel(dis)

                if (this%dy > 0) then
                    dis = (this%recr(i)%y - this%oy + 0.5*this%dy)/this%dy - (this%recr(i)%hy - 1)
                    this%recr(i)%interp_hy = kaiser_sinc_kernel(dis)
                end if

                if (depth_varying_dz) then
                    igx = ifelse(size(this%z_h, 1) == 1, 1, this%recr(i)%gx)
                    igy = ifelse(size(this%z_h, 2) == 1, 1, this%recr(i)%gy)
                    dis = this%recr(i)%z - this%oz - this%z_h(igx, igy, this%recr(i)%hz)
                    dis = dis/this%dz_h(igx, igy, this%recr(i)%hz)
                else
                    dis = (this%recr(i)%z - this%oz + 0.5*this%dz)/this%dz - (this%recr(i)%hz - 1)
                end if
                this%recr(i)%interp_hz = kaiser_sinc_kernel(dis)

            end if

        end do
        !        !$omp end parallel do

    contains

        function kaiser_sinc_kernel(x) result(h)

            ! Target fractional offset
            real, intent(in) :: x

            ! Output interpolation weights
            real, allocatable, dimension(:) :: h

            integer :: i
            real :: t, r

            call alloc_array(h, [-nkw, nkw])

            if (x == 0) then
                h(0) = 1.0
                return
            end if

            do i = -nkw, nkw

                ! Distance from interpolation point
                t = i - x

                ! Sinc
                h(i) = sinc(const_pi*t)

                ! Kaiser-windowed sinc
                if (t >= -nkw .and. t <= nkw) then
                    r = sqrt(1.0 - (t/nkw)**2)
                    r = bessel_i0(kaiser_b0*r)/bessel_i0(kaiser_b0)
                else
                    r = 0.0
                end if

                ! Interpolation coefficients
                h(i) = h(i)*r

            end do

        end function

    end subroutine

    !
    !> Unfiltered sinc wavelet -- may have long tails
    !
    function sinc_wavelet(f0, t) result(wavelet)

        real, intent(in) :: f0, t
        real :: wavelet

        wavelet = sinc(const_pi*f0*t)

    end function sinc_wavelet

    !
    !> Gaussian wavelet
    !
    function gaussian_wavelet(f0, t) result(wavelet)

        real, intent(in) :: f0, t
        real :: wavelet

        wavelet = exp(-(const_pi*f0*t)**2)

    end function gaussian_wavelet

    function gaussian_wavelet_hdur(hdur, t) result(wavelet)

        real :: t, hdur
        real :: wavelet

        real :: a

        a = 1.0/hdur**2
        wavelet = exp(-a*t**2)/(sqrt(const_pi)*hdur)

    end function gaussian_wavelet_hdur

    !
    !> Gaussian first derivative wavelet
    !
    function gaussian_deriv_wavelet(f0, t) result(wavelet)

        real, intent(in) :: f0, t
        real :: wavelet

        wavelet = -t*exp(-(const_pi*f0*t)**2)

    end function gaussian_deriv_wavelet

    !
    !> Gaussian second derivative wavelet, a.k.a. Ricker wavelet
    !
    function ricker_wavelet(f0, t) result(wavelet)

        real, intent(in) :: f0, t
        real :: wavelet

        wavelet = (const_pi*f0*t)**2
        wavelet = (1 - 2.0*wavelet)*exp(-wavelet)

    end function ricker_wavelet

    !
    !> First-order derivative of Ricker wavelet
    !
    function ricker_deriv_wavelet(f0, t) result(wavelet)

        real, intent(in) :: f0, t
        real :: wavelet

        wavelet = (const_pi*f0*t)**2
        wavelet = (2.0*wavelet - 3.0)*exp(-wavelet)*t

    end function ricker_deriv_wavelet

    !
    !> Omsby wavelet
    !
    function ormsby_wavelet(f, t) result(wavelet)

        real, dimension(1:4), intent(in) :: f
        real, intent(in) :: t
        real :: wavelet

        real :: f1, f2, f3, f4

        f1 = f(1)
        f2 = f(2)
        f3 = f(3)
        f4 = f(4)

        wavelet = const_pi &
            *(f4**2/(f4 - f3)*sinc(const_pi*f4*t)**2 &
            - f3**2/(f4 - f3)*sinc(const_pi*f3*t)**2 &
            - f2**2/(f2 - f1)*sinc(const_pi*f2*t)**2 &
            + f1**2/(f2 - f1)*sinc(const_pi*f1*t)**2)

    end function ormsby_wavelet

    !
    !> Morlet wavelet
    !
    function morlet_wavelet(f0, t, n) result(wavelet)

        real :: f0, t
        integer, optional :: n
        real :: wavelet

        integer :: ncycle
        real :: sigma

        if (present(n)) then
            ncycle = n
        else
            ncycle = 3
        end if

        sigma = ncycle/(2*const_pi*f0)
        wavelet = cos(2*const_pi*f0*t)*exp(-t**2/(2*sigma**2))

    end function morlet_wavelet

    !
    !> Cosine wavelet
    !
    function cosine_wavelet(f0, t) result(wavelet)

        real, intent(in) :: f0, t
        real :: wavelet

        wavelet = cos(2.0*const_pi*f0*t)*exp(-2*(f0*t)**2)

    end function cosine_wavelet

    !
    !> Gammatone wavelet
    !
    function gammatone_wavelet(f0, t, n, b, phi) result(wavelet)

        real :: f0, t
        integer, optional :: n
        real, optional :: b, phi
        real :: wavelet

        real :: order, bandwidth, phase

        if (present(n)) then
            order = n*1.0
        else
            order = 3.0
        end if

        if (present(b)) then
            bandwidth = b
        else
            bandwidth = 2*f0
        end if

        if (present(phi)) then
            phase = phi
        else
            phase = 1.0/f0
        end if

        wavelet = t**(order - 1.0)*exp(-2*const_pi*f0*t*bandwidth)*cos(2*const_pi*f0*t + phase)

    end function gammatone_wavelet

    !
    !> Prepare geometry
    !
    subroutine prepare_stf(this)

        class(source_receiver_geometry), intent(inout) :: this

        integer :: i, t, nw, funit

        ! source wavelet
        do i = 1, this%ns

            if (this%srcr(i)%valid) then

                if (this%srcr(i)%wavelet /= 'custom') then
                    ! For analytic wavelet, compute wavelet length and signature

                    ! Number of samples in the source time function
                    if (this%srcr(i)%nt == 0) then
                        select case (this%srcr(i)%wavelet)
                            case ('morlet')
                                ! Morlet wavelet
                                this%srcr(i)%nt = ceiling(this%srcr(i)%morlet_cycle*this%f0factor/this%srcr(i)%f0/this%dt)
                            case ('ormsby')
                                ! Ormsby wavelet has a longer duration
                                this%srcr(i)%nt = ceiling(5.0*this%f0factor/this%srcr(i)%f0/this%dt)
                            case default
                                ! Other types of wavelet has a 2/f0 duration
                                this%srcr(i)%nt = ceiling(2.0*this%f0factor/this%srcr(i)%f0/this%dt)
                        end select
                    end if

                    ! wavelet signature
                    this%srcr(i)%stf = zeros(this%srcr(i)%nt)

                    do t = 1, this%srcr(i)%nt

                        select case (this%srcr(i)%wavelet)
                            case ('gaussian')
                                this%srcr(i)%stf(t) = &
                                    gaussian_wavelet(this%srcr(i)%f0, (t - 1)*this%dt - this%f0factor/this%srcr(i)%f0)
                            case ('gaussian_deriv', 'gaussian_deriv1')
                                this%srcr(i)%stf(t) = &
                                    gaussian_deriv_wavelet(this%srcr(i)%f0, (t - 1)*this%dt - this%f0factor/this%srcr(i)%f0)
                            case ('gaussian_deriv_deriv', 'gaussian_deriv2', 'ricker')
                                this%srcr(i)%stf(t) = &
                                    ricker_wavelet(this%srcr(i)%f0, (t - 1)*this%dt - this%f0factor/this%srcr(i)%f0)
                            case ('gaussian_deriv_deriv_deriv', 'gaussian_deriv3', 'ricker_deriv', 'ricker_deriv1')
                                this%srcr(i)%stf(t) = &
                                    ricker_deriv_wavelet(this%srcr(i)%f0, (t - 1)*this%dt - this%f0factor/this%srcr(i)%f0)
                            case ('morlet')
                                this%srcr(i)%stf(t) = &
                                    morlet_wavelet(this%srcr(i)%f0, (t - 1)*this%dt - this%srcr(i)%morlet_cycle/2.0*this%f0factor/this%srcr(i)%f0, this%srcr(i)%morlet_cycle)
                            case ('cosine')
                                this%srcr(i)%stf(t) = &
                                    cosine_wavelet(this%srcr(i)%f0, (t - 1)*this%dt - this%f0factor/this%srcr(i)%f0)
                            case ('ormsby')
                                this%srcr(i)%stf(t) = &
                                    ormsby_wavelet(this%srcr(i)%ormsby_freqs, (t - 1)*this%dt - 5.0/2.0*this%f0factor/this%srcr(i)%f0)
                            case default
                                this%srcr(i)%stf(t) = &
                                    ricker_wavelet(this%srcr(i)%f0, (t - 1)*this%dt - this%f0factor/this%srcr(i)%f0)
                        end select

                    end do

                    this%srcr(i)%stf = this%srcr(i)%stf/maxval(abs(this%srcr(i)%stf))

                else
                    ! For custom wavelet, read in and resample if necessary

                    nw = count_nonempty_lines(this%srcr(i)%file_custom_stf)
                    call alloc_array(this%srcr(i)%stf_t, [1, nw])
                    call alloc_array(this%srcr(i)%stf_amp, [1, nw])
                    open (newunit=funit, file=tidy(this%srcr(i)%file_custom_stf), action='read', status='old')
                    do t = 1, nw
                        read (funit, *) this%srcr(i)%stf_t(t), this%srcr(i)%stf_amp(t)
                    end do
                    close (funit)
                    if (this%srcr(i)%stf_t(2) - this%srcr(i)%stf_t(1) /= this%dt) then
                        nw = nint((this%srcr(i)%stf_t(nw) - this%srcr(i)%stf_t(1))/this%dt) + 1
                        this%srcr(i)%stf = ginterp(this%srcr(i)%stf_t, this%srcr(i)%stf_amp, regspace(0.0, this%dt, (nw - 1)*this%dt), method='sinc')
                        this%srcr(i)%nt = nw
                    end if

                end if

                ! Filter the source time function in the frequency domain
                if (allocated(this%filtering_freqs) .and. allocated(this%filtering_coefs)) then
                    if (this%filtering_freqs(1) /= -1 .and. this%filtering_coefs(1) /= -1) then
                        ! Do the filtering
                        this%srcr(i)%stf = fourier_filt(this%srcr(i)%stf, this%dt, this%filtering_freqs, this%filtering_coefs)
                        ! Tapering head and tail to avoid nonzero values
                        this%srcr(i)%stf = taper(this%srcr(i)%stf, len=[nint(this%srcr(i)%nt*0.05), nint(this%srcr(i)%nt*0.05)])
                    end if
                end if

                ! Half length of the source time function
                this%srcr(i)%hnt = nint(this%srcr(i)%nt/2.0)

                ! Integrate or compute derivative of source wavelet if necessary
                ! For instance, for explosive source or moment tensor source, the mt is
                ! added to the stress component, therefore compute the derivative of the signature
                ! See:
                !       Pitarka, 1999, 3D Elastic Finite-Difference Modeling of Seismic Motion Using
                !       Staggered Grids with Nonuniform Spacing
                !
                if (this%srcr(i)%time_integration > 0) then
                    do t = 1, this%srcr(i)%time_integration
                        this%srcr(i)%stf = integ(this%srcr(i)%stf)*this%dt
                    end do
                end if
                if (this%srcr(i)%time_integration < 0) then
                    do t = 1, abs(this%srcr(i)%time_integration)
                        this%srcr(i)%stf = deriv(this%srcr(i)%stf)/this%dt
                    end do
                end if

                if (this%srcr(i)%time_integration /= 0) then
                    if (abs(this%srcr(i)%stf(this%srcr(i)%nt))/maxval(abs(this%srcr(i)%stf)) > 0.1) then
                        ! If after integration, the last sample has a non-trivial amplitude,
                        ! then pad the wavelet to the entire simulation length
                        call pad_array(this%srcr(i)%stf, [0, this%nt - this%srcr(i)%nt])
                        this%srcr(i)%nt = this%nt
                    end if
                end if

            end if

        end do

    end subroutine

    !
    !> Add a source value a, spread by the interpolation weights wx(irx)*wz(irz), to the 2D field w
    !> around the grid point (gx, gz).
    !>
    !> With a free surface, the source comes with its mirror image of value a_image (method of
    !> images): stencil rows above the surface act through their mirror rows below it, and on the
    !> surface row the source and its image coincide. Nothing is written above the surface. Fields
    !> on integer z rows have the surface on row 1 and mirror row r to 2 - r; fields on half-integer
    !> z rows (half = .true.) have the surface between rows 1 and 2 and mirror row r to 3 - r.
    !>
    !> On the vertically refined mesh of free-surface modeling, a and a_image are densities on a
    !> cell of the nominal height dz, and the value on each row r is scaled by zscale(r) = dz/dz(r)
    !> (zscale is used if present and allocated), so that each row receives the density of its own
    !> cell height. Together with receivers that interpolate with the same weights, this makes the
    !> injection the adjoint of the receiver sampling (source-receiver reciprocity).
    !
    subroutine add_source_value_2d(w, gx, gz, wx, wz, a, a_image, half, free_surface, zscale)

        real, allocatable, dimension(:, :), intent(inout) :: w
        integer, intent(in) :: gx, gz
        real, dimension(-nkw:nkw), intent(in) :: wx, wz
        real, intent(in) :: a, a_image
        logical, intent(in) :: half, free_surface
        real, allocatable, dimension(:), intent(in), optional :: zscale

        integer :: irx, irz, r, r1, rr
        real :: ar, ai
        logical :: scaled

        scaled = .false.
        if (present(zscale)) scaled = allocated(zscale)

        ! First row on or below the free surface, and the sum of a row and its mirror row
        r1 = merge(2, 1, half)
        rr = merge(3, 2, half)

        do irz = -nkw, nkw
            r = gz + irz
            if (.not. free_surface .or. r >= r1) then
                ar = a
                if (scaled) ar = a*zscale(r)
                do irx = -nkw, nkw
                    w(gx + irx, r) = w(gx + irx, r) + ar*wx(irx)*wz(irz)
                end do
            end if
            if (free_surface .and. rr - r >= r1) then
                ai = a_image
                if (scaled) ai = a_image*zscale(rr - r)
                do irx = -nkw, nkw
                    w(gx + irx, rr - r) = w(gx + irx, rr - r) + ai*wx(irx)*wz(irz)
                end do
            end if
        end do

    end subroutine add_source_value_2d

    !
    !> The 3D version of add_source_value_2d; only points inside the owned index range
    !> block = [x1, x2, y1, y2, z1, z2] of this MPI rank are written.
    !
    subroutine add_source_value_3d(w, gx, gy, gz, wx, wy, wz, a, a_image, half, free_surface, block, zscale)

        real, allocatable, dimension(:, :, :), intent(inout) :: w
        integer, intent(in) :: gx, gy, gz
        real, dimension(-nkw:nkw), intent(in) :: wx, wy, wz
        real, intent(in) :: a, a_image
        logical, intent(in) :: half, free_surface
        integer, dimension(6), intent(in) :: block
        real, allocatable, dimension(:), intent(in), optional :: zscale

        integer :: irx, iry, irz, i, j, r, r1, rr
        real :: ar, ai
        logical :: scaled, direct, image

        scaled = .false.
        if (present(zscale)) scaled = allocated(zscale)

        ! First row on or below the free surface, and the sum of a row and its mirror row
        r1 = merge(2, 1, half)
        rr = merge(3, 2, half)

        do irz = -nkw, nkw
            r = gz + irz
            direct = (.not. free_surface .or. r >= r1) .and. r >= block(5) .and. r <= block(6)
            image = free_surface .and. rr - r >= r1 .and. rr - r >= block(5) .and. rr - r <= block(6)
            ar = a
            ai = a_image
            if (scaled .and. direct) ar = a*zscale(r)
            if (scaled .and. image) ai = a_image*zscale(rr - r)
            do iry = -nkw, nkw
                j = gy + iry
                do irx = -nkw, nkw
                    i = gx + irx
                    if (i < block(1) .or. i > block(2) .or. j < block(3) .or. j > block(4)) cycle
                    if (direct) then
                        w(i, j, r) = w(i, j, r) + ar*wx(irx)*wy(iry)*wz(irz)
                    end if
                    if (image) then
                        w(i, j, rr - r) = w(i, j, rr - r) + ai*wx(irx)*wy(iry)*wz(irz)
                    end if
                end do
            end do
        end do

    end subroutine add_source_value_3d

    !
    !> Transpose of add_source_value_2d: the sums of the 2D field w weighted by the interpolation
    !> weights wx(irx)*wz(irz) (and zscale) over the points where add_source_value_2d adds the source
    !> (v(1)) and its mirror image (v(2)), as needed for the gradient with respect to source parameters
    !
    function source_value_2d(w, gx, gz, wx, wz, half, free_surface, zscale) result(v)

        real, allocatable, dimension(:, :), intent(in) :: w
        integer, intent(in) :: gx, gz
        real, dimension(-nkw:nkw), intent(in) :: wx, wz
        logical, intent(in) :: half, free_surface
        real, allocatable, dimension(:), intent(in), optional :: zscale
        real, dimension(2) :: v

        integer :: irx, irz, r, r1, rr
        real :: sr, si
        logical :: scaled

        scaled = .false.
        if (present(zscale)) scaled = allocated(zscale)

        r1 = merge(2, 1, half)
        rr = merge(3, 2, half)

        v = 0
        do irz = -nkw, nkw
            r = gz + irz
            if (.not. free_surface .or. r >= r1) then
                sr = 1.0
                if (scaled) sr = zscale(r)
                do irx = -nkw, nkw
                    v(1) = v(1) + w(gx + irx, r)*sr*wx(irx)*wz(irz)
                end do
            end if
            if (free_surface .and. rr - r >= r1) then
                si = 1.0
                if (scaled) si = zscale(rr - r)
                do irx = -nkw, nkw
                    v(2) = v(2) + w(gx + irx, rr - r)*si*wx(irx)*wz(irz)
                end do
            end if
        end do

    end function source_value_2d

    !
    !> The 3D version of source_value_2d; only points inside the owned index range
    !> block = [x1, x2, y1, y2, z1, z2] of this MPI rank are summed
    !
    function source_value_3d(w, gx, gy, gz, wx, wy, wz, half, free_surface, block, zscale) result(v)

        real, allocatable, dimension(:, :, :), intent(in) :: w
        integer, intent(in) :: gx, gy, gz
        real, dimension(-nkw:nkw), intent(in) :: wx, wy, wz
        logical, intent(in) :: half, free_surface
        integer, dimension(6), intent(in) :: block
        real, allocatable, dimension(:), intent(in), optional :: zscale
        real, dimension(2) :: v

        integer :: irx, iry, irz, i, j, r, r1, rr
        real :: sr, si
        logical :: scaled, direct, image

        scaled = .false.
        if (present(zscale)) scaled = allocated(zscale)

        r1 = merge(2, 1, half)
        rr = merge(3, 2, half)

        v = 0
        do irz = -nkw, nkw
            r = gz + irz
            direct = (.not. free_surface .or. r >= r1) .and. r >= block(5) .and. r <= block(6)
            image = free_surface .and. rr - r >= r1 .and. rr - r >= block(5) .and. rr - r <= block(6)
            sr = 1.0
            si = 1.0
            if (scaled .and. direct) sr = zscale(r)
            if (scaled .and. image) si = zscale(rr - r)
            do iry = -nkw, nkw
                j = gy + iry
                do irx = -nkw, nkw
                    i = gx + irx
                    if (i < block(1) .or. i > block(2) .or. j < block(3) .or. j > block(4)) cycle
                    if (direct) then
                        v(1) = v(1) + w(i, j, r)*sr*wx(irx)*wy(iry)*wz(irz)
                    end if
                    if (image) then
                        v(2) = v(2) + w(i, j, rr - r)*si*wx(irx)*wy(iry)*wz(irz)
                    end if
                end do
            end do
        end do

    end function source_value_3d

end module
