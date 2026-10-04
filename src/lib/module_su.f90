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


module mod_su

    use libflit
    use iso_fortran_env

    implicit none

    !
    !> SU trace header
    !
    type su_trace_header

        ! byte 1
        integer(4) :: TraceSequenceLine
        integer(4) :: TraceSequenceFile
        integer(4) :: FieldRecordNumber
        integer(4) :: TraceNumber
        integer(4) :: EnergySourcePoint
        integer(4) :: cdp
        integer(4) :: cdpTrace
        integer(2) :: TraceIdenitifactionCode
        integer(2) :: NSummedTraces
        integer(2) :: NStackedTraces
        integer(2) :: DataUse
        ! byte 41
        integer(4) :: offset
        integer(4) :: ReceiverGroupElevation
        integer(4) :: SourceSurfaceElevation
        integer(4) :: SourceDepth
        integer(4) :: ReceiverDatumElevation
        integer(4) :: SourceDatumElevation
        integer(4) :: SourceWaterDepth
        integer(4) :: GroupWaterDepth
        integer(2) :: ElevationScalar
        integer(2) :: SourceGroupScalar
        integer(4) :: SourceX
        integer(4) :: SourceY
        integer(4) :: GroupX
        integer(4) :: GroupY
        ! byte 89
        integer(2) :: CoordinateUnits
        integer(2) :: WeatheringVelocity
        integer(2) :: SubWeatheringVelocity
        integer(2) :: SourceUpholeTime
        integer(2) :: GroupUpholeTime
        integer(2) :: SourceStaticCorrection
        integer(2) :: GroupStaticCorrection
        integer(2) :: TotalStaticApplied
        integer(2) :: LagTimeA
        integer(2) :: LagTimeB
        integer(2) :: DelayRecordingTime
        integer(2) :: MuteTimeStart
        integer(2) :: MuteTimeEnd
        ! byte 115
        ! ns is a 4-byte integer in memory; see encode_su_ns for its encoding in the file
        integer(4) :: ns
        integer(2) :: dt
        integer(2) :: GainType
        integer(2) :: InstrumentGainConstant
        integer(2) :: InstrumentInitialGain
        integer(2) :: Correlated
        integer(2) :: SweepFrequenceStart
        integer(2) :: SweepFrequenceEnd
        integer(2) :: SweepLength
        integer(2) :: SweepType
        integer(2) :: SweepTraceTaperLengthStart
        integer(2) :: SweepTraceTaperLengthEnd
        integer(2) :: TaperType
        integer(2) :: AliasFilterFrequency
        integer(2) :: AliasFilterSlope
        integer(2) :: NotchFilterFrequency
        integer(2) :: NotchFilterSlope
        integer(2) :: LowCutFrequency
        integer(2) :: HighCutFrequency
        integer(2) :: LowCutSlope
        integer(2) :: HighCutSlope
        integer(2) :: YearDataRecorded
        integer(2) :: DayOfYear
        integer(2) :: HourOfDay
        integer(2) :: MinuteOfHour
        integer(2) :: SecondOfMinute
        integer(2) :: TimeBaseCode
        integer(2) :: TraceWeightningFactor
        integer(2) :: GeophoneGroupNumberRoll1
        integer(2) :: GeophoneGroupNumberFirstTraceOrigField
        integer(2) :: GeophoneGroupNumberLastTraceOrigField
        integer(2) :: GapSize
        integer(2) :: OverTravel

        ! begin SU/SEGY differences
        ! byte 181
        real(4) :: d1
        real(4) :: f1
        real(4) :: d2
        real(4) :: f2
        real(4) :: ungpow
        real(4) :: unscale
        ! byte 205
        integer(4) :: ntr
        integer(2) :: mark
        integer(2) :: shortpad
        ! byte 213; OWL stores ns > 65535 at bytes 233-236, see encode_su_ns
        integer(2), dimension(1:14) :: unass
        ! end SU/SEGY differences

    contains
        procedure, public :: print => print_su_trace_header

    end type su_trace_header

    type su_trace
        type(su_trace_header) :: header
        real(kind=4), allocatable, dimension(:) :: data
    end type su_trace

    ! Endianess of su file, =-1 little, =0 native, =1 big
    integer :: inendian = 0
    integer :: outendian = 0

    ! SU gather
    type su
        integer :: nr = 0, nt = 0
        real :: dt = 0, ot = 0
        type(su_trace), allocatable, dimension(:) :: trace
        character(len=24) :: sort = 'source'
        integer :: inendian = -1
        integer :: outendian = -1
        integer :: verbose = 0
        logical, allocatable, dimension(:) :: mask
        integer, allocatable, dimension(:, :) :: trace_range_group
        integer, allocatable, dimension(:, :) :: trace_range_global
    contains
        !        procedure, private, pass(this) :: copy_su
        !        generic :: assignment(=) => copy_su
        procedure, public :: init => init_su
        procedure, public :: load => load_su
        procedure, public :: output => su_output_su
        procedure, public :: to_array => su_to_array
        procedure, public :: from_array => array_to_su
        procedure, public :: select_from_mask => select_su_from_mask
        procedure, public :: sort_key => sort_su
        procedure, public :: shift => shift_su
        procedure, public :: resamp => resample_su
        procedure, public :: clean => clean_su
        procedure, public :: min => get_min_su
        procedure, public :: max => get_max_su
        procedure, public :: freqfilt => freqfilt_su
        procedure, public :: stdout => standard_output_su
        procedure, public :: collect_group => collect_su_group
        procedure, public :: get_key => get_key_su
        procedure, public :: zero_foreign_rank_traces_group => zero_foreign_rank_traces_su_group
    end type su

    ! Basic arithmetic operations for SU
    interface operator(+)
        module procedure :: su_add_su
        module procedure :: su_add_int
        module procedure :: int_add_su
        module procedure :: su_add_float
        module procedure :: float_add_su
        module procedure :: su_add_double
        module procedure :: double_add_su
    end interface

    interface operator(-)
        module procedure :: minus_su
        module procedure :: su_minus_su
        module procedure :: su_minus_int
        module procedure :: int_minus_su
        module procedure :: su_minus_float
        module procedure :: float_minus_su
        module procedure :: su_minus_double
        module procedure :: double_minus_su
    end interface

    interface operator(*)
        module procedure :: su_x_int
        module procedure :: int_x_su
        module procedure :: su_x_float
        module procedure :: float_x_su
        module procedure :: su_x_double
        module procedure :: double_x_su
    end interface

    interface operator(/)
        module procedure :: su_divide_int
        module procedure :: int_divide_su
        module procedure :: su_divide_float
        module procedure :: float_divide_su
        module procedure :: su_divide_double
        module procedure :: double_divide_su
    end interface

contains

    subroutine print_su_trace_header(this)

        class(su_trace_header), intent(in) :: this

        print *, 'TraceSequenceLine                       = ', this%TraceSequenceLine
        print *, 'TraceSequenceFile                       = ', this%TraceSequenceFile
        print *, 'FieldRecordNumber                       = ', this%FieldRecordNumber
        print *, 'TraceNumber                             = ', this%TraceNumber
        print *, 'EnergySourcePoint                       = ', this%EnergySourcePoint
        print *, 'cdp                                     = ', this%cdp
        print *, 'cdpTrace                                = ', this%cdpTrace
        print *, 'TraceIdenitifactionCode                 = ', this%TraceIdenitifactionCode
        print *, 'NSummedTraces                           = ', this%NSummedTraces
        print *, 'NStackedTraces                          = ', this%NStackedTraces
        print *, 'DataUse                                 = ', this%DataUse
        print *, 'offset                                  = ', this%offset

        print *, 'ReceiverGroupElevation                  = ', this%ReceiverGroupElevation
        print *, 'SourceSurfaceElevation                  = ', this%SourceSurfaceElevation
        print *, 'SourceDepth                             = ', this%SourceDepth
        print *, 'ReceiverDatumElevation                  = ', this%ReceiverDatumElevation
        print *, 'SourceDatumElevation                    = ', this%SourceDatumElevation
        print *, 'SourceWaterDepth                        = ', this%SourceWaterDepth
        print *, 'GroupWaterDepth                         = ', this%GroupWaterDepth
        print *, 'ElevationScalar                         = ', this%ElevationScalar
        print *, 'SourceGroupScalar                       = ', this%SourceGroupScalar
        print *, 'SourceX                                 = ', this%SourceX
        print *, 'SourceY                                 = ', this%SourceY
        print *, 'GroupX                                  = ', this%GroupX
        print *, 'GroupY                                  = ', this%GroupY

        print *, 'CoordinateUnits                         = ', this%CoordinateUnits
        print *, 'WeatheringVelocity                      = ', this%WeatheringVelocity
        print *, 'SubWeatheringVelocity                   = ', this%SubWeatheringVelocity
        print *, 'SourceUpholeTime                        = ', this%SourceUpholeTime
        print *, 'GroupUpholeTime                         = ', this%GroupUpholeTime
        print *, 'SourceStaticCorrection                  = ', this%SourceStaticCorrection
        print *, 'GroupStaticCorrection                   = ', this%GroupStaticCorrection
        print *, 'TotalStaticApplied                      = ', this%TotalStaticApplied
        print *, 'LagTimeA                                = ', this%LagTimeA
        print *, 'LagTimeB                                = ', this%LagTimeB
        print *, 'DelayRecordingTime                      = ', this%DelayRecordingTime
        print *, 'MuteTimeStart                           = ', this%MuteTimeStart
        print *, 'MuteTimeEnd                             = ', this%MuteTimeEnd

        print *, 'ns                                      = ', this%ns
        print *, 'dt                                      = ', this%dt
        print *, 'GainType                                = ', this%GainType
        print *, 'InstrumentGainConstant                  = ', this%InstrumentGainConstant
        print *, 'InstrumentInitialGain                   = ', this%InstrumentInitialGain
        print *, 'Correlated                              = ', this%Correlated
        print *, 'SweepFrequenceStart                     = ', this%SweepFrequenceStart
        print *, 'SweepFrequenceEnd                       = ', this%SweepFrequenceEnd
        print *, 'SweepLength                             = ', this%SweepLength
        print *, 'SweepType                               = ', this%SweepType
        print *, 'SweepTraceTaperLengthStart              = ', this%SweepTraceTaperLengthStart
        print *, 'SweepTraceTaperLengthEnd                = ', this%SweepTraceTaperLengthEnd
        print *, 'TaperType                               = ', this%TaperType
        print *, 'AliasFilterFrequency                    = ', this%AliasFilterFrequency
        print *, 'AliasFilterSlope                        = ', this%AliasFilterSlope
        print *, 'NotchFilterFrequency                    = ', this%NotchFilterFrequency
        print *, 'NotchFilterSlope                        = ', this%NotchFilterSlope
        print *, 'LowCutFrequency                         = ', this%LowCutFrequency
        print *, 'HighCutFrequency                        = ', this%HighCutFrequency
        print *, 'LowCutSlope                             = ', this%LowCutSlope
        print *, 'HighCutSlope                            = ', this%HighCutSlope
        print *, 'YearDataRecorded                        = ', this%YearDataRecorded
        print *, 'DayOfYear                               = ', this%DayOfYear
        print *, 'HourOfDay                               = ', this%HourOfDay
        print *, 'MinuteOfHour                            = ', this%MinuteOfHour
        print *, 'SecondOfMinute                          = ', this%SecondOfMinute
        print *, 'TimeBaseCode                            = ', this%TimeBaseCode
        print *, 'TraceWeightningFactor                   = ', this%TraceWeightningFactor
        print *, 'GeophoneGroupNumberRoll1                = ', this%GeophoneGroupNumberRoll1
        print *, 'GeophoneGroupNumberFirstTraceOrigField  = ', this%GeophoneGroupNumberFirstTraceOrigField
        print *, 'GeophoneGroupNumberLastTraceOrigField   = ', this%GeophoneGroupNumberLastTraceOrigField
        print *, 'GapSize                                 = ', this%GapSize
        print *, 'OverTravel                              = ', this%OverTravel

        print *, 'd1                                      = ', this%d1
        print *, 'f1                                      = ', this%f1
        print *, 'd2                                      = ', this%d2
        print *, 'f2                                      = ', this%f2
        print *, 'ungpow                                  = ', this%ungpow
        print *, 'unscale                                 = ', this%unscale
        print *, 'ntr                                     = ', this%ntr
        print *, 'mark                                    = ', this%mark
        print *, 'shortpad                                = ', this%shortpad

    end subroutine

    subroutine collect_su_group(this)

        class(su), intent(inout) :: this

        real, allocatable, dimension(:, :) :: d

        if (nrank_group > 1) then
            d = this%to_array()
            call allreduce_array_group(d)
            call this%from_array(d)
        end if

    end subroutine

    subroutine zero_foreign_rank_traces_su_group(this)

        class(su), intent(inout) :: this

        integer :: i

        if (nrank_group > 1) then

            !$omp parallel do private(i)
            do i = 1, this%nr
                if (i < this%trace_range_group(rankid_group, 1) .or. i > this%trace_range_group(rankid_group, 2)) then
                    this%trace(i)%data = 0
                end if
            end do
            !$omp end parallel do
        end if

    end subroutine

    function get_key_su(this, keys) result(chart)

        class(su), intent(in) :: this
        character(len=*), dimension(:), intent(in) :: keys
        real, allocatable, dimension(:, :) :: chart

        integer :: nc, i, j
        real :: xyscale, zscale

        nc = size(keys)
        chart = zeros(this%nr, nc)

        xyscale = 10.0d0**order_of_magnitude(this%trace(i)%header%SourceGroupScalar*1.0d0)
        zscale = 10.0d0**order_of_magnitude(this%trace(i)%header%ElevationScalar*1.0d0)

        !$omp parallel do private(i)
        do i = 1, this%nr

            do j = 1, nc
                select case (keys(j))
                    case ('sx', 'SourceX')
                        chart(i, j) = this%trace(i)%header%SourceX*xyscale
                    case ('sy', 'SourceY')
                        chart(i, j) = this%trace(i)%header%SourceY*xyscale
                    case ('sdepth', 'SourceDepth')
                        chart(i, j) = this%trace(i)%header%SourceDepth*zscale
                    case ('selev', 'SourceDatumElevation')
                        chart(i, j) = this%trace(i)%header%SourceDatumElevation*zscale
                    case ('gx', 'GroupX')
                        chart(i, j) = this%trace(i)%header%GroupX*xyscale
                    case ('gy', 'GroupY')
                        chart(i, j) = this%trace(i)%header%GroupY*xyscale
                    case ('gelev', 'ReceiverGroupElevation')
                        chart(i, j) = this%trace(i)%header%ReceiverGroupElevation*zscale
                    case ('offset', 'off')
                        chart(i, j) = this%trace(i)%header%offset*xyscale
                    case ('ep', 'EnergySourcePoint')
                        chart(i, j) = this%trace(i)%header%EnergySourcePoint
                    case ('trid', 'TraceIdenitifactionCode')
                        chart(i, j) = this%trace(i)%header%TraceIdenitifactionCode
                    case ('fldr', 'FieldRecordNumber')
                        chart(i, j) = this%trace(i)%header%FieldRecordNumber
                    case ('tracf', 'TraceSequenceFile')
                        chart(i, j) = this%trace(i)%header%TraceSequenceFile
                    case ('tracl', 'TraceSequenceLine')
                        chart(i, j) = this%trace(i)%header%TraceSequenceLine
                end select
            end do

        end do
        !$omp end parallel do

    end function get_key_su

    function get_min_su(this) result(ms)

        class(su), intent(in) :: this
        real :: ms

        ms = minval(this%to_array())

    end function get_min_su

    function get_max_su(this) result(ms)

        class(su), intent(in) :: this
        real :: ms

        ms = maxval(this%to_array())

    end function get_max_su

    subroutine clean_su(this)

        class(su), intent(inout) :: this

        integer :: i

        !$omp parallel do private(i)
        do i = 1, this%nr
            this%trace(i)%data = return_normal(this%trace(i)%data)
        end do
        !$omp end parallel do

    end subroutine clean_su

    !
    !> Resample trace data
    !
    subroutine resample_su(this, nnt, ddt, oot)

        class(su), intent(inout) :: this
        integer, intent(in) :: nnt
        real, intent(in) :: ddt
        real, intent(in), optional :: oot

        real :: resample_ot
        integer :: i
        real, allocatable, dimension(:, :) :: w

        if (present(oot)) then
            resample_ot = oot
        else
            resample_ot = 0
        end if

        if (this%nt == nnt .and. this%dt == ddt .and. resample_ot == 0) then
            return
        end if

        w = this%to_array()
        do i = 1, this%nr
            deallocate(this%trace(i)%data)
            this%trace(i)%data = zeros(nnt)
        end do

        !$omp parallel do private(i)
        do i = 1, this%nr

            if (maxval(abs(w(:, i))) > 0) then
                this%trace(i)%data = interp(w(:, i), this%nt, this%dt, 0.0, nnt, ddt, resample_ot, 'cubic')
            end if

            this%trace(i)%header%ns = nnt
            this%trace(i)%header%dt = int(ddt*1.0e6)
            this%trace(i)%header%d1 = ddt

        end do
        !$omp end parallel do

        this%nt = nnt
        this%dt = ddt

    end subroutine resample_su

    !
    !> Filtering the trace data based on given frequency band
    !
    subroutine freqfilt_su(this, f, a)

        class(su), intent(inout) :: this
        real, dimension(:), intent(in) :: f, a

        integer :: i

        call assert(size(f) > 1, ' <freqfilt_su> Error: size(f) must > 1')
        call assert(size(f) == size(a), ' <freqfilt_su> Error: size(f) must = size(a)')
        call assert(all(f >= 0), ' <freqfilt_su> Error: all f must >= 0')
        call assert(all(a >= 0), ' <freqfilt_su> Error: all a must >= 0')

        !$omp parallel do private(i)
        do i = 1, this%nr
            if (norm2(this%trace(i)%data) > 0) then
                this%trace(i)%data = fourier_filt(this%trace(i)%data, this%dt, f, a, method='hann')
            end if
        end do
        !$omp end parallel do

    end subroutine freqfilt_su

    !
    !> Shift trace data in time
    !
    subroutine shift_su(this, tshift)

        class(su), intent(inout) :: this
        real, intent(in) :: tshift

        integer :: st, i

        st = min(nint(abs(tshift)/this%dt), this%nt)

        !$omp parallel do private(i)
        do i = 1, this%nr
            if (tshift > 0) then
                this%trace(i)%data(st + 1:this%nt) = this%trace(i)%data(1:this%nt - st)
                this%trace(i)%data(1:st) = 0.0
            else
                this%trace(i)%data(1:this%nt - st) = this%trace(i)%data(st + 1:this%nt)
                this%trace(i)%data(this%nt - st + 1:this%nt) = 0.0
            end if
        end do
        !$omp end parallel do

    end subroutine shift_su

    !
    !> Sort traces based on keyword
    !
    subroutine sort_su(this, keyword, order)

        class(su), intent(inout) :: this
        character(len=*), intent(in) :: keyword
        integer, intent(in), optional :: order

        type(su) :: sorted_this
        integer :: i
        integer, allocatable, dimension(:) :: keyword_value, keyword_index
        integer :: sort_order

        sorted_this = this

        if (present(order)) then
            sort_order = order
        else
            sort_order = 1
        end if

        keyword_index = zeros(this%nr)
        keyword_value = zeros(this%nr)

        select case (keyword)
            case ('SourceX', 'sx')
                !$omp parallel do private(i)
                do i = 1, this%nr
                    keyword_value(i) = this%trace(i)%header%SourceX
                end do
                !$omp end parallel do
            case ('SourceY', 'sy')
                !$omp parallel do private(i)
                do i = 1, this%nr
                    keyword_value(i) = this%trace(i)%header%SourceY
                end do
                !$omp end parallel do
            case ('SourceElevation', 'selev')
                !$omp parallel do private(i)
                do i = 1, this%nr
                    keyword_value(i) = this%trace(i)%header%SourceSurfaceElevation
                end do
                !$omp end parallel do
            case ('SourceDepth', 'sdepth')
                !$omp parallel do private(i)
                do i = 1, this%nr
                    keyword_value(i) = this%trace(i)%header%SourceDepth
                end do
                !$omp end parallel do
            case ('GroupX', 'ReceiverX', 'gx', 'rx')
                !$omp parallel do private(i)
                do i = 1, this%nr
                    keyword_value(i) = this%trace(i)%header%GroupX
                end do
                !$omp end parallel do
            case ('GroupY', 'ReceiverY', 'gy', 'ry')
                !$omp parallel do private(i)
                do i = 1, this%nr
                    keyword_value(i) = this%trace(i)%header%GroupY
                end do
                !$omp end parallel do
            case ('GroupElevation', 'ReceiverElevation', 'gelev', 'relev')
                !$omp parallel do private(i)
                do i = 1, this%nr
                    keyword_value(i) = this%trace(i)%header%ReceiverGroupElevation
                end do
                !$omp end parallel do
            case ('Offset', 'offset', 'off', 'srdist')
                !$omp parallel do private(i)
                do i = 1, this%nr
                    keyword_value(i) = this%trace(i)%header%offset
                end do
                !$omp end parallel do
        end select

        ! Sort based on keyword values
        call sort_index(keyword_value, keyword_index, sort_order)

        ! Sort and modify
        !$omp parallel do private(i)
        do i = 1, this%nr
            this%trace(i) = sorted_this%trace(keyword_index(i))
        end do
        !$omp end parallel do

    end subroutine sort_su

    !
    !> su + su
    !
    function su_add_su(a, b) result(c)

        type(su), intent(in) :: a, b
        type(su) :: c

        integer :: i

        call assert(a%nt == b%nt, 'Error: To add two SU data, nt must be same.')
        call assert(a%dt == b%dt, 'Error: To add two SU data, dt must be same.')
        call assert(a%nr == b%nr, 'Error: To add two SU data, nr must be same.')

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data + b%trace(i)%data
        end do

    end function su_add_su

    !
    !> su - su
    !
    function su_minus_su(a, b) result(c)

        type(su), intent(in) :: a, b
        type(su) :: c

        integer :: i

        call assert(a%nt == b%nt, 'Error: To subtract two SU data, nt must be same.')
        call assert(a%dt == b%dt, 'Error: To subtract two SU data, dt must be same.')
        call assert(a%nr == b%nr, 'Error: To subtract two SU data, nr must be same.')

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data - b%trace(i)%data
        end do

    end function su_minus_su

    !
    !> -su
    !
    function minus_su(a) result(c)

        type(su), intent(in) :: a
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = -a%trace(i)%data
        end do

    end function minus_su

    !
    !> su + integer
    !
    function su_add_int(a, s) result(c)

        type(su), intent(in) :: a
        integer, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data + s
        end do

    end function su_add_int

    !
    !> integer + su
    !
    function int_add_su(s, a) result(c)

        type(su), intent(in) :: a
        integer, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data + s
        end do

    end function int_add_su

    !
    !> su - integer
    !
    function su_minus_int(a, s) result(c)

        type(su), intent(in) :: a
        integer, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data - s
        end do

    end function su_minus_int

    !
    !> integer - su
    !
    function int_minus_su(s, a) result(c)

        type(su), intent(in) :: a
        integer, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = s - a%trace(i)%data
        end do

    end function int_minus_su

    !
    !> su*integer
    !
    function su_x_int(a, s) result(c)

        type(su), intent(in) :: a
        integer, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data*s
        end do

    end function su_x_int

    !
    !> integer*su
    !
    function int_x_su(s, a) result(c)

        type(su), intent(in) :: a
        integer, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data*s
        end do

    end function int_x_su

    !
    !> su/integer
    !
    function su_divide_int(a, s) result(c)

        type(su), intent(in) :: a
        integer, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data/s
        end do

    end function su_divide_int

    !
    !> zz
    !
    function int_divide_su(s, a) result(c)

        type(su), intent(in) :: a
        integer, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = s/a%trace(i)%data
        end do

    end function int_divide_su

    function su_add_float(a, s) result(c)

        type(su), intent(in) :: a
        real, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data + s
        end do

    end function su_add_float

    function float_add_su(s, a) result(c)

        type(su), intent(in) :: a
        real, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data + s
        end do

    end function float_add_su

    function su_minus_float(a, s) result(c)

        type(su), intent(in) :: a
        real, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data - s
        end do

    end function su_minus_float

    function float_minus_su(s, a) result(c)

        type(su), intent(in) :: a
        real, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = s - a%trace(i)%data
        end do

    end function float_minus_su

    function su_x_float(a, s) result(c)

        type(su), intent(in) :: a
        real, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data*s
        end do

    end function su_x_float

    function float_x_su(s, a) result(c)

        type(su), intent(in) :: a
        real, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data*s
        end do

    end function float_x_su

    function su_divide_float(a, s) result(c)

        type(su), intent(in) :: a
        real, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data/s
        end do

    end function su_divide_float

    function float_divide_su(s, a) result(c)

        type(su), intent(in) :: a
        real, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = s/a%trace(i)%data
        end do

    end function float_divide_su

    function su_add_double(a, s) result(c)

        type(su), intent(in) :: a
        double precision, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data + s
        end do

    end function su_add_double

    function double_add_su(s, a) result(c)

        type(su), intent(in) :: a
        double precision, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data + s
        end do

    end function double_add_su

    function su_minus_double(a, s) result(c)

        type(su), intent(in) :: a
        double precision, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data - s
        end do

    end function su_minus_double

    function double_minus_su(s, a) result(c)

        type(su), intent(in) :: a
        double precision, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = s - a%trace(i)%data
        end do

    end function double_minus_su

    function su_x_double(a, s) result(c)

        type(su), intent(in) :: a
        double precision, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data*s
        end do

    end function su_x_double

    function double_x_su(s, a) result(c)

        type(su), intent(in) :: a
        double precision, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data*s
        end do

    end function double_x_su

    function su_divide_double(a, s) result(c)

        type(su), intent(in) :: a
        double precision, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = a%trace(i)%data/s
        end do

    end function su_divide_double

    function double_divide_su(s, a) result(c)

        type(su), intent(in) :: a
        double precision, intent(in) :: s
        type(su) :: c

        integer :: i

        c = a

        do i = 1, c%nr
            c%trace(i)%data = s/a%trace(i)%data
        end do

    end function double_divide_su

    function select_su_from_mask(this) result(s)

        class(su), intent(in) :: this
        type(su) :: s

        integer :: nr, i, l

        ! First copy
        s = this

        ! Then assign
        nr = count(this%mask)
        s%nr = nr
        deallocate (s%trace)
        allocate (s%trace(s%nr))

        l = 1
        do i = 1, this%nr
            if (this%mask(i)) then
                s%trace(l) = this%trace(i)
                l = l + 1
            end if
        end do

    end function select_su_from_mask

    subroutine array_to_su(this, x)

        class(su), intent(inout) :: this
        real, dimension(:, :), intent(in) :: x

        integer :: nt, nr, i

        nt = clip(size(x, 1), 1, this%nt)
        nr = clip(size(x, 2), 1, this%nr)

        do i = 1, nr
            this%trace(i)%data = x(:, i)
            this%trace(i)%header%ns = nt
            !            this%trace(i)%header%nsample = nt
        end do

    end subroutine array_to_su

    function su_to_array(this) result(x)

        class(su), intent(in) :: this
        real, allocatable, dimension(:, :) :: x

        integer :: i

        x = zeros(this%nt, this%nr)
        do i = 1, this%nr
            x(:, i) = this%trace(i)%data
        end do

    end function su_to_array

    !
    !> Initialize an SU
    !
    subroutine init_su(this, nt, dt, nr)

        class(su), intent(inout) :: this
        integer, intent(in), optional :: nt
        real, intent(in), optional :: dt
        integer, intent(in), optional :: nr

        integer :: i

        if (present(nt)) then
            this%nt = nt
        end if
        if (present(dt)) then
            this%dt = dt
        end if
        if (present(nr)) then
            this%nr = nr
        end if

        call assert(this%nt >= 1, ' <init_su> Error: nt must >= 1')
        call assert(this%dt > 0, ' <init_su> Error: dt must > 0')
        call assert(this%nr >=1, '<init_su> Error: nr must >= 1')

        ! Allocate memory for traces
        if (allocated(this%trace)) then
            deallocate (this%trace)
        end if
        allocate (this%trace(1:this%nr))

        do i = 1, this%nr
            this%trace(i)%header%ns = this%nt
            !            this%trace(i)%header%nsample = this%nt
            this%trace(i)%header%dt = int(this%dt*1.0e6)
            this%trace(i)%header%d1 = this%dt
            allocate (this%trace(i)%data(1:this%nt))
            this%trace(i)%data = 0.0
        end do

        if (allocated(this%mask)) then
            deallocate (this%mask)
        end if
        allocate (this%mask(1:this%nr))
        this%mask = .true.

        ! Divide traces into groups or in global domain
        call alloc_array(this%trace_range_group, [0, nrank_group - 1, 1, 2])
        call cut(1, this%nr, nrank_group, this%trace_range_group)

        call alloc_array(this%trace_range_global, [0, nrank - 1, 1, 2])
        call cut(1, this%nr, nrank, this%trace_range_global)

    end subroutine init_su

    !
    !> Encode the number of samples into the SU trace header. SU defines ns at
    !> bytes 115-116 as an unsigned 16-bit integer, so this field holds ns <= 65535.
    !> For a larger ns, this field holds the lower 16 bits of ns, and bytes 233-236,
    !> which are unassigned in SU, hold the full ns as a 4-byte integer.
    !
    pure subroutine encode_su_ns(ns, ns16, ns32)

        integer, intent(in) :: ns
        integer(2), intent(out) :: ns16
        integer(4), intent(out) :: ns32

        integer :: low

        ! Lower 16 bits of ns, with the bit pattern of an unsigned 16-bit integer
        low = iand(ns, 65535)
        if (low > 32767) then
            low = low - 65536
        end if
        ns16 = int(low, kind=2)

        if (ns > 65535) then
            ns32 = ns
        else
            ns32 = 0
        end if

    end subroutine encode_su_ns

    !
    !> Decode the number of samples from the SU trace header; inverse of encode_su_ns
    !
    pure function decode_su_ns(ns16, ns32) result(ns)

        integer(2), intent(in) :: ns16
        integer(4), intent(in) :: ns32
        integer :: ns

        ns = iand(int(ns16), 65535)

        ! Bytes 233-236 count only when they agree with bytes 115-116, because
        ! SU files from other programs may hold unrelated values there
        if (ns32 > 65535 .and. iand(ns32, 65535) == ns) then
            ns = ns32
        end if

    end function decode_su_ns

    subroutine load_su(this, infile, nt, dt, nr, header_only, stdin)

        class(su), intent(inout) :: this
        character(len=*), intent(in) :: infile
        integer, intent(in), optional :: nt
        real, intent(in), optional :: dt
        integer, intent(in), optional :: nr
        logical, intent(in), optional :: header_only, stdin

        integer(kind=2) :: ns16, ds
        integer(kind=4) :: ns32
        integer :: ns, funit, i
        integer(kind=8) :: trcbegpos
        logical :: read_header_only
        logical :: from_stdin

        if (allocated(this%mask)) then
            deallocate (this%mask)
        end if
        if (allocated(this%trace)) then
            deallocate (this%trace)
        end if
        if (present(stdin)) then
            from_stdin = stdin
        else
            from_stdin = .false.
        end if

        ! Check file existence
        if (from_stdin) then
            funit = input_unit
        else
            call assert(file_exists(infile), ' <load_su> Error: Input file '//tidy(infile)//' not found. ')
        end if

        ! Read textual header
        if (from_stdin) then
            ! get ns and ntr information
            select case (this%inendian)
                case (1)
                    open (funit, access='stream', status='old', action='read', &
                        convert='big_endian')
                case (-1)
                    open (funit, access='stream', status='old', action='read', &
                        convert='little_endian')
                case (0)
                    open (funit, access='stream', status='old', action='read', &
                        convert='native')
            end select
        else
            select case (this%inendian)
                case (1)
                    open (newunit=funit, file=tidy(infile), access='stream', status='old', action='read', &
                        convert='big_endian')
                case (-1)
                    open (newunit=funit, file=tidy(infile), access='stream', status='old', action='read', &
                        convert='little_endian')
                case (0)
                    open (newunit=funit, file=tidy(infile), access='stream', status='old', action='read', &
                        convert='native')
            end select
        end if


        ! Number of samples -- assume all traces have same number of samples
        if (present(nt)) then
            this%nt = nt
        else
            read (funit, pos=115) ns16
            read (funit, pos=233) ns32
            ns = decode_su_ns(ns16, ns32)
            if (this%nt == 0) then
                this%nt = ns
            else
                if (this%nt /= ns) then
                    call warn(' <load_su> Warning: Inconsistent number of samples, trace header nt = '//num2str(ns)//', preset nt = '//num2str(this%nt))
                end if
            end if
        end if
        call assert(this%nt >= 1, ' <load_su> Error: nt must >= 1.')

        ! Sampling interval -- assume all traces have same value
        if (present(dt)) then
            this%dt = dt
        else
            read (funit, pos=117) ds
            if (this%dt == 0) then
                this%dt = ds/1.0e6
            else
                if (nint(this%dt*1e6) /= ds) then
                    call warn(' <load_su> Warning: Inconsistent sample interval, trace header dt = '//num2str(ds/1.0e6)//', preset dt = '//num2str(this%dt))
                end if
            end if
        end if
        call assert(this%dt > 0, ' <load_su> Error: dt must >= 0.')

        ! Number of traces
        if (present(nr)) then
            this%nr = nr
        else
            if (this%nr == 0) then
                if (present(stdin)) then
                    if (stdin) then
                        this%nr = floor(get_stdin_size()/(240.0d0 + 4.0d0*this%nt))
                    end if
                else
                    this%nr = floor(get_file_size(infile)/(240.0d0 + 4.0d0*this%nt))
                end if
            end if
        end if
        call assert(this%nr >= 1, ' <load_su> Error: nr must >= 1.')

        if (present(header_only)) then
            read_header_only = header_only
        else
            read_header_only = .false.
        end if

        ! Allocate
        allocate (this%trace(1:this%nr))
        allocate (this%mask(1:this%nr))
        this%mask = .true.

        ! Divide traces into groups or in global domain
        call alloc_array(this%trace_range_group, [0, nrank_group - 1, 1, 2])
        call cut(1, this%nr, nrank_group, this%trace_range_group)

        call alloc_array(this%trace_range_global, [0, nrank - 1, 1, 2])
        call cut(1, this%nr, nrank, this%trace_range_global)

        ! Read
        do i = 1, this%nr

            ! read full header information
            trcbegpos = int((i - 1)*(240 + int(this%nt, 8)*4) + 1, kind=8)
            read (funit, pos=trcbegpos) this%trace(i)%header%TraceSequenceLine
            read (funit, pos=trcbegpos + 4) this%trace(i)%header%TraceSequenceFile
            read (funit, pos=trcbegpos + 8) this%trace(i)%header%FieldRecordNumber
            read (funit, pos=trcbegpos + 12) this%trace(i)%header%TraceNumber
            read (funit, pos=trcbegpos + 16) this%trace(i)%header%EnergySourcePoint
            read (funit, pos=trcbegpos + 20) this%trace(i)%header%cdp
            read (funit, pos=trcbegpos + 24) this%trace(i)%header%cdpTrace
            read (funit, pos=trcbegpos + 28) this%trace(i)%header%TraceIdenitifactionCode
            read (funit, pos=trcbegpos + 30) this%trace(i)%header%NSummedTraces
            read (funit, pos=trcbegpos + 32) this%trace(i)%header%NStackedTraces
            read (funit, pos=trcbegpos + 34) this%trace(i)%header%DataUse
            read (funit, pos=trcbegpos + 36) this%trace(i)%header%offset
            read (funit, pos=trcbegpos + 40) this%trace(i)%header%ReceiverGroupElevation
            read (funit, pos=trcbegpos + 44) this%trace(i)%header%SourceSurfaceElevation
            read (funit, pos=trcbegpos + 48) this%trace(i)%header%SourceDepth
            read (funit, pos=trcbegpos + 52) this%trace(i)%header%ReceiverDatumElevation
            read (funit, pos=trcbegpos + 56) this%trace(i)%header%SourceDatumElevation
            read (funit, pos=trcbegpos + 60) this%trace(i)%header%SourceWaterDepth
            read (funit, pos=trcbegpos + 64) this%trace(i)%header%GroupWaterDepth
            read (funit, pos=trcbegpos + 68) this%trace(i)%header%ElevationScalar
            read (funit, pos=trcbegpos + 70) this%trace(i)%header%SourceGroupScalar
            read (funit, pos=trcbegpos + 72) this%trace(i)%header%SourceX
            read (funit, pos=trcbegpos + 76) this%trace(i)%header%SourceY
            read (funit, pos=trcbegpos + 80) this%trace(i)%header%GroupX
            read (funit, pos=trcbegpos + 84) this%trace(i)%header%GroupY
            read (funit, pos=trcbegpos + 88) this%trace(i)%header%CoordinateUnits
            read (funit, pos=trcbegpos + 90) this%trace(i)%header%WeatheringVelocity
            read (funit, pos=trcbegpos + 92) this%trace(i)%header%SubWeatheringVelocity
            read (funit, pos=trcbegpos + 94) this%trace(i)%header%SourceUpholeTime
            read (funit, pos=trcbegpos + 96) this%trace(i)%header%GroupUpholeTime
            read (funit, pos=trcbegpos + 98) this%trace(i)%header%SourceStaticCorrection
            read (funit, pos=trcbegpos + 100) this%trace(i)%header%GroupStaticCorrection
            read (funit, pos=trcbegpos + 102) this%trace(i)%header%TotalStaticApplied
            read (funit, pos=trcbegpos + 104) this%trace(i)%header%LagTimeA
            read (funit, pos=trcbegpos + 106) this%trace(i)%header%LagTimeB
            read (funit, pos=trcbegpos + 108) this%trace(i)%header%DelayRecordingTime
            read (funit, pos=trcbegpos + 110) this%trace(i)%header%MuteTimeStart
            read (funit, pos=trcbegpos + 112) this%trace(i)%header%MuteTimeEnd
            read (funit, pos=trcbegpos + 114) ns16
            read (funit, pos=trcbegpos + 116) this%trace(i)%header%dt
            read (funit, pos=trcbegpos + 118) this%trace(i)%header%GainType
            read (funit, pos=trcbegpos + 120) this%trace(i)%header%InstrumentGainConstant
            read (funit, pos=trcbegpos + 122) this%trace(i)%header%InstrumentInitialGain
            read (funit, pos=trcbegpos + 124) this%trace(i)%header%Correlated
            read (funit, pos=trcbegpos + 126) this%trace(i)%header%SweepFrequenceStart
            read (funit, pos=trcbegpos + 128) this%trace(i)%header%SweepFrequenceEnd
            read (funit, pos=trcbegpos + 130) this%trace(i)%header%SweepLength
            read (funit, pos=trcbegpos + 132) this%trace(i)%header%SweepType
            read (funit, pos=trcbegpos + 134) this%trace(i)%header%SweepTraceTaperLengthStart
            read (funit, pos=trcbegpos + 136) this%trace(i)%header%SweepTraceTaperLengthEnd
            read (funit, pos=trcbegpos + 138) this%trace(i)%header%TaperType
            read (funit, pos=trcbegpos + 140) this%trace(i)%header%AliasFilterFrequency
            read (funit, pos=trcbegpos + 142) this%trace(i)%header%AliasFilterSlope
            read (funit, pos=trcbegpos + 144) this%trace(i)%header%NotchFilterFrequency
            read (funit, pos=trcbegpos + 146) this%trace(i)%header%NotchFilterSlope
            read (funit, pos=trcbegpos + 148) this%trace(i)%header%LowCutFrequency
            read (funit, pos=trcbegpos + 150) this%trace(i)%header%HighCutFrequency
            read (funit, pos=trcbegpos + 152) this%trace(i)%header%LowCutSlope
            read (funit, pos=trcbegpos + 154) this%trace(i)%header%HighCutSlope
            read (funit, pos=trcbegpos + 156) this%trace(i)%header%YearDataRecorded
            read (funit, pos=trcbegpos + 158) this%trace(i)%header%DayOfYear
            read (funit, pos=trcbegpos + 160) this%trace(i)%header%HourOfDay
            read (funit, pos=trcbegpos + 162) this%trace(i)%header%MinuteOfHour
            read (funit, pos=trcbegpos + 164) this%trace(i)%header%SecondOfMinute
            read (funit, pos=trcbegpos + 166) this%trace(i)%header%TimeBaseCode
            read (funit, pos=trcbegpos + 168) this%trace(i)%header%TraceWeightningFactor
            read (funit, pos=trcbegpos + 170) this%trace(i)%header%GeophoneGroupNumberRoll1
            read (funit, pos=trcbegpos + 172) this%trace(i)%header%GeophoneGroupNumberFirstTraceOrigField
            read (funit, pos=trcbegpos + 174) this%trace(i)%header%GeophoneGroupNumberLastTraceOrigField
            read (funit, pos=trcbegpos + 176) this%trace(i)%header%GapSize
            read (funit, pos=trcbegpos + 178) this%trace(i)%header%OverTravel

            read (funit, pos=trcbegpos + 180) this%trace(i)%header%d1
            read (funit, pos=trcbegpos + 184) this%trace(i)%header%f1
            read (funit, pos=trcbegpos + 188) this%trace(i)%header%d2
            read (funit, pos=trcbegpos + 192) this%trace(i)%header%f2
            read (funit, pos=trcbegpos + 204) this%trace(i)%header%ntr
            !            read (funit, pos=trcbegpos + 208) this%trace(i)%header%nsample
            read (funit, pos=trcbegpos + 232) ns32
            this%trace(i)%header%ns = decode_su_ns(ns16, ns32)

            ! read trace data
            if (.not. read_header_only) then
                allocate (this%trace(i)%data(1:this%nt))
                read (funit, pos=trcbegpos + 240) this%trace(i)%data
            end if

            !            ! Overwrite
            !            this%trace(i)%header%ntr = this%nr
            !            this%trace(i)%header%ns = int(this%nt, kind=2)
            !            this%trace(:)%header%dt = int(this%dt*1.0e6, kind=2)
            !            this%trace(i)%header%nsample = this%nt
            !            this%trace(:)%header%d1 = this%dt

            if (this%verbose /= 0) then
                if (mod(i, this%verbose) == 0 .or. i == this%nr) then
                    write (error_unit, *) date_time_compact()//' >> Reading trace '//num2str(i)//' of '//num2str(this%nr)
                end if
            end if

        end do
        close (funit)

        if (this%verbose /= 0) then
            write (error_unit, *)
            write (error_unit, *) date_time_compact()//' Traces input finished'
        end if

        !        if (present(dt)) then
        !            this%dt = dt
        !            this%trace(:)%header%dt = int(dt*1.0e6, kind=2)
        !            this%trace(:)%header%d1 = dt
        !        else
        !            i = 1
        !            if (this%trace(i)%header%dt /= 0 .and. this%trace(i)%header%d1 == 0) then
        !                this%trace(:)%header%d1 = this%trace(:)%header%dt/1.0e6
        !            else if (this%trace(i)%header%dt == 0 .and. this%trace(i)%header%d1 /= 0) then
        !                this%trace(:)%header%dt = int(this%trace(:)%header%d1*1.0e6)
        !                ! int(this%trace(:)%header%d1*1.0e6, kind=2)
        !            else if (this%trace(i)%header%dt == 0 .and. this%trace(i)%header%d1 == 0) then
        !                call warn(date_time_compact()//' <load_su> Warning: Sampling interval set to 1e-3 s.')
        !                this%trace(:)%header%dt = 1000
        !                this%trace(:)%header%d1 = 1.0e-3
        !            end if
        !        end if
        !        this%dt = this%trace(1)%header%d1

    end subroutine load_su

    subroutine su_output_su(this, outfile, select, append)

        class(su), intent(in) :: this
        character(len=*), intent(in) :: outfile
        integer, dimension(:), intent(in), optional :: select
        logical, intent(in), optional :: append

        integer :: ntr, funit, i, ir
        integer(8) :: trcbegpos, existing_size
        integer(2) :: ns16
        integer(4) :: ns32
        integer, allocatable, dimension(:) :: trcrange
        character(len=24) :: convert_type
        logical :: append_to_existing

        if (present(append)) then
            append_to_existing = append
        else
            append_to_existing = .false.
        end if

        ! Endianness of the file
        if (present(select)) then
            trcrange = select
        else
            trcrange = regspace(1, 1, this%nr)
        end if
        ntr = size(trcrange)
        call assert(ntr >= 1, ' <su_output_su> Error: Number of output traces <= 0. ')

        ! Endianness of the file
        select case (this%outendian)
            case (-1)
                convert_type = 'little_endian'
            case (1)
                convert_type = 'big_endian'
            case (0)
                convert_type = 'native'
        end select
        if (append_to_existing) then
            open (newunit=funit, file=tidy(outfile), access='stream', status='old', action='write', &
                convert=tidy(convert_type), position='append')
            existing_size = get_file_size(outfile)
        else
            open (newunit=funit, file=tidy(outfile), access='stream', status='replace', action='write', &
                convert=tidy(convert_type))
            existing_size = 0
        end if

        ! Write SU traces
        do ir = 1, ntr

            i = trcrange(ir)

            ! Trace start byte location
            trcbegpos = existing_size + int(int(ir - 1, kind=8)*(240 + int(this%nt, kind=8)*4) + 1, kind=8)

            ! Number of samples in the header fields
            call encode_su_ns(this%trace(i)%header%ns, ns16, ns32)

            ! Write trace header
            write (funit, pos=trcbegpos) this%trace(i)%header%TraceSequenceLine
            write (funit, pos=trcbegpos + 4) this%trace(i)%header%TraceSequenceFile
            write (funit, pos=trcbegpos + 8) this%trace(i)%header%FieldRecordNumber
            write (funit, pos=trcbegpos + 12) this%trace(i)%header%TraceNumber
            write (funit, pos=trcbegpos + 16) this%trace(i)%header%EnergySourcePoint
            write (funit, pos=trcbegpos + 20) this%trace(i)%header%cdp
            write (funit, pos=trcbegpos + 24) this%trace(i)%header%cdpTrace
            write (funit, pos=trcbegpos + 28) this%trace(i)%header%TraceIdenitifactionCode
            write (funit, pos=trcbegpos + 30) this%trace(i)%header%NSummedTraces
            write (funit, pos=trcbegpos + 32) this%trace(i)%header%NStackedTraces
            write (funit, pos=trcbegpos + 34) this%trace(i)%header%DataUse
            write (funit, pos=trcbegpos + 36) this%trace(i)%header%offset
            write (funit, pos=trcbegpos + 40) this%trace(i)%header%ReceiverGroupElevation
            write (funit, pos=trcbegpos + 44) this%trace(i)%header%SourceSurfaceElevation
            write (funit, pos=trcbegpos + 48) this%trace(i)%header%SourceDepth
            write (funit, pos=trcbegpos + 52) this%trace(i)%header%ReceiverDatumElevation
            write (funit, pos=trcbegpos + 56) this%trace(i)%header%SourceDatumElevation
            write (funit, pos=trcbegpos + 60) this%trace(i)%header%SourceWaterDepth
            write (funit, pos=trcbegpos + 64) this%trace(i)%header%GroupWaterDepth
            write (funit, pos=trcbegpos + 68) this%trace(i)%header%ElevationScalar
            write (funit, pos=trcbegpos + 70) this%trace(i)%header%SourceGroupScalar
            write (funit, pos=trcbegpos + 72) this%trace(i)%header%SourceX
            write (funit, pos=trcbegpos + 76) this%trace(i)%header%SourceY
            write (funit, pos=trcbegpos + 80) this%trace(i)%header%GroupX
            write (funit, pos=trcbegpos + 84) this%trace(i)%header%GroupY
            write (funit, pos=trcbegpos + 88) this%trace(i)%header%CoordinateUnits
            write (funit, pos=trcbegpos + 90) this%trace(i)%header%WeatheringVelocity
            write (funit, pos=trcbegpos + 92) this%trace(i)%header%SubWeatheringVelocity
            write (funit, pos=trcbegpos + 94) this%trace(i)%header%SourceUpholeTime
            write (funit, pos=trcbegpos + 96) this%trace(i)%header%GroupUpholeTime
            write (funit, pos=trcbegpos + 98) this%trace(i)%header%SourceStaticCorrection
            write (funit, pos=trcbegpos + 100) this%trace(i)%header%GroupStaticCorrection
            write (funit, pos=trcbegpos + 102) this%trace(i)%header%TotalStaticApplied
            write (funit, pos=trcbegpos + 104) this%trace(i)%header%LagTimeA
            write (funit, pos=trcbegpos + 106) this%trace(i)%header%LagTimeB
            write (funit, pos=trcbegpos + 108) this%trace(i)%header%DelayRecordingTime
            write (funit, pos=trcbegpos + 110) this%trace(i)%header%MuteTimeStart
            write (funit, pos=trcbegpos + 112) this%trace(i)%header%MuteTimeEnd
            write (funit, pos=trcbegpos + 114) ns16
            write (funit, pos=trcbegpos + 116) this%trace(i)%header%dt
            write (funit, pos=trcbegpos + 118) this%trace(i)%header%GainType
            write (funit, pos=trcbegpos + 120) this%trace(i)%header%InstrumentGainConstant
            write (funit, pos=trcbegpos + 122) this%trace(i)%header%InstrumentInitialGain
            write (funit, pos=trcbegpos + 124) this%trace(i)%header%Correlated
            write (funit, pos=trcbegpos + 126) this%trace(i)%header%SweepFrequenceStart
            write (funit, pos=trcbegpos + 128) this%trace(i)%header%SweepFrequenceEnd
            write (funit, pos=trcbegpos + 130) this%trace(i)%header%SweepLength
            write (funit, pos=trcbegpos + 132) this%trace(i)%header%SweepType
            write (funit, pos=trcbegpos + 134) this%trace(i)%header%SweepTraceTaperLengthStart
            write (funit, pos=trcbegpos + 136) this%trace(i)%header%SweepTraceTaperLengthEnd
            write (funit, pos=trcbegpos + 138) this%trace(i)%header%TaperType
            write (funit, pos=trcbegpos + 140) this%trace(i)%header%AliasFilterFrequency
            write (funit, pos=trcbegpos + 142) this%trace(i)%header%AliasFilterSlope
            write (funit, pos=trcbegpos + 144) this%trace(i)%header%NotchFilterFrequency
            write (funit, pos=trcbegpos + 146) this%trace(i)%header%NotchFilterSlope
            write (funit, pos=trcbegpos + 148) this%trace(i)%header%LowCutFrequency
            write (funit, pos=trcbegpos + 150) this%trace(i)%header%HighCutFrequency
            write (funit, pos=trcbegpos + 152) this%trace(i)%header%LowCutSlope
            write (funit, pos=trcbegpos + 154) this%trace(i)%header%HighCutSlope
            write (funit, pos=trcbegpos + 156) this%trace(i)%header%YearDataRecorded
            write (funit, pos=trcbegpos + 158) this%trace(i)%header%DayOfYear
            write (funit, pos=trcbegpos + 160) this%trace(i)%header%HourOfDay
            write (funit, pos=trcbegpos + 162) this%trace(i)%header%MinuteOfHour
            write (funit, pos=trcbegpos + 164) this%trace(i)%header%SecondOfMinute
            write (funit, pos=trcbegpos + 166) this%trace(i)%header%TimeBaseCode
            write (funit, pos=trcbegpos + 168) this%trace(i)%header%TraceWeightningFactor
            write (funit, pos=trcbegpos + 170) this%trace(i)%header%GeophoneGroupNumberRoll1
            write (funit, pos=trcbegpos + 172) this%trace(i)%header%GeophoneGroupNumberFirstTraceOrigField
            write (funit, pos=trcbegpos + 174) this%trace(i)%header%GeophoneGroupNumberLastTraceOrigField
            write (funit, pos=trcbegpos + 176) this%trace(i)%header%GapSize
            write (funit, pos=trcbegpos + 178) this%trace(i)%header%OverTravel

            write (funit, pos=trcbegpos + 180) this%trace(i)%header%d1
            write (funit, pos=trcbegpos + 184) this%trace(i)%header%f1
            write (funit, pos=trcbegpos + 188) this%trace(i)%header%d2
            write (funit, pos=trcbegpos + 192) this%trace(i)%header%f2

            write (funit, pos=trcbegpos + 204) ntr
            !            write (funit, pos=trcbegpos + 208) this%nt
            write (funit, pos=trcbegpos + 232) ns32

            ! Write trace data
            write (funit, pos=trcbegpos + 240) this%trace(i)%data

            if (this%verbose /= 0) then
                if(mod(ir, min(this%verbose, ntr)) == 0) then
                    write (error_unit, *) date_time_compact()//' >> Writing trace '//num2str(ir) &
                        //' of '//num2str(ntr)
                end if
            end if

        end do
        close (funit)

        if (this%verbose /= 0 .and. ntr >= 1) then
            write (error_unit, *) date_time_compact()//' Traces output finished'
        end if

    end subroutine su_output_su

    !
    !> Write data and header information to a SU file
    !
    subroutine standard_output_su(this, select, append)

        class(su), intent(in) :: this
        integer, dimension(:), intent(in), optional :: select
        logical, intent(in), optional :: append

        integer :: ntr, i, ir
        integer(8) :: trcbegpos
        integer(2) :: ns16
        integer(4) :: ns32
        integer, allocatable, dimension(:) :: trcrange
        character(len=24) :: convert_type
        logical :: append_to_existing

        if (present(append)) then
            append_to_existing = append
        else
            append_to_existing = .false.
        end if

        ! Endianness of the file
        if (present(select)) then
            trcrange = select
        else
            trcrange = regspace(1, 1, this%nr)
        end if
        ntr = size(trcrange)
        call assert(ntr >= 1, ' <su_output_su> Error: Number of output traces <= 0. ')

        ! Endianness of the file
        select case (this%outendian)
            case (-1)
                convert_type = 'little_endian'
            case (1)
                convert_type = 'big_endian'
            case (0)
                convert_type = 'native'
        end select

        close (output_unit)
        if (append_to_existing) then
            open (output_unit, access='stream', status='old', action='write', &
                convert=tidy(convert_type), position='append')
        else
            open (output_unit, access='stream', status='replace', action='write', &
                convert=tidy(convert_type))
        end if

        ! Write SU traces
        do ir = 1, ntr

            i = trcrange(ir)

            ! Trace start byte location
            trcbegpos = int(int(ir - 1, kind=8)*(240 + int(this%nt, kind=8)*4) + 1, kind=8)

            ! Number of samples in the header fields
            call encode_su_ns(this%trace(i)%header%ns, ns16, ns32)

            ! Write trace header
            write (output_unit, pos=trcbegpos) this%trace(i)%header%TraceSequenceLine
            write (output_unit, pos=trcbegpos + 4) this%trace(i)%header%TraceSequenceFile
            write (output_unit, pos=trcbegpos + 8) this%trace(i)%header%FieldRecordNumber
            write (output_unit, pos=trcbegpos + 12) this%trace(i)%header%TraceNumber
            write (output_unit, pos=trcbegpos + 16) this%trace(i)%header%EnergySourcePoint
            write (output_unit, pos=trcbegpos + 20) this%trace(i)%header%cdp
            write (output_unit, pos=trcbegpos + 24) this%trace(i)%header%cdpTrace
            write (output_unit, pos=trcbegpos + 28) this%trace(i)%header%TraceIdenitifactionCode
            write (output_unit, pos=trcbegpos + 30) this%trace(i)%header%NSummedTraces
            write (output_unit, pos=trcbegpos + 32) this%trace(i)%header%NStackedTraces
            write (output_unit, pos=trcbegpos + 34) this%trace(i)%header%DataUse
            write (output_unit, pos=trcbegpos + 36) this%trace(i)%header%offset
            write (output_unit, pos=trcbegpos + 40) this%trace(i)%header%ReceiverGroupElevation
            write (output_unit, pos=trcbegpos + 44) this%trace(i)%header%SourceSurfaceElevation
            write (output_unit, pos=trcbegpos + 48) this%trace(i)%header%SourceDepth
            write (output_unit, pos=trcbegpos + 52) this%trace(i)%header%ReceiverDatumElevation
            write (output_unit, pos=trcbegpos + 56) this%trace(i)%header%SourceDatumElevation
            write (output_unit, pos=trcbegpos + 60) this%trace(i)%header%SourceWaterDepth
            write (output_unit, pos=trcbegpos + 64) this%trace(i)%header%GroupWaterDepth
            write (output_unit, pos=trcbegpos + 68) this%trace(i)%header%ElevationScalar
            write (output_unit, pos=trcbegpos + 70) this%trace(i)%header%SourceGroupScalar
            write (output_unit, pos=trcbegpos + 72) this%trace(i)%header%SourceX
            write (output_unit, pos=trcbegpos + 76) this%trace(i)%header%SourceY
            write (output_unit, pos=trcbegpos + 80) this%trace(i)%header%GroupX
            write (output_unit, pos=trcbegpos + 84) this%trace(i)%header%GroupY
            write (output_unit, pos=trcbegpos + 88) this%trace(i)%header%CoordinateUnits
            write (output_unit, pos=trcbegpos + 90) this%trace(i)%header%WeatheringVelocity
            write (output_unit, pos=trcbegpos + 92) this%trace(i)%header%SubWeatheringVelocity
            write (output_unit, pos=trcbegpos + 94) this%trace(i)%header%SourceUpholeTime
            write (output_unit, pos=trcbegpos + 96) this%trace(i)%header%GroupUpholeTime
            write (output_unit, pos=trcbegpos + 98) this%trace(i)%header%SourceStaticCorrection
            write (output_unit, pos=trcbegpos + 100) this%trace(i)%header%GroupStaticCorrection
            write (output_unit, pos=trcbegpos + 102) this%trace(i)%header%TotalStaticApplied
            write (output_unit, pos=trcbegpos + 104) this%trace(i)%header%LagTimeA
            write (output_unit, pos=trcbegpos + 106) this%trace(i)%header%LagTimeB
            write (output_unit, pos=trcbegpos + 108) this%trace(i)%header%DelayRecordingTime
            write (output_unit, pos=trcbegpos + 110) this%trace(i)%header%MuteTimeStart
            write (output_unit, pos=trcbegpos + 112) this%trace(i)%header%MuteTimeEnd
            write (output_unit, pos=trcbegpos + 114) ns16
            write (output_unit, pos=trcbegpos + 116) this%trace(i)%header%dt
            write (output_unit, pos=trcbegpos + 118) this%trace(i)%header%GainType
            write (output_unit, pos=trcbegpos + 120) this%trace(i)%header%InstrumentGainConstant
            write (output_unit, pos=trcbegpos + 122) this%trace(i)%header%InstrumentInitialGain
            write (output_unit, pos=trcbegpos + 124) this%trace(i)%header%Correlated
            write (output_unit, pos=trcbegpos + 126) this%trace(i)%header%SweepFrequenceStart
            write (output_unit, pos=trcbegpos + 128) this%trace(i)%header%SweepFrequenceEnd
            write (output_unit, pos=trcbegpos + 130) this%trace(i)%header%SweepLength
            write (output_unit, pos=trcbegpos + 132) this%trace(i)%header%SweepType
            write (output_unit, pos=trcbegpos + 134) this%trace(i)%header%SweepTraceTaperLengthStart
            write (output_unit, pos=trcbegpos + 136) this%trace(i)%header%SweepTraceTaperLengthEnd
            write (output_unit, pos=trcbegpos + 138) this%trace(i)%header%TaperType
            write (output_unit, pos=trcbegpos + 140) this%trace(i)%header%AliasFilterFrequency
            write (output_unit, pos=trcbegpos + 142) this%trace(i)%header%AliasFilterSlope
            write (output_unit, pos=trcbegpos + 144) this%trace(i)%header%NotchFilterFrequency
            write (output_unit, pos=trcbegpos + 146) this%trace(i)%header%NotchFilterSlope
            write (output_unit, pos=trcbegpos + 148) this%trace(i)%header%LowCutFrequency
            write (output_unit, pos=trcbegpos + 150) this%trace(i)%header%HighCutFrequency
            write (output_unit, pos=trcbegpos + 152) this%trace(i)%header%LowCutSlope
            write (output_unit, pos=trcbegpos + 154) this%trace(i)%header%HighCutSlope
            write (output_unit, pos=trcbegpos + 156) this%trace(i)%header%YearDataRecorded
            write (output_unit, pos=trcbegpos + 158) this%trace(i)%header%DayOfYear
            write (output_unit, pos=trcbegpos + 160) this%trace(i)%header%HourOfDay
            write (output_unit, pos=trcbegpos + 162) this%trace(i)%header%MinuteOfHour
            write (output_unit, pos=trcbegpos + 164) this%trace(i)%header%SecondOfMinute
            write (output_unit, pos=trcbegpos + 166) this%trace(i)%header%TimeBaseCode
            write (output_unit, pos=trcbegpos + 168) this%trace(i)%header%TraceWeightningFactor
            write (output_unit, pos=trcbegpos + 170) this%trace(i)%header%GeophoneGroupNumberRoll1
            write (output_unit, pos=trcbegpos + 172) this%trace(i)%header%GeophoneGroupNumberFirstTraceOrigField
            write (output_unit, pos=trcbegpos + 174) this%trace(i)%header%GeophoneGroupNumberLastTraceOrigField
            write (output_unit, pos=trcbegpos + 176) this%trace(i)%header%GapSize
            write (output_unit, pos=trcbegpos + 178) this%trace(i)%header%OverTravel

            write (output_unit, pos=trcbegpos + 180) this%trace(i)%header%d1
            write (output_unit, pos=trcbegpos + 184) this%trace(i)%header%f1
            write (output_unit, pos=trcbegpos + 188) this%trace(i)%header%d2
            write (output_unit, pos=trcbegpos + 192) this%trace(i)%header%f2

            write (output_unit, pos=trcbegpos + 204) ntr
            !            write (output_unit, pos=trcbegpos + 208) this%nt
            write (output_unit, pos=trcbegpos + 232) ns32

            ! Write trace data
            write (output_unit, pos=trcbegpos + 240) this%trace(i)%data

            if (this%verbose /= 0) then
                if(mod(ir, min(this%verbose, ntr)) == 0) then
                    write (error_unit, *) date_time_compact()//' >> Writing trace '//num2str(ir) &
                        //' of '//num2str(ntr)
                end if
            end if

        end do
        close (output_unit)

        if (this%verbose /= 0 .and. ntr >= 1) then
            write (error_unit, *) date_time_compact()//' Traces output finished'
        end if

    end subroutine standard_output_su

end module
