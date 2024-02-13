!>  \brief  Useful functions which are not methods and don't belong anywhere else
module UsefulFunctions
    use Globals
    implicit none

    interface FileSize
        module procedure FileSizeFromFilename
        module procedure FileSizeFromUnitNumber
    end interface

    interface RadiansToDegrees
        module procedure DegreesSingle
        module procedure DegreesDouble
    end interface

    interface DegreesToRadians
        module procedure RadiansSingle
        module procedure RadiansDouble
    end interface

    contains

        !>  Convert degrees to radians
        elemental pure real function RadiansSingle(degrees)
            real, intent(in) :: degrees
            RadiansSingle = degrees * pi / 180.
        end function RadiansSingle

        !>  convert radians to degrees
        elemental pure real function DegreesSingle(radians)
            real, intent(in) :: radians
            DegreesSingle = radians / pi * 180.
        end function DegreesSingle

        !>  convert degrees to radians (double precision)
        pure elemental function RadiansDouble(degrees)
            real(kind=8), intent(in) :: degrees
            real(kind=8)    ::  RadiansDouble
            RadiansDouble = degrees * dpi / 180.
        end function RadiansDouble

        !>  convert radians to degrees (double precision)
        pure elemental function DegreesDouble(radians)
            real(kind=8), intent(in) :: radians
            real(kind=8)    ::  DegreesDouble
            DegreesDouble = radians / dpi * 180.
        end function DegreesDouble

        pure elemental function DegreesBetween0And360(degrees)
            real(kind=8),   intent(in)  ::  degrees
            real(kind=8)                ::  DegreesBetween0And360
            DegreesBetween0And360 = mod(degrees,360.0d0)
            if (DegreesBetween0And360 .lt. 0.0d0) DegreesBetween0And360 = DegreesBetween0And360 + 360.0d0
        end function DegreesBetween0And360

        !>  \brief  Normalised cross-correlation between two arrays of reals
        !!  \todo   Optimize  - possibly replace with library function from say BLAS or MKL
        pure function NormalizedCrossCorrelation(array_1,array_2) result(cc)
            real(kind=8),   intent(in)  ::  array_1(:)
            real(kind=8),   intent(in)  ::  array_2(:)
            ! result
            real(kind=8)                ::  cc
            ! private variables
            real(kind=8)                ::  array_1_mean, array_2_mean
            real(kind=8)                ::  array_1_sigma, array_2_sigma
            ! start work
            array_1_mean = sum(array_1) / real(size(array_1))
            array_2_mean = sum(array_2) / real(size(array_2))
            array_1_sigma = sum(array_1**2) / real(size(array_1)) - (sum(array_1) / real(size(array_1)))**2
            array_2_sigma = sum(array_2**2) / real(size(array_2)) - (sum(array_2) / real(size(array_2)))**2
            if (array_1_sigma .gt. 0.0d0 .and. array_2_sigma .gt. 0.0d0) then
                cc = sum((array_1-array_1_mean)*(array_2-array_2_mean)) &
                     / sqrt(array_1_sigma*array_2_sigma) / real(size(array_1))
            else
                cc = 0.0d0
            endif
        end function NormalizedCrossCorrelation

        !>  \brief  Rank cross-correlation
        pure function RankCrossCorrelation(array_1,array_2) result(cc)
            use napack_sort2, only : sort2
            ! arguments
            real(kind=8),   intent(in)  ::  array_1(:)
            real(kind=8),   intent(in)  ::  array_2(:)
            ! result
            real(kind=8)                ::  cc
            ! private variables
            real(kind=8)                ::  work_array(size(array_1))
            real(kind=8)                ::  order_1(size(array_1)), order_2(size(array_2))
            real(kind=8)                ::  sorted_array_1(size(array_1)), sorted_array_2(size(array_2))
            ! start work
            sorted_array_1 = array_1
            sorted_array_2 = array_2
            call sort2(sorted_array_1,order_1,work_array,size(array_1))
            call sort2(sorted_array_2,order_2,work_array,size(array_2))

            cc = NormalizedCrossCorrelation(order_1,order_2)
        end function RankCrossCorrelation



        !>  \brief  Check a file exists on disk
        logical function FileExists(filename)
            character(len=*), intent(in)    ::  filename
            inquire(file=trim(adjustl(filename)), exist=FileExists)
        end function FileExists

        !>  \brief  Copy a file by reading in and writing out its content. Do not use for large files.
        subroutine FileCopy(filename_src,filename_dest)
            character(len=*),   intent(in)  ::  filename_src,filename_dest
            ! Private variables
            integer                         ::  lun_src, lun_dest ! logical unit numbers
            integer                         ::  number_of_bytes
            integer(kind=1), allocatable    ::  bytes(:)
            integer                         ::  io_status
            ! Start work
            if (.not. FileExists(filename_src)) then
                call this_program%TerminateWithFatalError('file_copy','Source file does not exist: '//trim(adjustl(filename_src)))
            endif
            ! Open the files
            lun_src = this_program%GetAvailableUnit()
            open(unit=lun_src, file=filename_src, action='read', status='old', access='stream', iostat=io_status)

            if (io_status .ne. 0) then
                close(lun_src)
                call this_program%ReleaseUnit(lun_src)
                call this_program%TerminateWithFatalError('UsefulFunctions::FileCopy', &
                                                  'File exists but cannot be opened: '//trim(adjustl(filename_src)))
            endif

            lun_dest = this_program%GetAvailableUnit()
            open(unit=lun_dest,file=filename_dest,action='write',status='replace',access='stream', iostat=io_status)

            if (io_status .ne. 0) then
                close(lun_dest)
                call this_program%ReleaseUnit(lun_dest)
            endif

            ! Find out size of file to copy
            number_of_bytes = FileSize(lun_src)
            if (number_of_bytes .gt. 0) then
                ! Read file data in
                allocate(bytes(number_of_bytes))
                read(unit=lun_src,pos=1,iostat=io_status) bytes
                if (io_status .ne. 0) then
                    write(*,'(a,i0,2a)') '**error(file_copy): io error ', io_status, ' when reading from: ', filename_src
                    call this_program%TerminateWithFatalError('file_copy','Read error')
                endif
                write(unit=lun_dest,pos=1,iostat=io_status) bytes
                if (io_status .ne. 0) then
                    write(*,'(a,i0,2a)') '**error(file_copy): io error ', io_status, ' when writing to: ', filename_dest
                    call this_program%TerminateWithFatalError('file_copy','Write error')
                endif
            endif

            ! close files..

            if (UnitIsOpen(lun_src)) then
                close(lun_src)
                call this_program%ReleaseUnit(lun_src)
            endif

            if (UnitIsOpen(lun_dest)) then
                close(lun_dest)
                call this_program%ReleaseUnit(lun_dest)
            endif



            ! Deallocate memory
            if (allocated(bytes)) deallocate(bytes)



        end subroutine FileCopy

        !> \brief   Find file size in bytes
        function FileSizeFromFilename(filename) result(file_size)
            character(len=*), intent(in)    ::  filename
            integer(kind=8)                 ::  file_size
            inquire(file=trim(adjustl(filename)),size=file_size)
        end function FileSizeFromFilename

        !> \brief   Find file size in bytes
        function FileSizeFromUnitNumber(lun) result(file_size)
            integer, intent(in)    ::  lun
            integer(kind=8)                 ::  file_size
            inquire(unit=lun,size=file_size)
        end function FileSizeFromUnitNumber

        !>  \brief  Check whether a IO unit is current open
        logical function UnitIsOpen(unit_number)
            integer, intent(in)      ::  unit_number
            integer :: io_status
            character(len=100) :: io_message
            io_status = 0
            inquire(unit=unit_number, opened=UnitIsOpen,iostat=io_status,iomsg=io_message)
            if (io_status .ne. 0) then
                print *, 'UnitIsOpen: IO error ', io_status, ': ', trim(adjustl(io_message))
                call this_program%TerminateWithFatalError('UnitIsOpen','IO error: '//trim(adjustl(io_message)))
            endif
        end function UnitIsOpen

        !>    \brief    returns true if the argument is even
        pure elemental logical function IsEven(int)
            integer,    intent(in)    ::    int
            ! test bit 0 of number. if it is 0, the number is even
            IsEven = .not. btest(int,0)
        end function IsEven

        !>    \brief    returns true if the argument is odd
        pure elemental logical function IsOdd(int)
            integer,    intent(in)    ::    int
            ! start work
            IsOdd = btest(int,0)
        end function IsOdd

        !>  \brief  Compute the phase to be applied to an element distance_from_origin away from the origin so that the real-space shift
        !!          real_space_shift is applied
        pure function ReturnPhaseFromShift(real_space_shift, distance_from_origin, dimension_size)

            ! Arguments

            real,           intent(in)      ::  real_space_shift          !<  Real space shift (in pixels)
            integer,        intent(in)      ::  distance_from_origin      !<  Distance from origin in Fourier transform (in pixels) along the axis of the shift
            integer,        intent(in)      ::  dimension_size            !<  The total extent of the image in the dimension relevant to the shift

            ! Result
            real                            ::  ReturnPhaseFromShift
            ! Private variables

            ! Start work
            ReturnPhaseFromShift = real_space_shift * real(distance_from_origin) * 2.0e0 * pi / real(dimension_size)

        end function ReturnPhaseFromShift

        !>  \brief  Compute the phase shift, given the phase in each dimension
        pure complex function Return3DPhaseFromIndividualDimensions(phase_x, phase_y, phase_z)
            !use LookUpTables
            ! Arguments
            real,   intent(in)        ::  phase_x
            real,   intent(in)        ::  phase_y
            real,   intent(in)        ::  phase_z

            ! Private Variables
            real                      ::  temp_phase
            logical,       parameter  ::  use_look_up_tables = .false.

            ! Start work
            temp_phase  = - phase_x - phase_y - phase_z
            !if (use_look_up_tables) then
            !    Return3DPhaseFromIndividualDimensions = cmplx(cos_lookup(temp_phase), sin_lookup(temp_phase))
            !else
                Return3DPhaseFromIndividualDimensions = cmplx(cos(temp_phase), sin(temp_phase))
            !endif

        end function Return3DPhaseFromIndividualDimensions

        !>  \brief  Sort the elements of an array using the quick sort algorithm
        recursive subroutine QuickSort(array)
            real(kind=8),   intent(inout)   ::  array(:)
            ! private variables
            integer ::  pivot_index
            ! start work
            if (size(array) .gt. 1) then
                call Partition(array,pivot_index)
                call QuickSort(array(:pivot_index-1))
                call QuickSort(array(pivot_index:))
            endif
        end subroutine QuickSort

        !>  \brief  Partition an array so that all elements to "the left" and to "the right" are below and above the pivot point
        pure subroutine Partition(array,pivot_index)
            real(kind=8),   intent(inout)   ::  array(:)
            integer,        intent(out)     ::  pivot_index
            ! private variables
            integer         ::  i,j
            real(kind=8)    ::  temp
            real(kind=8)    ::  pivot_value
            ! start work
            pivot_value = array(1)
            i = 0
            j = size(array) + 1

            do
                j=j-1
                do
                    if(array(j) .le. pivot_value) exit
                    j=j-1
                enddo
                i=i+1
                do
                    if (array(i) .ge. pivot_value) exit
                    i=i+1
                enddo
                if (i .lt. j) then
                    ! exchange
                    temp = array(i)
                    array(i) = array(j)
                    array(j) = temp
                else if (i .eq. j) then
                    pivot_index = i+1
                    return
                else
                    pivot_index = i
                    return
                endif
            enddo
        end subroutine Partition



end module UsefulFunctions
