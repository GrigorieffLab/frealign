module ProgramInstances

    implicit none
    private

    integer,    parameter   ::  lowest_unit_number = 10
    integer,    parameter   ::  highest_unit_number = 200

    !>  \brief Hold data about a program instance
    type, public :: ProgramInstance
        private

        logical, public                         ::  running_interactively     !<  whether program is being run interactively
        logical, public                         ::  reading_from_terminal        !<  whether program is reading input from standard input

        character(len=:), public, allocatable   ::  program_name
        character(len=:), public, allocatable   ::  control_filename

        ! start and finish time
        integer                                 ::  date_start(8)
        integer                                 ::  date_finish(8)


        ! Unit numbers - we keep track of the ones we're currently using
        logical                                 ::  unit_is_available(lowest_unit_number:highest_unit_number) = .true.

        contains
            procedure                           ::  Init
            procedure                           ::  SetStartDate
            procedure                           ::  ParseCommandLineArguments
            procedure                           ::  TerminateWithFatalError
            procedure                           ::  Terminate
            procedure                           ::  GetAvailableUnit
            procedure                           ::  ReleaseUnit
    end type


    contains

    !>  \brief  Find an available unit and reserve it for use
    integer function GetAvailableUnit(self)
        ! arguments
        class(ProgramInstance), intent(inout)   ::  self
        ! private variables
        logical     ::  success
        integer     ::  current_unit
        ! start work
        !$omp critical (ProgramInstance_UnitManagement)
        success = .false.
        do current_unit=lowest_unit_number,highest_unit_number
            if (self%unit_is_available(current_unit)) then
                success = .true.
                GetAvailableUnit = current_unit
                exit
            endif
        enddo
        if (success) then
            self%unit_is_available(GetAvailableUnit) = .false.
        else
            call self%TerminateWithFatalError('ProgramInstance::GetAvailableUnit','Could not find an unused io unit')
        endif
        !$omp end critical (ProgramInstance_UnitManagement)
    end function GetAvailableUnit

    !>  \brief  Let the world know that a particular unit is not being used any more
    subroutine ReleaseUnit(self,unit_number)
        class(ProgramInstance),     intent(inout)   ::  self
        integer,                    intent(in)      ::  unit_number
        ! start work
        !$omp critical (ProgramInstance_UnitManagement)
        if (unit_number .ge. lbound(self%unit_is_available,1) .and. &
            unit_number .le. ubound(self%unit_is_available,1)) then
            self%unit_is_available(unit_number) = .true.
        endif
        !$omp end critical (ProgramInstance_UnitManagement)
    end subroutine ReleaseUnit

    !>  \brief  Set the start date
    subroutine SetStartDate(self)
        use DatesAndTimes
        class(ProgramInstance), intent(inout)   ::  self
        type(DateAndTime)   ::  my_date_and_time
        self%date_start = my_date_and_time%DateAndTimeAsIntegers()
    end subroutine SetStartDate

    !>  \brief  Initialise a program instance
    subroutine Init(self, wanted_program_name, wanted_version, wanted_copyright_year)
        use DatesAndTimes
#ifdef __INTEL_COMPILER
        use ifport         ! if we are using the fortran compiler include IFPORT to check if we are reading from tty.
#else
#if defined(NAGFOR_COMPILER)
        use f90_unix_env
#endif
        !
#endif
        ! arguments
        class(ProgramInstance), intent(inout)  ::  self
        character(len=*), intent(in)           ::  wanted_program_name
        character(len=*), intent(in)           ::  wanted_version
        character(len=*), intent(in)           ::  wanted_copyright_year


        ! private variables
        type(DateAndTime)                      ::  my_date_and_time
        character(len=2),    parameter         ::  text_width    =    '30'
        integer                                ::  counter
        integer                                ::  length_of_program_name
        integer                                ::  number_of_needed_spaces
        integer                                ::  length_of_version_string


#ifdef __INTEL_COMPILER
        if (this_image() .eq. 1) then
#endif

        allocate(character(len=len_trim(adjustl(wanted_program_name))) :: self%program_name)

        ! start work

        self%date_start = my_date_and_time%DateAndTimeAsIntegers()
        self%program_name = wanted_program_name

        !Set all unit slots to available
        self%unit_is_available = .true.

        ! Parse command line arguments

        write(*,*)
        call self%ParseCommandLineArguments()

        ! Check to see if we are connected to a terminal
#ifdef __INTEL_COMPILER
        self%reading_from_terminal = isatty(5)
#elif defined(NAGFOR_COMPILER)
        call isatty(5,self%reading_from_terminal)
#endif

        length_of_program_name = len(trim(adjustl(wanted_program_name)))
        length_of_program_name = length_of_program_name + 19
        number_of_needed_spaces = 30 - (length_of_program_name / 2)

        do counter = 0, number_of_needed_spaces
            write(*, "(a)", advance="no")' '
        enddo

        write(*,'(3a)') '**  Welcome to ', trim(adjustl(wanted_program_name)), '  **'

        length_of_version_string = len(trim(adjustl(wanted_version)))
        length_of_version_string = length_of_version_string + 10
        number_of_needed_spaces = 30 - (length_of_version_string / 2)

        do counter = 0, number_of_needed_spaces
            write(*, "(a)", advance="no")' '
        enddo


        write(*,'(2a)')'Version - ', trim(adjustl(wanted_version))
        write(*,'(a)')                         ' '

        if (self%running_interactively .and. self%reading_from_terminal) then
            write(*,'(a'//text_width//',1x,a)')        'Input Mode:', 'Interactive'
        else if (.not. self%reading_from_terminal) then
            write(*,'(a'//text_width//',1x,a)')        'Input Mode:', 'Batch'
        else
            write(*,'(a'//text_width//',1x,a)')        'Input Mode:', 'Control File'
            write(*,'(a'//text_width//',1x,a)')        'Control Filename:',    trim(adjustl(self%control_filename))
        endif

        write(*,'(a'//text_width//',1x,a)')        'Date & Time:', my_date_and_time%DateAndTimeAsString()
        write(*,'(a)')

        write(*,'(/3a)') 'Copyright ', wanted_copyright_year, ' Howard Hughes Medical Institute. All rights reserved.'
        write(*,'(a)') 'Use is subject to Janelia Farm Research Campus Software Copyright 1.1'
        write(*,'(a/)') 'license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )'

#ifdef __INTEL_COMPILER
        endif
        flush(6)
        sync all
#endif

    end subroutine Init

    !>  \brief  terminate the program
    subroutine Terminate(self)
        use DatesAndTimes

        ! arguments
        class(ProgramInstance), intent(inout)  ::  self
        ! private variables
        type(DateAndTime)   ::  my_date_and_time

        ! start work
#ifdef __INTEL_COMPILER
        flush(6)
        sync all
        if (this_image() .eq. 1) then
#endif
        self%date_finish = my_date_and_time%DateAndTimeAsIntegers()
        write(*,*)
        write(*,'(2a)') 'Total execution time : ', my_date_and_time%DurationAsString(  &
                                real(my_date_and_time%SecondsBetweenDates(self%date_start, self%date_finish)))
        write(*,'(4a)') my_date_and_time%DateAndTimeAsString(), '  : ', self%program_name, ' finished cleanly.'
        write(*,*)
#ifdef __INTEL_COMPILER
        endif
        sync all
#endif
        if (allocated(self%program_name)) deallocate(self%program_name)
        if (allocated(self%control_filename)) deallocate(self%control_filename)
    end subroutine Terminate



    subroutine TerminateWithFatalError(self,name_of_calling_function, wanted_error)
        use DatesAndTimes
#ifdef __INTEL_COMPILER
        use ifcore
#endif

        ! arguments
        class(ProgramInstance), intent(inout)   ::  self
        character(len=*),       intent(in)      ::  name_of_calling_function
        character(len=*),       intent(in)      ::  wanted_error

        ! private variables
        type(DateAndTime)                       ::  my_date_and_time

        ! start work
#if defined(__INTEL_COMPILER)
#if defined(WITH_DEBUG_SYMBOLS)
        call tracebackqq(user_exit_code=-1)
#endif
#endif
        self%date_finish = my_date_and_time%DateAndTimeAsIntegers()
        write(*,*)
        write(*,'(2a)') 'Total execution time : ', my_date_and_time%DurationAsString(  &
                                real(my_date_and_time%SecondsBetweenDates(self%date_start, self%date_finish)))
        write(*,'(/4a)',advance='no') my_date_and_time%DateAndTimeAsString(), ': Fatal error (', name_of_calling_function, '): '
        write(*,'(a//)') wanted_error

        ! Do some deallocation
        if (allocated(self%program_name)) deallocate(self%program_name)
        if (allocated(self%control_filename)) deallocate(self%control_filename)

#ifdef __INTEL_COMPILER
        error stop
#else
        stop
#endif

    end subroutine TerminateWithFatalError

    subroutine ParseCommandLineArguments(self)

        class(ProgramInstance), intent(inout)  ::  self

        ! local variables
        character(200)      ::  command_line    ! complete command line
        integer             ::  command_line_len    ! number of characters in command line
        integer             ::  arg_count   ! number of command-line arguments
        integer             ::  max_arg_count   ! maximum number of arguments
        integer             ::  stat    ! status variable for error-checking
        character(200)      ::  current_arg ! argument currently being parsed
        integer             ::  i
        logical             ::  file_exists

        max_arg_count   =   2

        self%control_filename = ''


        call get_command(command_line, command_line_len, stat)
        if (stat .eq. -1) then
            call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                              'command_line variable not long enough to hold the command line!')
        elseif (stat .gt. 0) then
            call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                              'command line could not be retrieved!')
        elseif (stat .lt. -1) then
            call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                              'fatal error in parse_command_line_arguments!')
        endif

        arg_count   =   command_argument_count()

        if (arg_count .gt. max_arg_count) then
            write(*,'(a,i1,a,i1,a,i1,a)')   '**warning (parse_command_line_arguments): ', arg_count,                            &
                                            ' command-line arguments were supplied, but only expecting ', max_arg_count,    &
                                            '. ', arg_count - max_arg_count, ' arguments will be ignored.'
        elseif (arg_count .eq. 0) then
            !write(*,'(a)') '**debug (parse_command_line_arguments): no arguments found, setting interactive mode'
            self%running_interactively = .true.
        endif

        ! if there is more than one argument, we need to parse them
        if (arg_count .gt. 0) then
            do i=1,arg_count

                call get_command_argument(i,current_arg,status=stat)
                if(stat .eq. -1) then
                    write(*,'(3a)')     '**error: argument supplied is longer', &
                                        ' than maximum allowable. either give a shorter argument, or recompile ',   &
                                        'program with larger character string.'
                    call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                      'fatal error in parse_command_line_arguments!')
                elseif(stat .ne. 0) then
                    call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                      'unknown fatal error in parse_command_line_arguments!')
                endif

                ! check if this is a flag
                if (current_arg(1:1) .eq. '-') then
                    if (current_arg(2:2) .eq. 'i') then
                        self%running_interactively = .true.
                    elseif (current_arg(1:18) .eq. '--old-school-input') then
                        ! ugly workaround to allow grandfathering of input for ctffind
                        self%running_interactively = .true.
                    else
                        write(*,'(2a)') '**error (parse_command_line_arguments): unrecognised flag: ', current_arg
                    endif
                else ! not a flag
                    inquire(file=trim(adjustl(current_arg)), exist=file_exists)
                    if (file_exists) then
                        self%control_filename = current_arg
                    else
                        write(*,'(2a)') '**warning(parse_command_line_arguments): this file does not exist: ', &
                                        trim(adjustl(current_arg))
                        write(*,'(3a)')  '**warning(parse_command_line_arguments): will run ', self%program_name, &
                                            ' in interactive mode'
                        self%control_filename = ''
                        self%running_interactively = .true.
                    endif
                endif
            enddo
        endif

        if ( self%control_filename .eq. '' .and. .not. self%running_interactively ) then
           call self%TerminateWithFatalError('ProgramInstance::ParseCommandLineArguments', &
                                                      'Not in interactive mode, but no control filename either')
        endif

    end subroutine ParseCommandLineArguments

end module ProgramInstances
