!>  \brief  class & methods defining run-time parameters and their behaviour
module UserSuppliedParameters

    implicit none

    private

    type,   abstract,   public  ::  UserSuppliedParameter
        character(len=:), allocatable       ::  keyword                 !<  keyword used in command files to identify the parameter
        character(len=:), allocatable       ::  description             !<  description
        logical                             ::  set_by_user =   .false. !<  flag to indicate whether the value has been set by the user, or is default

    contains

        procedure                           ::  PrintInfo
    end type

    type, extends(UserSuppliedParameter), public :: UserSuppliedReal
        real                           ::  value
    end type

    type, extends(UserSuppliedParameter), public :: UserSuppliedInteger
        integer                        ::  value
     end type

    type, extends(UserSuppliedParameter), public :: UserSuppliedLogical
        logical                        ::  value
    end type

    type, extends(UserSuppliedParameter), public :: UserSuppliedFilename
        character(len=:), allocatable  ::  value
    end type


    contains

    !>  \brief  print info about the parameter, over multiple lines if necessary
    subroutine PrintInfo(self)

        ! arguments
        class(UserSuppliedParameter),     intent(in)  ::  self
        ! private variables
        !integer                       ::  descr_length
        !integer                       ::  descr_lead_blanks
        !integer                       ::  num_lines
        !integer                       ::  current_line
        !integer                       ::  num_spaces
        !integer                       ::  begin
        !integer                       ::  end

        !character(len=5)              ::  num_spaces_here
        !character(len=:), allocatable ::  tmp_str



        integer,            parameter ::  keyword_column_width = 16
        character(len=4),   parameter ::  keyword_column_width_c   =   '16'
        integer,            parameter ::  value_column_width = 45
        character(len=4),   parameter ::  value_column_width_c   =   '45'
        integer,            parameter ::  setbyuser_column_width = 7
        character(len=4),   parameter ::  setbyuser_column_width_c   =   '7'
        integer,            parameter ::  description_column_width = 64
        character(len=4),   parameter ::  description_column_width_c   =   '64'


        !logical ::  no_space

       ! start work

        ! First off just write out the keyword..

        write(*,'(a'//keyword_column_width_c//' ,1x)',advance='no') trim(adjustl(self%keyword))

        select type (self)
!            type is (UserSuppliedParameter)
!            ! Do Nothing - this should never be called
            class is (UserSuppliedReal)
                write(*,'(f'//value_column_width_c//'.4  ,1x)',advance='no') self%value
            class is (UserSuppliedInteger)
                write(*,'(i'//value_column_width_c//'  ,1x)',advance='no') self%value
            class is (UserSuppliedLogical)
                write(*,'(l'//value_column_width_c//'  ,1x)',advance='no') self%value
            class is (UserSuppliedFilename)
                write(*,'(a'//value_column_width_c//'  ,1x)',advance='no') self%value
        end select

        write(*,'(l'//setbyuser_column_width_c//'  ,1x)',advance='no') self%set_by_user

!        descr_length = len_trim(adjustl(self%description))
!       descr_lead_blanks = len(self%description) - len_trim(self%description)
!       num_lines = ceiling(real(descr_length)/real(description_column_width))
!       num_spaces = keyword_column_width + 1 + value_column_width + 1 + setbyuser_column_width + 1
!       allocate(character(description_column_width) :: tmp_str)
!       do current_line=1,num_lines
!           if (current_line .gt. 1) then
!               ! we need a new line
!               write(*,'(a)') ' '
!           endif
!           !
!           begin = ((current_line-1)* description_column_width) + 1
!           end   =  current_line   * description_column_width
!           !
!           !begin = begin + descr_lead_blanks
!           !end   = end   + descr_lead_blanks
!           !
!           end   = min(end,descr_lead_blanks+descr_length)
!           tmp_str = self%description(begin:end)
!           no_space = .false.
!           if (current_line .gt. 1) then
!               write(num_spaces_here,'(i5)') num_spaces
!           else
!               num_spaces_here = '0'
!               no_space = .true.
!           endif
!           if (no_space) then
!               write(*,                       '(a'//description_column_width_c//',1x)',advance='no') adjustl(tmp_str)
!           else
!               write(*,'('//num_spaces_here//'x,a'//description_column_width_c//',1x)',advance='no') adjustl(tmp_str)
!           endif
!       enddo
        write(*,'(a)') ' '
    end subroutine PrintInfo
end module UserSuppliedParameters

!
!
!    type, extends(rtp), public  ::  rtp_real
!        !real    ::  default_value   !<  default value for the parameter
!        real            ::  value           !<  current value for the parameter
!        contains
!            procedure   ::  init                =>  rtp_init_real
!            procedure   ::  print_info          =>  rtp_real_print_info
!            procedure   ::  print_to_command_file   =>  rtp_real_print_to_command_file
!    end type
!
!    type, extends(rtp), public  ::  rtp_intg
!        !integer ::  default_value   !<  default value for the parameter
!        integer ::  value           !<  current value for the parameter
!        contains
!            procedure   ::  init                =>  rtp_init_intg
!            procedure   ::  print_info          =>  rtp_intg_print_info
!            procedure   ::  print_to_command_file   =>  rtp_intg_print_to_command_file
!    end type
!
!    type, extends(rtp), public  ::  rtp_logi
!        !logical ::  default_value   !<  default value for the parameter
!        logical ::  value           !<  current value for the parameter
!        contains
!            procedure   ::  init                =>  rtp_init_logi
!            procedure   ::  print_info          =>  rtp_logi_print_info
!            procedure   ::  print_to_command_file   =>  rtp_logi_print_to_command_file
!    end type
!
!    type, extends(rtp), public  ::  rtp_char
!        !character(len=line_len)    ::  default_value   !<  default value for the parameter
!        character(len=line_len)    ::  value           !<  current value for the parameter
!        contains
!            procedure   ::  init                =>  rtp_init_char
!            procedure   ::  print_info          =>  rtp_char_print_info
!            procedure   ::  print_to_command_file   =>  rtp_char_print_to_command_file
!    end type
!
!    interface   rtp_assign
!        module procedure    rtp_assign_real
!        module procedure    rtp_assign_char
!        module procedure    rtp_assign_intg
!        module procedure    rtp_assign_logi
!    end interface
!
!    interface rtp_init
!        module procedure    rtp_init_real
!        module procedure    rtp_init_intg
!        module procedure    rtp_init_char
!        module procedure    rtp_init_logi
!    end interface rtp_init
!
!    contains
!
!    !>  \brief  print rtp information in format suitable for command file
!    subroutine rtp_real_print_to_command_file(self,lun)
!        implicit none
!        ! arguments
!        class(rtp_real),    intent(in)  ::  self
!        integer,            intent(in)  ::  lun
!        ! private variables
!
!        ! start work
!        write(lun,'(a'//com_colwidth_c(1)//',f0.6,10x,a,1x,a)')   adjustl(self%keyword), self%value,    &
!                                                                '! ', trim(adjustl(self%description))
!    end subroutine rtp_real_print_to_command_file
!
!    !>  \brief  print rtp information in format suitable for command file
!    subroutine rtp_intg_print_to_command_file(self,lun)
!        implicit none
!        ! arguments
!        class(rtp_intg),    intent(in)  ::  self
!        integer,            intent(in)  ::  lun
!        ! private variables
!
!        ! start work
!        write(lun,'(a'//com_colwidth_c(1)//',i0,10x,a,1x,a)') adjustl(self%keyword), self%value, '! ', &
!                                                            trim(adjustl(self%description))
!    end subroutine rtp_intg_print_to_command_file
!
!    !>  \brief  print rtp information in format suitable for command file
!    subroutine rtp_char_print_to_command_file(self,lun)
!        implicit none
!        ! arguments
!        class(rtp_char),    intent(in)  ::  self
!        integer,            intent(in)  ::  lun
!        ! private variables
!
!        ! start work
!        write(lun,'(a'//com_colwidth_c(1)//',a,10x,a,1x,a)') adjustl(self%keyword), trim(adjustl(self%value)), '! ', &
!                                                            trim(adjustl(self%description))
!    end subroutine rtp_char_print_to_command_file
!
!    !>  \brief  print rtp information in format suitable for command file
!    subroutine rtp_logi_print_to_command_file(self,lun)
!        implicit none
!        ! arguments
!        class(rtp_logi),    intent(in)  ::  self
!        integer,            intent(in)  ::  lun
!        ! private variables
!
!        ! start work
!        write(lun,'(a'//com_colwidth_c(1)//',l1,10x,a,1x,a)') adjustl(self%keyword), self%value, '! ',   &
!                                                            trim(adjustl(self%description))
!    end subroutine rtp_logi_print_to_command_file
!
!    !>  \brief  print information about a real rtp
!    subroutine rtp_real_print_info(self)
!        implicit none
!        ! arguments
!        class(rtp_real), intent(in)  ::  self
!        ! private variables
!
!        ! start work
!        call rtp_print_info_slave_1(self)
!        write(*,'(f'//colwidth_c(2)//'.4,1x)',advance='no') self%value!, self%default_value
!        call rtp_print_info_slave_2(self)
!        call rtp_print_info_slave_3(self)
!        write(*,'(a)') ' '
!    end subroutine rtp_real_print_info
!
!    !>  \brief  print information about a real rtp
!    subroutine rtp_intg_print_info(self)
!        implicit none
!        ! arguments
!        class(rtp_intg), intent(in)  ::  self
!        ! private variables
!
!        ! start work
!        call rtp_print_info_slave_1(self)
!        write(*,'(i'//colwidth_c(2)//',1x)',advance='no') self%value!, self%default_value
!        call rtp_print_info_slave_2(self)
!        call rtp_print_info_slave_3(self)
!        write(*,'(a)') ' '
!    end subroutine rtp_intg_print_info
!
!    !>  \brief  print information about a real rtp
!    subroutine rtp_logi_print_info(self)
!        implicit none
!        ! arguments
!        class(rtp_logi), intent(in)  ::  self
!        ! private variables
!
!        ! start work
!        call rtp_print_info_slave_1(self)
!        write(*,'(l'//colwidth_c(2)//',1x)',advance='no') self%value!, self%default_value
!        call rtp_print_info_slave_2(self)
!        call rtp_print_info_slave_3(self)
!        write(*,'(a)') ' '
!    end subroutine rtp_logi_print_info
!
!    !>  \brief  print information about a real rtp
!    subroutine rtp_char_print_info(self)
!        implicit none
!        ! arguments
!        class(rtp_char), intent(in)  ::  self
!        ! private variables
!
!        ! start work
!        call rtp_print_info_slave_1(self)
!        write(*,'(a'//colwidth_c(2)//',1x)',advance='no') trim(self%value)
!        call rtp_print_info_slave_2(self)
!        call rtp_print_info_slave_3(self)
!        write(*,'(a)') ' '
!    end subroutine rtp_char_print_info
!
!    !>  \brief  print basic information (not the value) of an rtp
!    subroutine rtp_print_info_slave_1(self)
!        implicit none
!        ! arguments
!        class(rtp),     intent(in)  ::  self
!        ! private variables
!
!        ! start work
!
!        write(*,'(a'//colwidth_c(1)//',1x)',advance='no') trim(adjustl(self%keyword))
!    end subroutine rtp_print_info_slave_1
!
!    !>  \brief  print basic information (not the value) of an rtp
!    subroutine rtp_print_info_slave_2(self)
!        implicit none
!        ! arguments
!        class(rtp),     intent(in)  ::  self
!        ! private variables
!
!        ! start work
!
!        write(*,'(l'//colwidth_c(3)//',1x)',advance='no') self%set_by_user
!    end subroutine rtp_print_info_slave_2
!
!    !>  \brief  print the description of an rtp, over multiple lines if necessary
!    subroutine rtp_print_info_slave_3(self)
!        implicit none
!        ! arguments
!        class(rtp),     intent(in)  ::  self
!        ! private variables
!        integer ::  descr_length, descr_lead_blanks
!        integer :: num_lines, current_line, num_spaces
!        character(len=5)    ::  num_spaces_here
!        integer ::  begin,end
!        character(len=:), allocatable :: tmp_str
!        logical ::  no_space
!        ! start work
!
!!        write(*,'(a<colwidth(4)>,x)',advance='no') adjustl(self%description)
!
!        descr_length = len_trim(adjustl(self%description))
!        descr_lead_blanks = len(self%description) - len_trim(self%description)
!        num_lines = ceiling(real(descr_length)/real(colwidth(4)))
!        num_spaces = colwidth(1) + 1 + colwidth(2) + 1 + colwidth(3) + 1
!
!        allocate(character(len=colwidth(4)) :: tmp_str)
!
!
!        do current_line=1,num_lines
!            if (current_line .gt. 1) then
!                ! we need a new line
!                write(*,'(a)') ' '
!            endif
!            !
!            begin = (current_line-1)* colwidth(4) + 1
!            end   =  current_line   * colwidth(4)
!            !
!            !begin = begin + descr_lead_blanks
!            !end   = end   + descr_lead_blanks
!            !
!            end   = min(end,descr_lead_blanks+descr_length)
!            tmp_str = self%description(begin:end)
!            no_space = .false.
!            if (current_line .gt. 1) then
!                write(num_spaces_here,'(i5)') num_spaces
!            else
!                num_spaces_here = '0'
!                no_space = .true.
!            endif
!            if (no_space) then
!                write(*,                       '(a'//colwidth_c(4)//',1x)',advance='no') adjustl(tmp_str)
!            else
!                write(*,'('//num_spaces_here//'x,a'//colwidth_c(4)//',1x)',advance='no') adjustl(tmp_str)
!            endif
!        enddo
!    end subroutine rtp_print_info_slave_3
!
!    !>  \brief  print out header for rtp descriptions
!    subroutine rtp_print_info_header()
!        implicit none
!        ! arguments
!
!        ! private variables
!        character(len=colwidth(4))  ::  description
!        character(len=colwidth(5))  ::  help
!        ! start work
!        description = 'description'
!        help = 'help'
!        write(*,'(a'//colwidth_c(1)//',1x,a'//colwidth_c(2)//',1x,a'//colwidth_c(3)//',1x,a'//colwidth_c(4)//')') &
!                                                                                'key', 'value', 'set by user?', adjustl(description)
!    end subroutine rtp_print_info_header
!
!    !>  \brief  initialise a runtime parameter
!    subroutine rtp_init_real(self,default_value,keyword,description,help)
!        implicit none
!        ! arguments
!        class(rtp_real),                intent(inout)   ::  self
!        real,                           intent(in)      ::  default_value   !<  default value for the parameter
!        character(len=*),               intent(in)      ::  keyword         !<  keyword used in frealix command files to identify the parameter
!        character(len=*),   optional,   intent(in)      ::  description     !<  description used in the question and as comment in frealix command files
!        character(len=*),   optional,   intent(in)      ::  help            !<  help to be given to user to describe what the parameter means / does
!        ! start work
!        call rtp_init_slave(self,keyword,description=description,help=help)
!        !self%default_value  =   default_value
!        self%value          =   default_value
!        self%set_by_user    =   .false.
!    end subroutine rtp_init_real
!
!    !>  \brief  initialise a runtime parameter
!    subroutine rtp_init_intg(self,default_value,keyword,description,help)
!        implicit none
!        ! arguments
!        class(rtp_intg),                intent(inout)   ::  self
!        integer,                        intent(in)      ::  default_value   !<  default value for the parameter
!        character(len=*),               intent(in)      ::  keyword         !<  keyword used in frealix command files to identify the parameter
!        character(len=*),   optional,   intent(in)      ::  description     !<  description used in the question and as comment in frealix command files
!        character(len=*),   optional,   intent(in)      ::  help            !<  help to be given to user to describe what the parameter means / does
!        ! start work
!        call rtp_init_slave(self,keyword,description=description,help=help)
!        !self%default_value  =   default_value
!        self%value          =   default_value
!        self%set_by_user    =   .false.
!    end subroutine rtp_init_intg
!
!    !>  \brief  initialise a runtime parameter
!    subroutine rtp_init_logi(self,default_value,keyword,description,help)
!        implicit none
!        ! arguments
!        class(rtp_logi),                intent(inout)   ::  self
!        logical,                        intent(in)      ::  default_value   !<  default value for the parameter
!        character(len=*),               intent(in)      ::  keyword         !<  keyword used in frealix command files to identify the parameter
!        character(len=*),   optional,   intent(in)      ::  description     !<  description used in the question and as comment in frealix command files
!        character(len=*),   optional,   intent(in)      ::  help            !<  help to be given to user to describe what the parameter means / does
!        ! start work
!        call rtp_init_slave(self,keyword,description=description,help=help)
!        !self%default_value  =   default_value
!        self%value          =   default_value
!        self%set_by_user    =   .false.
!    end subroutine rtp_init_logi
!
!    !>  \brief  initialise a runtime parameter
!    subroutine rtp_init_char(self,default_value,keyword,description,help)
!        implicit none
!        ! arguments
!        class(rtp_char),                intent(inout)   ::  self
!        character(len=*),               intent(in)      ::  default_value   !<  default value for the parameter
!        character(len=*),               intent(in)      ::  keyword         !<  keyword used in frealix command files to identify the parameter
!        character(len=*),   optional,   intent(in)      ::  description     !<  description used in the question and as comment in frealix command files
!        character(len=*),   optional,   intent(in)      ::  help            !<  help to be given to user to describe what the parameter means / does
!        ! start work
!        call rtp_init_slave(self,keyword,description=description,help=help)
!        !self%default_value  =   default_value
!        self%value          =   default_value
!        self%set_by_user    =   .false.
!    end subroutine rtp_init_char
!
!
!
!    !>  \brief  initialise a runtime parameter
!    subroutine rtp_init_slave(self,keyword,description,help)
!        implicit none
!        ! arguments
!        class(rtp),                     intent(inout)   ::  self
!        character(len=*),               intent(in)      ::  keyword     !<  keyword used in frealix command files to identify the parameter
!        character(len=*),   optional,   intent(in)      ::  description !<  description used in the question and as comment in frealix command files
!        character(len=*),   optional,   intent(in)      ::  help        !<  help to be given to user to describe what the parameter means / does
!        ! private variables
!
!        ! start work
!        self%keyword    =   keyword
!        if (present(description)) then
!            self%description = description
!        endif
!
!    end subroutine rtp_init_slave
!
!    !>  \brief  the rtp's value pointer is pointed at the relevant variable
!    subroutine rtp_assign_real(self,var)
!        implicit none
!        ! arguments
!        class(rtp),         intent(inout)   ::  self
!        real,               intent(in)      ::  var
!        ! private variables
!        ! start work
!        select type(self)
!            type is (rtp_real)
!                self%value = var
!            class default
!                write(*,'(a)') '**error(rtp_assign_real): rtp / variable type mismatch'
!                call terminate('fatal error in rtp_assign_real')
!        end select
!    end subroutine rtp_assign_real
!
!    !>  \brief  the rtp's value pointer is pointed at the relevant variable
!    subroutine rtp_assign_intg(self,var)
!        implicit none
!        ! arguments
!        class(rtp),             intent(inout)   ::  self
!        integer,        target, intent(in)      ::  var
!        ! private variables
!        ! start work
!        select type(self)
!            type is (rtp_intg)
!                self%value = var
!            class default
!                write(*,'(a)') '**error(rtp_assign_intg): rtp / variable type mismatch'
!                call terminate('fatal error in rtp_assign_intg')
!        end select
!    end subroutine rtp_assign_intg
!
!    !>  \brief  the rtp's value pointer is pointed at the relevant variable
!    subroutine rtp_assign_logi(self,var)
!        implicit none
!        ! arguments
!        class(rtp),             intent(inout)   ::  self
!        logical,        target, intent(in)      ::  var
!        ! private variables
!        ! start work
!        select type(self)
!            type is (rtp_logi)
!                self%value = var
!            class default
!                write(*,'(a)') '**error(rtp_assign_logi): rtp / variable type mismatch'
!                call terminate('fatal error in rtp_assign_logi')
!        end select
!    end subroutine rtp_assign_logi
!
!    !>  \brief  the rtp's value pointer is pointed at the relevant variable
!    subroutine rtp_assign_char(self,var)
!        implicit none
!        ! arguments
!        class(rtp),                     intent(inout)   ::  self
!        character(len=*),       target, intent(in)      ::  var
!        ! private variables
!        ! start work
!        select type(self)
!            type is (rtp_char)
!                self%value = var
!            class default
!                write(*,'(a)') '**error(rtp_assign_char): rtp / variable type mismatch'
!                call terminate('fatal error in rtp_assign_char')
!        end select
!    end subroutine rtp_assign_char
!
!
!
!    !>  \brief  fetch the value for a runtime parameter, either via user interface or by parsing command file
!    subroutine rtp_get_value(self)
!        use ui
!        use parsers
!        implicit none
!        ! arguments
!        class(rtp),                 intent(inout)   ::  self
!        ! private variables
!        character(len=line_len) ::  question
!        ! start work
!        question = 'please give ' // trim(self%keyword)
!        if (interactive) then
!            ! eventually, the ui routines should be rewritten to take rtp objects as arguments, so that the select type will be done at a lower level
!            select type(self)
!                type is (rtp_char)
!                    call ui_get_char(question,self%value)
!                type is (rtp_logi)
!                    call ui_get_logi(question,self%value)
!                type is (rtp_real)
!                    call ui_get_real(question,self%value)
!                type is (rtp_intg)
!                    call ui_get_intg(question,self%value)
!            end select
!            self%set_by_user = .true.
!        else
!            select type(self)
!                type is (rtp_char)
!                    call file_get_char(control_filename,self%keyword,self%value,found=self%set_by_user)
!                type is (rtp_logi)
!                    call file_get_logi(control_filename,self%keyword,self%value,found=self%set_by_user)
!                type is (rtp_real)
!                    call file_get_real(control_filename,self%keyword,self%value,found=self%set_by_user)
!                type is (rtp_intg)
!                    call file_get_intg(control_filename,self%keyword,self%value,found=self%set_by_user)
!            end select
!        endif
!    end subroutine rtp_get_value
!
!    !>  \brief  given an array of rtps, and a string of characters, find the rtp value in the string if possible
!    subroutine  rtp_from_string(rtps,string)
!        use parsers
!        implicit none
!        ! arguments
!        class(rtp),         allocatable,    intent(inout)   ::  rtps(:) !<  array of runtime parameters to check the line against
!        character(len=word_len),            intent(in)      ::  string  !<  a line (string) of characters
!        ! private variables
!        character(len=word_len),    allocatable ::  words(:)
!        character(len=word_len),    allocatable ::  label, value
!        integer ::  i
!        ! start work
!
!        !
!        ! analyse the line. the first word is the label, the second is the value
!        !
!        words = split(string)
!        label = words(1)
!        value = words(2)
!
!        !
!        ! loop over the array of rtps, to find the one, if any, which is designated by the string
!        !
!        do i=1,size(rtps)
!            if (string_equal(label,rtps(i)%keyword)) then
!                call rtp_value_from_string(rtps(i),value)
!                rtps(i)%set_by_user = .true.
!            endif
!        enddo
!    end subroutine rtp_from_string
!
!    !>  \brief  set rtp's value to that described by the supplied string
!    subroutine rtp_value_from_string(self,string)
!        implicit none
!        ! arguments
!        class(rtp),                         intent(inout)   ::  self    !<  runtime parameter for which the value should be set
!        character(len=word_len),            intent(in)      ::  string  !<  character string with value for the rtp
!        ! private variables
!
!        ! start work
!        select type(self)
!            type is (rtp_char)
!                read(string,*) self%value
!            type is (rtp_logi)
!                read(string,*) self%value
!            type is (rtp_real)
!                read(string,*) self%value
!            type is (rtp_intg)
!                read(string,*) self%value
!        end select
!    end subroutine rtp_value_from_string

