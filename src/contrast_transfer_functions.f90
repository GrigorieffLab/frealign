!>\brief    Contrast transfer functions for electron micrographs
module ContrastTransferFunctions
    use Globals
    use Units
    implicit none
    private
    public :: ContrastTransferFunction

    type :: ContrastTransferFunction
        private
        real    ::  spherical_aberration                =   0.0             !<  Spherical aberration
        integer ::  spherical_aberration_units          =   millimeters
        real    ::  wavelength                          =   0.0             !<  Wavelength of the electrons
        integer ::  wavelength_units                    =   angstroms
        real    ::  amplitude_contrast                  =   0.0             !<  Fraction of amplitude contrast (from 0. to 1., usually btw .07 and .14, see Mindell 03)
        real    ::  defocus_1                           =   0.0             !< Underfocus along first axis  (positive values for underfocus; larger value = more underfocus)
        real    ::  defocus_2                           =   0.0             !< Underfocus along second axis
        integer ::  defocus_1_units                     =   microns
        integer ::  defocus_2_units                     =   microns
        real    ::  astigmatism_azimuth                 =   0.0             !< Azimuth of first axis. 0.0 means axis is at 3 o'clock
        integer ::  astigmatism_azimuth_units           =   degrees
        ! Fitting parameters
        real    ::  lowest_frequency_for_fitting        =   0.0
        integer ::  lowest_frequency_for_fitting_units  =   reciprocal_pixels
        real    ::  highest_frequency_for_fitting       =   0.5
        integer ::  highest_frequency_for_fitting_units =   reciprocal_pixels
        real    ::  astigmatism_tolerance               =   0.0                 !<  Expected (tolerated) astigmatism. During parameter search, astigmatism values much larger than this will be penalised against. Set to 0.0 to ignore this restraint.
        integer ::  astigmatism_tolerance_units         =   angstroms
        contains
            procedure,  public      ::  Init
            procedure,  public      ::  PrintInfo
            procedure,  public      ::  EvaluateAtSquaredSpatialFrequency
            procedure,  public      ::  CountNumberOfExtremaBeforeSquaredSpatialFrequency
            !
            procedure               ::  SetDefocusScalars
            procedure               ::  SetDefocusArray
            generic,    public      ::  SetDefocus => SetDefocusArray, SetDefocusScalars
            procedure,  public      ::  HasPixelUnits
            procedure,  public      ::  SetLowestFrequencyForFitting
            procedure,  public      ::  GetLowestFrequencyForFitting
            procedure,  public      ::  GetHighestFrequencyForFitting
            procedure,  public      ::  GetAstigmatismTolerance
            procedure,  public      ::  GetAstigmatism
            procedure,  public      ::  GetDefocusParametersInAngstromsAndDegrees
            procedure,  public      ::  SetDefocusParametersInAngstromsAndDegrees
            procedure,  public      ::  GetDefocusParameters
            procedure               ::  ConvertToPixelsAndRadians
            procedure,  public      ::  GetDefocus1InAngstroms
            procedure,  public      ::  GetDefocus2InAngstroms
            procedure,  public      ::  GetAstigmatismAzimuthInDegrees
            procedure,  public      ::  GetAstigmatismAzimuthInRadians
    end type

    contains

    !>  \brief  Initialise a CTF object
    subroutine Init(self,acceleration_voltage,spherical_aberration,amplitude_contrast, &
                         defocus_1,defocus_2,astigmatism_azimuth, &
                         lowest_frequency_for_fitting,highest_frequency_for_fitting, &
                         astigmatism_tolerance, &
                         pixel_size)
        ! arguments
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  acceleration_voltage            !<  Kilovolts
        real,                               intent(in)      ::  spherical_aberration            !<  mm
        real,                               intent(in)      ::  amplitude_contrast              !<  fraction
        real,                               intent(in)      ::  defocus_1                       !<  um
        real,                               intent(in)      ::  defocus_2                       !<  um
        real,                               intent(in)      ::  astigmatism_azimuth             !<  degrees
        real,                               intent(in)      ::  lowest_frequency_for_fitting    !<  1/A
        real,                               intent(in)      ::  highest_frequency_for_fitting   !<  1/A
        real,                               intent(in)      ::  astigmatism_tolerance           !<  A
        real,                               intent(in)      ::  pixel_size                      !<  A
        ! start work
        self%wavelength = akv_to_wl(acceleration_voltage)
        self%wavelength_units = angstroms
        self%spherical_aberration = spherical_aberration
        self%spherical_aberration_units = millimeters
        self%amplitude_contrast = amplitude_contrast
        self%defocus_1 = defocus_1
        self%defocus_2 = defocus_2
        self%defocus_1_units = microns
        self%defocus_2_units = microns
        self%astigmatism_azimuth = astigmatism_azimuth
        !if (present(lowest_frequency_for_fitting)) then
            self%lowest_frequency_for_fitting = lowest_frequency_for_fitting
            self%lowest_frequency_for_fitting_units = reciprocal_angstroms
        !endif
        !if (present(highest_frequency_for_fitting)) then
            self%highest_frequency_for_fitting = highest_frequency_for_fitting
            self%highest_frequency_for_fitting_units = reciprocal_angstroms
        !endif
        !if (present(astigmatism_tolerance)) then
            self%astigmatism_tolerance = astigmatism_tolerance
            self%astigmatism_tolerance_units = angstroms
        !endif
        !if (present(pixel_size)) then
            call self%ConvertToPixelsAndRadians(pixel_size)
        !endif
    end subroutine Init

    subroutine SetLowestFrequencyForFitting(self,lowest_frequency_for_fitting)
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  lowest_frequency_for_fitting    !<  1/pixels

        self%lowest_frequency_for_fitting = lowest_frequency_for_fitting
        self%lowest_frequency_for_fitting_units = reciprocal_pixels
    end subroutine SetLowestFrequencyForFitting

    !>  \brief  Set defocus and astigmatism
    subroutine SetDefocusScalars(self,defocus_1,defocus_2,astigmatism_azimuth)
        ! arguments
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  defocus_1               !<  pixels
        real,                               intent(in)      ::  defocus_2               !<  pixels
        real,                               intent(in)      ::  astigmatism_azimuth     !<  radians
        !
        self%defocus_1 = defocus_1
        self%defocus_1_units = pixels
        self%defocus_2 = defocus_2
        self%defocus_2_units = pixels
        self%astigmatism_azimuth = astigmatism_azimuth
        self%astigmatism_azimuth_units = radians
    end subroutine SetDefocusScalars
    subroutine SetDefocusArray(self,my_array)
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  my_array(:)
        call self%SetDefocusScalars(my_array(1),my_array(2),my_array(3))
    end subroutine SetDefocusArray
    subroutine SetDefocusParametersInAngstromsAndDegrees(self,defocus_parameters,pixel_size)
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  defocus_parameters(3)
        real,                               intent(in)      ::  pixel_size
        self%defocus_1 = defocus_parameters(1) / pixel_size
        self%defocus_1_units = pixels
        self%defocus_2 = defocus_parameters(2) / pixel_size
        self%defocus_2_units = pixels
        self%astigmatism_azimuth = convert(defocus_parameters(3),degrees,radians)
        self%astigmatism_azimuth_units = radians
    end subroutine

    function GetDefocusParameters(self) result(defocus)
        class(ContrastTransferFunction),    intent(in)      ::  self
        real                                                ::  defocus(3)
        defocus = [self%defocus_1,self%defocus_2,self%astigmatism_azimuth]
    end function
    function GetDefocusParametersInAngstromsAndDegrees(self,pixel_size) result (defocus)
        class(ContrastTransferFunction),    intent(in)      ::  self
        real                                                ::  defocus(3)
        real,                               intent(in)      ::  pixel_size
        defocus(1) = convert(self%defocus_1,self%defocus_1_units,angstroms,pixel_size)
        defocus(2) = convert(self%defocus_2,self%defocus_2_units,angstroms,pixel_size)
        defocus(3) = convert(self%astigmatism_azimuth,self%astigmatism_azimuth_units,degrees)
    end function
    real function GetDefocus1InAngstroms(self,pixel_size)
        class(ContrastTransferFunction),    intent(in)  ::  self
        real,   optional,                   intent(in)  ::  pixel_size
        GetDefocus1InAngstroms = convert(self%defocus_1,self%defocus_1_units,angstroms,pixel_size)
    end function
    real function GetDefocus2InAngstroms(self, pixel_size)
        class(ContrastTransferFunction),    intent(in)  ::  self
        real,   optional,                   intent(in)  ::  pixel_size
        GetDefocus2InAngstroms = convert(self%defocus_2,self%defocus_2_units,angstroms,pixel_size)
    end function
    real function GetAstigmatismAzimuthInDegrees(self) result(azimuth)
        class(ContrastTransferFunction),    intent(in)  ::  self
        azimuth = convert(self%astigmatism_azimuth,self%astigmatism_azimuth_units,degrees)
        azimuth = azimuth - 180.0e0 * nint(azimuth/180.0e0)
    end function
    real function GetAstigmatismAzimuthInRadians(self) result(azimuth)
        class(ContrastTransferFunction),    intent(in)  ::  self
        azimuth = convert(self%astigmatism_azimuth,self%astigmatism_azimuth_units,radians)
    end function

    subroutine PrintInfo(self)
        ! arguments
        class(ContrastTransferFunction),    intent(in)      ::  self
        ! start work
        write(*,'(a)')              '** ContrastTransferFunction **'
        write(*,'(a,f0.3,1x,a)')    'Wavelength = ', self%wavelength, unit_to_string(self%wavelength_units)
        write(*,'(a,f0.3,1x,a)')    'Spherical aberration = ', self%spherical_aberration, &
                                                                unit_to_string(self%spherical_aberration_units)
        write(*,'(a,f0.3)')         'Amplitude contrast = ', self%amplitude_contrast
        write(*,'(a,f0.3,1x,a)')    'Defocus 1 = ', self%defocus_1, unit_to_string(self%defocus_1_units)
        write(*,'(a,f0.3,1x,a)')    'Defocus 2 = ', self%defocus_2, unit_to_string(self%defocus_2_units)
        write(*,'(a,f0.3,1x,a)')    'Azimuth of astigmatism = ', self%astigmatism_azimuth, &
                                                                unit_to_string(self%astigmatism_azimuth_units)
        write(*,'(2(a,f0.3,a))')    'Frequencies for fitting: from ', self%lowest_frequency_for_fitting, &
                                    unit_to_string(self%lowest_frequency_for_fitting_units), ' to ', &
                                    self%highest_frequency_for_fitting, unit_to_string(self%highest_frequency_for_fitting_units)
        write(*,'(a,f0.1,1x,a)')    'Astigmatism tolerance = ', self%astigmatism_tolerance, &
                                                                unit_to_string(self%astigmatism_tolerance_units)
    end subroutine PrintInfo

    !>  \brief  Check that all units are in pixels / radians
    pure logical function HasPixelUnits(self)
        class(ContrastTransferFunction),    intent(in)      ::  self
        HasPixelUnits =         self%spherical_aberration_units == pixels &
                        .and.   self%wavelength_units == pixels &
                        .and.   self%defocus_1_units == pixels &
                        .and.   self%defocus_2_units == pixels &
                        .and.   self%astigmatism_azimuth_units == radians &
                        .and.   self%highest_frequency_for_fitting_units == reciprocal_pixels &
                        .and.   self%lowest_frequency_for_fitting_units == reciprocal_pixels &
                        .and.   self%astigmatism_tolerance_units == pixels
    end function HasPixelUnits

    pure real function GetLowestFrequencyForFitting(self) result(lowest_frequency)
        class(ContrastTransferFunction),    intent(in)  ::  self
        lowest_frequency = self%lowest_frequency_for_fitting
    end function

    pure real function GetHighestFrequencyForFitting(self) result(highest_frequency)
        class(ContrastTransferFunction),    intent(in)  ::  self
        highest_frequency = self%highest_frequency_for_fitting
    end function

    pure real function GetAstigmatismTolerance(self) result(astigmatism_tolerance)
        class(ContrastTransferFunction),    intent(in)  ::  self
        astigmatism_tolerance = self%astigmatism_tolerance
    end function

    pure real function GetAstigmatism(self) result(astigmatism)
        class(ContrastTransferFunction),    intent(in)  ::  self
        astigmatism = self%defocus_1 - self%defocus_2
    end function


    !>  \brief  Convert all values to pixel/radians units
    subroutine ConvertToPixelsAndRadians(self,pixel_size)
        class(ContrastTransferFunction),    intent(inout)   ::  self
        real,                               intent(in)      ::  pixel_size
        ! start work
        call unit_conversion(self%spherical_aberration,self%spherical_aberration_units,pixels,pixel_size)
        call unit_conversion(self%wavelength,self%wavelength_units,pixels,pixel_size)
        call unit_conversion(self%defocus_1,self%defocus_1_units,pixels,pixel_size)
        call unit_conversion(self%defocus_2,self%defocus_2_units,pixels,pixel_size)
        call unit_conversion(self%astigmatism_azimuth,self%astigmatism_azimuth_units,radians)
        call unit_conversion(self%lowest_frequency_for_fitting,self%lowest_frequency_for_fitting_units, &
                            reciprocal_pixels,pixel_size)
        call unit_conversion(self%highest_frequency_for_fitting,self%highest_frequency_for_fitting_units, &
                            reciprocal_pixels,pixel_size)
        call unit_conversion(self%astigmatism_tolerance,self%astigmatism_tolerance_units,pixels,pixel_size)
    end subroutine ConvertToPixelsAndRadians

    !>  \brief wrapper around eval_ctf which takes a ctf object as argument
    pure elemental function EvaluateAtSquaredSpatialFrequency(self,squared_spatial_frequency,azimuth, &
                                                                return_sign_only) result(ctf)
        ! arguments
        class(ContrastTransferFunction),intent(in)      ::  self
        real,                           intent(in)      ::  squared_spatial_frequency   !<  squared reciprocal pixels (0.25 is Nyquist)
        real,                           intent(in)      ::  azimuth                     !<  radians. azimuth at which to evalulate the ctf
        logical,        optional,       intent(in)      ::  return_sign_only            !<  return the sign of the ctf rather than its value
        !
        real(kind=4)                                    ::  ctf
        !
        ctf      = eval_ctf_slave(  self%spherical_aberration,self%wavelength,self%amplitude_contrast,    &
                                    self%defocus_1,self%defocus_2,self%astigmatism_azimuth, &
                                    squared_spatial_frequency,azimuth,return_sign_only)
    end function EvaluateAtSquaredSpatialFrequency

    !>  \brief  Count how many extrema the CTF goes through before reaching the given spatial frequency at the given azimuth
    integer function CountNumberOfExtremaBeforeSquaredSpatialFrequency(self,squared_spatial_frequency,azimuth) &
                                                                result(number_of_extrema)
        class(ContrastTransferFunction),intent(in)      ::  self
        real,                           intent(in)      ::  squared_spatial_frequency   !<  squared reciprocal pixels (0.25 is Nyquist)
        real,                           intent(in)      ::  azimuth                     !<  radians. azimuth at which to evalulate the ctf
        ! private variable
        integer,    parameter   ::  maximum_number_of_rings = 128  ! this should be increased to something very large when debugging has been done
        real,       allocatable ::  args_of_minima(:) ! CTF arguments at which minima (-1.0) occur
        real,       allocatable ::  args_of_maxima(:) ! CTF arguments at which maxima (+1.0) occur
        ! start work

        ! Allocate memory for temporary results
        allocate(args_of_minima(maximum_number_of_rings))
        allocate(args_of_maxima(maximum_number_of_rings))

        ! Solve the CTF equation
        call ctf_solve_for_arg(self%amplitude_contrast,args_of_minima,-1.0e0)
        call ctf_solve_for_arg(self%amplitude_contrast,args_of_maxima,+1.0e0)

        ! Convert CTF arguments to squared spatial frequencies
        args_of_minima = ctf_sq_sf_from_arg(self%spherical_aberration,self%wavelength, &
                                            self%defocus_1,self%defocus_2,self%astigmatism_azimuth,args_of_minima,azimuth)
        args_of_maxima = ctf_sq_sf_from_arg(self%spherical_aberration,self%wavelength, &
                                            self%defocus_1,self%defocus_2,self%astigmatism_azimuth,args_of_maxima,azimuth)

        ! Let's count
        number_of_extrema = count(args_of_minima .le. squared_spatial_frequency .and. args_of_minima .gt. 0.0e0)
        number_of_extrema = count(args_of_maxima .le. squared_spatial_frequency .and. args_of_maxima .gt. 0.0e0) + number_of_extrema

        ! Can't remember why, but the ctf_solve_for_arg returns pairs of solutions that are very close to each other, except for the very first minimum
        if (number_of_extrema .gt. 0) then
            number_of_extrema = (number_of_extrema-1)/2 + 1
        endif

    end function CountNumberOfExtremaBeforeSquaredSpatialFrequency

    !>  \brief returns the ctf, based on ctffind3 subroutine (see mindell 2003)
    pure elemental real function eval_ctf_slave(cs,wl,ampl_cont,dfmid1,dfmid2,angast,spa_freq_sq,ang,sign_only)
        real,                   intent(in)  ::  cs          !< spherical aberation (pixels)
        real,                   intent(in)  ::  wl          !< electron wavelength (pixels)
        real,                   intent(in)  ::  ampl_cont   !< fraction of amplitude contrast (from .07 to .14, see mindell 2003)
        real,                   intent(in)  ::  dfmid1      !< defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2      !< defocus along second axis (for astigmatic ctf, dfmid1 .ne. dfmid2) (pixels)
        real,                   intent(in)  ::  angast      !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  spa_freq_sq !< squared spatial frequency at which to compute the ctf (1/pixels^2)
        real,                   intent(in)  ::  ang         !< angle at which to compute the ctf (radians)
        logical,    optional,   intent(in)  ::  sign_only   !<  return 1.0 or -1.0
        ! private variables
        real    ::  wgh1    ! phase contrast
        real    ::  wgh2    ! amplitude contrast
        real    ::  ctf_arg
        ! start work

        ! compute phase and amplitude contrast
        wgh1    =   sqrt(1-ampl_cont**2)
        wgh2    =   ampl_cont

        ctf_arg = eval_ctf_arg(cs,wl,dfmid1,dfmid2,angast,spa_freq_sq,ang)
        eval_ctf_slave = - wgh1 * sin(ctf_arg) - wgh2 * cos(ctf_arg)

        if (present(sign_only)) then
            if (sign_only) eval_ctf_slave = sign(1.0,eval_ctf_slave)
        endif
    end function eval_ctf_slave



    !>  \brief returns the argument (radians) to the sine and cosine terms of the ctf
    pure elemental real function eval_ctf_arg(cs,wl,dfmid1,dfmid2,angast,spa_freq_sq,ang) result(ctf_arg)
        real,                   intent(in)  ::  cs          !< spherical aberation (pixels)
        real,                   intent(in)  ::  wl          !< electron wavelength (pixels)
        real,                   intent(in)  ::  dfmid1      !< defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2      !< defocus along second axis (for astigmatic ctf, dfmid1 .ne. dfmid2) (pixels)
        real,                   intent(in)  ::  angast      !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  spa_freq_sq !< square of spatial frequency at which to compute the ctf (1/pixels^2)
        real,                   intent(in)  ::  ang         !< angle at which to compute the ctf (radians)
        ! private variables
        real    ::  df      ! defocus at point at which we're evaluating the ctf
        ! start work
        ! compute the defocus
        df = eval_df(ang,dfmid1,dfmid2,angast)
        ! compute the ctf argument
        ctf_arg = pi * wl * spa_freq_sq * (df - 0.5 * wl**2 * spa_freq_sq * cs)
    end function eval_ctf_arg

    !>  \brief  Return the effective defocus given the astigmatism parameters and the azimuth of interest
    pure elemental real function eval_df(ang,dfmid1,dfmid2,angast) result (df)
        real,                   intent(in)  ::  ang         !< angle at which to compute the defocus (radians)
        real,                   intent(in)  ::  dfmid1      !< defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2      !< defocus along second axis (for astigmatic ctf, dfmid1 .ne. dfmid2) (pixels)
        real,                   intent(in)  ::  angast      !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        !
        df = 0.5*(dfmid1+dfmid2+cos(2.0*(ang-angast))*(dfmid1-dfmid2))
    end function eval_df

    !>  \brief  Given the argument to the ctf, return the spatial frequency
    pure elemental function ctf_sf_from_arg(cs,wl,dfmid1,dfmid2,angast,arg,ang) result(sf)
        real,                   intent(in)  ::  cs          !< Spherical aberation (pixels)
        real,                   intent(in)  ::  wl          !< Electron wavelength (pixels)
        real,                   intent(in)  ::  dfmid1      !< Defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2      !< Defocus along second axis (for astigmatic ctf, dfmid1 .NE. dfmid2) (pixels)
        real,                   intent(in)  ::  angast      !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  arg         !< Radius at which to compute the ctf (pixels)
        real,                   intent(in)  ::  ang         !< Angle at which to compute the ctf (radians)
        real                                ::  sf          !< Spatial frequency
        ! Private variables
        !real    ::  df      ! Defocus
        ! Start work
        sf = ctf_sq_sf_from_arg(cs,wl,dfmid1,dfmid2,angast,arg,ang)
        sf = sqrt(sf)
    end function ctf_sf_from_arg

    !>  \brief  Given the argument to the ctf, return the squared spatial frequency
    !!
    !!  According to Maxima, there are two positive solutions
    pure elemental function ctf_sq_sf_from_arg(cs,wl,dfmid1,dfmid2,angast,arg,ang) result(sq_sf)
        real,                   intent(in)  ::  cs          !< Spherical aberation (pixels)
        real,                   intent(in)  ::  wl          !< Electron wavelength (pixels)
        real,                   intent(in)  ::  dfmid1      !< Defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2      !< Defocus along second axis (for astigmatic ctf, dfmid1 .NE. dfmid2) (pixels)
        real,                   intent(in)  ::  angast      !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  arg         !< Radius at which to compute the ctf (pixels)
        real,                   intent(in)  ::  ang         !< Angle at which to compute the ctf (radians)
        real                                ::  sq_sf       !< Squared spatial frequency (pixels^-2)
        ! Private variables
        real    ::  df      ! Defocus
        ! Start work

        ! Compute the defocus
        df = eval_df(ang,dfmid1,dfmid2,angast)
        !
        if (pi*df**2 .gt. 2.0e0*wl*arg*cs) then
            sq_sf = sqrt(pi)*df*wl - sqrt(wl**2*(pi*df**2 - 2.0e0*wl*arg*cs) )
            sq_sf = sq_sf / (sqrt(pi)*wl**3*cs)
        else
            sq_sf = 0.0e0
        endif
    end function ctf_sq_sf_from_arg

    !>  \brief  Given the argument to the ctf, return the defocus
    pure subroutine ctf_df_from_arg(cs,wl,arg,sf,df)
        real,                   intent(in)  ::  cs          !< Spherical aberation (pixels)
        real,                   intent(in)  ::  wl          !< Electron wavelength (pixels)
        real,                   intent(in)  ::  arg         !< ctf argument (radians)
        real,                   intent(in)  ::  sf          !< Spatial frequency (1/pixels)
        real,                   intent(out) ::  df          !< Defocus (pixels)
        ! Private variables

        ! Start work

        ! Compute the defocus
        df = (0.5 * cs * pi * sf**4 * wl**3 + arg) / (pi * sf**2 * wl)
    end subroutine ctf_df_from_arg

    !>  \brief  Find the ctf argument values at which CTF=0
    !!
    !!  According to Wolfram Alpha, the solutions to a*cos(t) + b*sin(t) = c are:
    !!  1)  t = 2 (atan((b-sqrt(a^2+b^2-c^2))/(a+c))+pi*n),   with n an integer
    !!  2)  t = 2 (atan((b+sqrt(a^2+b^2-c^2))/(a+c))+pi*n),   with n an integer
    !!
    !!  In our case, a = - ampl_const
    !!           and b = - sqrt(1-AMPL_CONT**2)
    !!  and t is the "argument" to the ctf, which can be "converted" to a spatial frequency
    subroutine ctf_solve_for_arg(ampl_cont,sols,ctf_value)
        real,                   intent(in)      ::  ampl_cont
        real,   allocatable,    intent(inout)   ::  sols(:)     !<  Solutions
        real,   optional,       intent(in)      ::  ctf_value   !<  Value of the ctf at the solutions
        ! Private variables
        integer ::  i,j
        real    ::  a,b
        real    ::  c
        ! Start work
        if (.not. allocated(sols)) then
            call this_program%TerminateWithFatalError('ctf_solve_for_arg','array of solutions is not allocated')
        endif
        if (size(sols) .lt. 1) then
            write(*,'(a,i0)') '**error(ctf_solve_for_arg): array of solutions has unexpected size ', size(sols)
            call this_program%TerminateWithFatalError('ctf_solve_for_arg','array of solutions has unexpected size')
        endif
        c = 0.0e0
        if (present(ctf_value)) c = ctf_value
        !
        a = -ampl_cont
        b = -sqrt(1-ampl_cont**2)
        ! Loop over the zeroes
        do i=1,size(sols)/2
            j = (i-1)*2+1
            sols(j) = 2.0*(atan((b-sqrt(a**2+b**2-c**2))/(a+c))+pi*real(i-1))
            j = (i-1)*2+2
            sols(j) = 2.0*(atan((b+sqrt(a**2+b**2-c**2))/(a+c))+pi*real(i))
        enddo
    end subroutine ctf_solve_for_arg

    !>  \brief returns the argument (radians) to the sine and cosine terms of the ctf
    pure subroutine eval_ctf_arg_array(cs,wl,dfmid1,dfmid2,angast,spa_freq_sq,ang,ctf_arg)
        real,                   intent(in)  ::  cs              !< spherical aberation (pixels)
        real,                   intent(in)  ::  wl              !< electron wavelength (pixels)
        real,                   intent(in)  ::  dfmid1          !< defocus along first axis (pixels)
        real,                   intent(in)  ::  dfmid2          !< defocus along second axis (for astigmatic ctf, dfmid1 .ne. dfmid2) (pixels)
        real,                   intent(in)  ::  angast          !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,                   intent(in)  ::  spa_freq_sq(:)  !< square of spatial frequency at which to compute the ctf (1/pixels^2)
        real,                   intent(in)  ::  ang             !< angle at which to compute the ctf (radians)
        real,                   intent(out) ::  ctf_arg(:)
        ! private variables
        real    ::  df      ! defocus at point at which we're evaluating the ctf
        !real    ::  ccos
        ! start work
        ! compute the defocus
        df = eval_df(ang,dfmid1,dfmid2,angast)
        ! compute the ctf argument
        ctf_arg(:) = pi * wl * spa_freq_sq * (df - 0.5 * wl**2 * spa_freq_sq(:) * cs)
    end subroutine eval_ctf_arg_array

    !>  \brief  convert acceleration voltage in kv into electron wavelength in angstroms
    pure real function akv_to_wl(akv)
        real, intent(in) :: akv !<  kV
        akv_to_wl = 12.26/sqrt(1000.0*akv+0.9784*(1000.0*akv)**2/(10.0**6.0))
    end function akv_to_wl


end module
