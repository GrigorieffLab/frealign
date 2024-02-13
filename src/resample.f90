program resample

    use Globals


    implicit none

    call main()

    contains

    ! The actual subroutine
    subroutine main()
        use UserSuppliedParameters
        use UserInputs
        use Images
        use ImageFiles
        use ProgressBars
        ! Variables associated with user input
        type(UserInput)             ::  my_user_input
        type(UserSuppliedFilename)  ::  input_filename
        type(UserSuppliedFilename)  ::  output_filename
        type(UserSuppliedLogical)   ::  volume_not_stack
        type(UserSuppliedLogical)   ::  real_space_binning
        type(UserSuppliedInteger)   ::  real_space_binning_factor
        type(UserSuppliedInteger)   ::  new_x_dimension
        type(UserSuppliedInteger)   ::  new_y_dimension
        type(UserSuppliedInteger)   ::  new_z_dimension
        ! Other private variables
        type(ImageFile)             ::  input_file
        type(ImageFile)             ::  output_file
        integer                     ::  output_file_dimensions(3)
        integer                     ::  output_image_dimensions(3)
        real                        ::  output_voxel_size(3)
        logical                     ::  output_voxel_is_cubic
        type(Image)                 ::  input_image
        type(Image)                 ::  output_image
        type(ProgressBar)           ::  my_progress_bar
        integer                     ::  current_loc, number_of_locations
        integer                     ::  locs_done
        ! Start user interaction

        ! Initialise this_program with the program name and version number
        call this_program%Init('Resample', '0.0.0','2014')

        call my_user_input%Init(this_program%program_name)

        input_filename                  =   my_user_input%GetFilenameFromUser(  'Input filename',   &
                                                                    'Input image filename. Can be stack of images or a volume.', &
                                                                    'input_filename','input.mrc',file_must_exist=.true.)

        ! Open the input file
        call input_file%Init(input_filename%value)

        output_filename                 =   my_user_input%GetFilenameFromUser(  'Output filename',  &
                                                                    'Output image filename', &
                                                                    'output_filename','output.mrc',file_must_exist=.false.)
        volume_not_stack                =   my_user_input%GetLogicalFromUser(   'Is the input image a volume?', &
                                                                    'If you answer NO here, we will assume the input'// &
                                                                    ' is a 2D image or a stack of 2D images', &
                                                                    'volume_not_stack','NO')
        real_space_binning              =   my_user_input%GetLogicalFromUser(   'Real space binning?', &
                                                                    'Average neighbouring voxels in real space rather than'// &
                                                                    ' clipping or padding the Fourier transform', &
                                                                    'real_space_binning','NO')
        if (real_space_binning%value) then
            real_space_binning_factor   =   my_user_input%GetIntegerFromUser(   'Binning factor', &
                                                                    'Number of neighbouring input voxels to average'//&
                                                                    ' together in each dimension', &
                                                                    'binning_factor','2')
        else

            if (volume_not_stack%value) then
                write(*,'(a,3(i0,1x))') 'Input file dimensions are: ', input_file%GetDimensions()
            else
                write(*,'(a,2(i0,1x))') 'Input file dimensions are: ', input_file%GetDimension(1), input_file%GetDimension(2)
            endif

            new_x_dimension             =   my_user_input%GetIntegerFromUser(   'New X dimension', &
                                                                    'X dimension of output image(s)', 'new_x_dimension','100')
            new_y_dimension             =   my_user_input%GetIntegerFromUser(   'New Y dimension', &
                                                                    'Y dimension of output image(s)', 'new_y_dimension','100')
            if (volume_not_stack%value) then
                new_z_dimension         =   my_user_input%GetIntegerFromUser(   'New Z dimension', &
                                                                    'Z dimension of output image(s)', 'new_z_dimension','100')
            endif
        endif

        call my_user_input%UpdateDefaults()

        ! End of user interaction


        ! Work out the output file & image dimensions
        if (real_space_binning%value) then
            output_image_dimensions = input_file%GetDimensions() / real(real_space_binning_factor%value)
        else
            output_image_dimensions = (/new_x_dimension%value,new_y_dimension%value,new_z_dimension%value/)
        endif
        if ( volume_not_stack%value) then
            output_file_dimensions = output_image_dimensions
        else
            output_file_dimensions(1:2) = output_image_dimensions(1:2)
            output_file_dimensions(3) = input_file%GetStackSize()
            output_image_dimensions(3) = 1
        endif

        ! Work out output pixel size
        output_voxel_size = input_file%GetPixelSize() * real(input_file%GetDimensions()) / real(output_file_dimensions)

        ! Is the output voxel cubic?
        output_voxel_is_cubic = abs(output_voxel_size(1) - output_voxel_size(2)) .lt. 0.0001
        if (volume_not_stack%value) then
            output_voxel_is_cubic = output_voxel_is_cubic .and. abs(output_voxel_size(1) - output_voxel_size(3)) .lt. 0.0001
        endif
        if (.not. output_voxel_is_cubic) then
            call this_program%TerminateWithFatalError('Resample','Output pixel/voxel would not be square/cubic')
        endif


        ! Open the output file
        call output_file%Init(output_filename%value,dim_x=output_file_dimensions(1), &
                                                    dim_y=output_file_dimensions(2), &
                                                    dim_z=output_file_dimensions(3), &
                                                    pixel_size=output_voxel_size(1),&
                                                    delete_if_already_exists=.true.)

        if (volume_not_stack%value) then
            call input_image%ReadFromDisk(input_file,read_volume_not_slice=.true.)
            if (real_space_binning%value) then
                call this_program%TerminateWithFatalError('Resample','Real-space binning not yet implemented')
            else
                call input_image%ForwardFFT()
                call output_image%Allocate(dims=output_image_dimensions,in_real_space=.false.)
                output_image = (0.0,0.0)
                call input_image%ClipInto(output_image)
                call output_image%ComputeAmplitudeSpectrum(input_image)
                call output_image%BackwardFFT()
                call output_image%WriteToDisk(output_file)
            endif
        else
            if (real_space_binning%value) then
                call this_program%TerminateWithFatalError('Resample','Real-space binning not yet implemented')
            endif
            number_of_locations = input_file%GetStackSize()
            call my_progress_bar%Begin(number_of_locations)
            locs_done = 0
            !$omp parallel default(shared) private(current_loc,input_image,output_image)
            call input_image%Reset()
            call output_image%Reset()
            !$omp do
            do current_loc=1,number_of_locations
                !$omp critical (crit_read)
                call input_image%ReadFromDisk(input_file,current_loc)
                !$omp end critical (crit_read)
                if (.not. output_image%IsAllocated()) call output_image%Allocate(dims=output_image_dimensions)
                call input_image%ForwardFFT()
                call input_image%ClipInto(output_image)
                call output_image%BackwardFFT()
                !$omp critical (crit_write)
                call output_image%WriteToDisk(output_file,current_loc)
                locs_done = locs_done + 1
                !$omp end critical (crit_write)
                call my_progress_bar%Update(locs_done)
            enddo
            !$omp end do
            call input_image%Deallocate()
            call output_image%Deallocate()
            !$omp end parallel
            call my_progress_bar%Finish()
        endif

        call output_file%Close()
        call input_file%Close()


        call this_program%Terminate()

    end subroutine main

end program resample
