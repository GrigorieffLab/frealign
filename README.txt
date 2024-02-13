New features in v9.11

- Asymmetric reconstructions (enable by using lowercase
  symmetry symbols).

New features in v9.10:

- Faster FFT algorithm.

New features in v9.09:

- The log likelihood values (LogP) in the parameter files 
  are now not negated anymore.

New features in v9.08:

- New scripts to run Frealign. Run INSTALL for instructions
  on how to install scripts. After installation, type

  frealign_help

  for information on available run scripts.

- New flag FDUMP (replacing old flag FSTAT): If set to T,
  Frealign will dump intermediate files from a 3D
  reconstruction into a file for later merging with new
  program merge_3d. This allows parallelization of 3D
  reconstruction over many CPUs with speedup factors of
  10 or more. See example script included with the
  distribution.

New features in v9.07:

- New Card 7 parameter RCLAS: high-resolution limit used
  for classification. It should typically be set to the
  same resolution limit used also for the refinement, or a
  bit lower.

New features in v9.06:

- Memory management: The old flag IBLOW has been replaced
  by a new flag IMEM. For refinement in Mode 1, 2, 3 or 4
  several instances of Frealign are usually run to refine
  several parts of a stack of particles in parallel. For
  the refinement, the non-mp version of Frealign should be
  used (no _nm at the end of the executable). In previous
  versions, setting IBLOW = 4 accelerated the processing
  but required additional memory. To get the same speedup
  and memory usage as before, set IMEM = 1 or IMEM = 3
  (no difference if no reconstruction is calculated).

  Reconstruction can also use parallelization (usually done
  with Mode 0). The script running Frealign should include a
  line like this:

  setenv NCPUS 8

  to set the number of parallel CPUs to be used. There is
  now a choice of how Frealign does the parallelization. If
  IMEM = 0 (or 1  but this makes no difference if no 
  refinement is done), it parallelizes the way it used to.
  For small particles, it is recommended to use 4 or 8 CPUs 
  as additional CPUs will probably not make it run faster.
  If IMEM = 2 (or 3, again no difference unless refinement
  is done) Frealign uses a new parallelization scheme that
  should run faster but requires more memory. For this
  scheme, 8 or 16 CPUs are recommended provided there is 
  sufficient memory available (try it and monitor usage).

- Maximum likelihood multi-reference refinement

- Helical symmetry operator HP: allows the processing of
  helical filaments with a seam, such as microtubules (the
  additional 'P' stands for 'pseudo').

- Weighted correlation coefficient now uses SSNR table
  from the previous iteration for more accuracy.

- Option for nearest-neighbor and trilinear interpolation
  used in 3D reconstruction (INTERP = 0 or 1).

- IMPORTANT: Shifts are now measured in Angstroms (see below).

- Each particle is now evaluated in terms of a score. This
  replaces the previous "phase residual". This means that a
  high number now indicates a good fit with the reference
  (previously, a high number indicated a bad fit). This
  should be considered when setting thresholds to exclude
  particles from a reconstruction.

- The user has to enter an approximate molecular mass of the
  particle (parameter MW in CARD 2) to calculate the 
  optimal filter for the 3D reconstruction (option FFILT).


New features in v8.10:

- Implemented single particle Wiener filter (use old FCREF
  flag to activate; FCREF now called FFILT).

- Implemented B-factor correction of the final output map.
  The B-factor is determined from a line fit to the log of
  the map power spectrum. Use new FBFACT to activate.


New features in v8.09:

- Added helical symmetry operator (H): If used, an
  additional input line right after CARD 5 is required 
  with 5 numbers (see helical refinement example included
  with the distribution):

  ALPHA, RISE, #SUBUNITS, #STARTS, STIFFNESS

  ALPHA        The rotation angle in deg around the helical
               axis between subunits.
               Example for TMV test data: 22.03

  RISE         The shift in Angstrom along the helical axis
               from one subunit to the next.
               Example for TMV test data: 1.408

  #SUBUNITS    The number of unique subunits per helical
               segment. This will depend on how the helical
               filaments were segmented. If segments are
               taken with a step size D, then

                  #SUBUNITS = #STARTS x D / RISE

               Example for TMV test data: 49

  #STARTS      The number of starts in the helix.
               Example for TMV test data: 1

  STIFFNESS    A parameter that determines how strongly the
               PSI and THETA angles of each helical segment
               are restrained to the average PSI and THETA
               angles of the helical filament. Reasonable
               values are between 1 (weak restraint) and 100
               (very strong restraint).
               Example for TMV test data: 20.0

  The included program set_polarity can be used to check that
  the PSI angles of each segment point all segments into a
  consistent direction to ensure uniform polarity along a
  filament. A tolerance for PSI deviations from the average
  PSI angle can be provided. set_polarity reads a parameter
  file and adjusts the PSI angle of segments that fall outside
  this tolerance. Also, filament polarity can be flipped
  depending on a phase residual (PRES) threshold.

New features in v8.08:

- Now requires Fortran 90 compiler (gfortran, pgf90)
- Dynamically allocated memory: Parameters for NN1 and
  NNBIG are now set automatically
- New order of control flags in Card 1
- New flags: FSTAT (T or F) to switch on/off additional
  statistics (requires more memory); IBLOW (1,2 or 4) to
  specify size of padded reference structure (IBLOW=4 
  requires the most memory but results in the fastest
  search & refinement

New features in v8.07:

- Removed FLIP option (Card 1)
- Added flag FPART for defocus refinement of individual
  particles (Card 1)
- Added restraining function for defocus refinement (active
  when FPART=T): Set expected defocus uncertainty with
  parameter DFSIG (Card 7)


New features in v8.06:

- Memory saving (IFSC flag on CARD1)
- Multiprocessor support (requires Portlan Group compiler)
  >>> Need to compiler with _mp Makefiles and set NCPUS
      environmental variable to the number of cores to be used
- Output of the two maps containing half the data each
- New I/O subroutines
- IMAGIC image format implemented
- x,y shift parameter distrubution function used for restraints
- Optional use of FFTW


FREALIGN notes for V8.00      Niko        1.2.08
FREALIGN notes for V7.00      Niko        1.2.06
FREALIGN notes for V6.00      Niko        9.3.02
FREALIGN notes for V5.00	RH	  2.9.01
FREALIGN notes for V4.06	RH	 27.8.01
FREALIGN notes for V4.02	RH	  9.8.01
FREALIGN notes for V4.00	RH	 22.4.01
FREALIGN notes for V3.05	RH	 25.2.01
FREALIGN notes for V3.04	RH	 29.6.00
FREALIGN notes for V3.03	RH	 10.5.00
FREALIGN notes for V3.02	RH	 14.3.00
FREALIGN notes for V3.01	RH	 27.8.99
FREALIGN notes for V3.00	RH	 27.7.99
FREALIGN notes for V2.07        RH        4.7.99
FREALIGN notes for V2.06	RH	  1.6.99
FREALIGN notes for V2.05	RH	 21.2.99
FREALIGN notes for V2.03	RH	  2.2.99
FREALIGN notes for V2.02	RH	 24.1.99
FREALIGN notes for V2.01	RH	  4.1.99
FREALIGN notes for V2.00	RH	24.12.98
================================================

FREALIGN carries out search and refinement of particle parameters, CTF 
correction (allowing for astigmatism), and 3D reconstruction using 
interpolation in Fourier space. The magnification of each dataset and the 
defocus and astigmatism values of particles grouped by their film number can 
also be refined. Reconstructions can be corrected for Ewald sphere curvature.
On output, there are also diagnostic data, such as Fourier Shell Correlation,
Fourier Shell Phase Residual, Q-factors in resolution zones, Average Particle
Phase Residual between particles and 3D reconstruction, variance of the
reconstruction, and point spread function indicating resolution in 3 dimensions.

FREALIGN uses MRC, SPIDER or IMAGIC image file formats. The MRC
mapformat is identical to that used by CCP4 crystallography programs. Current 
array dimensions provide space for 256x256 pixel images (set parameter NN1 in
frealign_v6.f) and the program is limited at present to transforms of even
dimensions.

Some of the principles of FREALIGN are explained in:

1) N. Grigorieff (1998), J. Mol. Biol. 277, 1033-1046.
2) Stewart, A. & Grigorieff, N. (2004), Ultramicroscopy 102, 67-84.
3) Wolf, M., DeRosier, D. J. & Grigorieff, N. (2006), Ultramicroscopy 106, 376-382.
4) Sindelar, C. V. & Grigorieff, N. (2012), J. Struct. Biol. 180, 26-38.
5) Lyumkis, D. & Brilot, A. F., Theobald, D. L. & Grigorieff, N. (2013), J. Struct. Biol. 183, 377-388.

********************************************************************************
Input cards : cards 1 to 18, including either 10a or 10b are always needed.
===========   cards 5a and 18a are used rarely only if appropriate flag is set 

Card 1 describes the overall program flow control logic cards
Cards 2 to 5 are global parameters governing all the data for the run 
Cards 6 to 12 form a repeatable group, one set describing each dataset which 
      may be made up of any number of particle images from different films
Cards 13 and 14 are filenames for the 3D structure I/O (overwrites input)
Cards 15 to 18 are filenames for 3D diagnostic output files (18a is backup)

Card 1          CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FDUMP,IMEM,INTERP
Card 2		RO,RI,PSIZE,MW,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
Card 3		PMASK,DMASK - [0/1] parameters to include in refinement [1 1 1 1 1]
Card 4		IFIRST,ILAST - FIRST, LAST PARTICLES
Card 5		ASYM - symmetry required Cn,Dn,T,O,I,I1,I2,N or H (can be zero)
   Card 5a	  only if N and N.ne.0, ((SYMOP(J,K,I), J=1,3), K=1,3), I=1,N
   Card 5b	  only if H, ALPHA, RISE, #SUBUNITS, #STARTS, STIFFNESS
		  If H is followed by P ("HP"), Frealign will not reset the shift
		  parameters to be within one helical asymetric unit. This can be
                  useful when processing pseudo-helical filaments such as
                  microtubules with a seam.

Cards 6 to 12 describe each dataset, terminating with RELMAG=0 (or negative)
T Card 6	RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
|      			Relative magnification to apply		(1.0)
|      			Densitometer step size in microns	(7.0)
|      			Target phase residual for search/refine	(15.0)
|      			Worst phase residual for inclusion	(90.0)
|      			CS, KV					(2.0, 120.0)
|      			Beam tilt in X, Y direction		(0.0, 0.0)
| Card 7	RREC,RMAX1,RMAX2,DFSTD,RBFACT - map resoln, refine low/high,
|			defocus uncertainty, B-factor
| Card 8	FINPAT1 - PARTICLE IMAGE STACK FILENAME
| Card 9	FINPAT2 - MATCHING PROJECTIONS STACK (for O/P if FMATCH=T)
| Card 10    - at least one of the following is required
|  Card 10a	FINPAR - INPUT PARAMETER FILE, required if mode key=0,1,2,3,4
|  Card 10b	NIN,ABSMAGPIN,IFILMIN,DFMID1IN,DFMID2IN,ANGASTIN,MORE if mode<0
|      		  number of particles, magnification, film number, defocus 
|                 parameters for this dataset.  If MORE='1', then more cards
|      		  of this type follow, until MORE=0 terminates the parameter 
|      		  data.  The information on these cards is used to create a 
|      		  new parameter file.
| Card 11	FOUTPAR - OUTPUT PARAMETER FILE
| Card 12	FOUTSH  - OUTPUT PARAMETER SHIFTS FILE, File name for file 
|      			containing the parameter shifts for this refinement 
|______     		cycle (should become smaller if refinement converges).

Card 13		F3D    - 3D MAP FILE FOR INPUT and OUTPUT, Input 3D reference 
			reconstruction, overwritten by output reconstruction!
Card 14		FWEIGH - 3D WEIGHTS FILE FOR INPUT and OUTPUT, Input 3D file
			containing the sum of weights (as defined in the JMB
			reference), will be overwritten with new file!
Card 15		MAP1   - 3D MAP FILE FOR OUTPUT, containing only odd-numbered
			particles (if IFSC = 0, 1), even-numbered particles
			(IFSC = 2) or all particles (IFSC = 3).
Card 16		MAP2   - 3D MAP FILE FOR OUTPUT, containing only even-numbered
			particles (if IFSC = 0), or no output (if IFSC = 1,2,3).
Card 17         FPHA   - 3D PHASE RESIDUAL FILE FOR OUTPUT, Output 3D file with
			average phase differences for each voxel (in Fourier
			space) between images and reference.
Card 18		FPOI   - 3D POINT SPREAD FUNCTION FOR OUTPUT, Output 3D file 
			with point spread function indicating anisotropies 
			in resolution (in real space).
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Parameters as read in on above cards 1 to 7 and 10b (other cards are filenames)
================================================================================

CFORM	M/S/I	Input/Output image file format: M=MRC, S=Spider, I-Imagic

IFLAG	0/1/2/3	Mode key.      -4:bootstrap parameter file, then IFLAG=4
      			       -3:bootstrap parameter file, then IFLAG=3
      			0:Reconstruction only parameters as read in
      			1:Refinement & Reconstruction
      			2:Random Search & Refinement
      			3:Simple search & Refine
      			4:Search,Refine,Randomise & extend to RREC
      	       -4 = create a parameter file from scratch, then IFLAG=4
      	       -3 = create a parameter file from scratch, then IFLAG=3
      		0 = use previously determined particle parameters 
      			and calculate a 3D reconstruction; 
      		1 = carry out refinement and reconstruction starting with 
      			previous roughly determined parameters
                        NOTE: If RELMAG in the termination line in card 6 
                        is set to a negative value instead of 0.0 no 3D
                        reconstruction is calculated (refinement only). This
                        is useful when processing smaller parts of the data
                        stack on a computer cluster. A reconstruction is then
                        done with another run and CFORM=0 (see example scripts).
      		2 = carry out refinement with randomly assigned particle
      			parameters (to check if current particle parameters 
      			correct; if correct then all other parameters should 
      			give worse phase residual)
      		3 = carry out systematic parameter search for initial 
      			assignment, with subsequent refinement.
      		4 = systematic search & refinement of particle orientation
                        parameters, with a randomised ITMAX loop to speed up 
                        convergence (only tries randomisation until phase
                        residual goes below TARGET).  Then, resolution is 
      			subsequently extended step-by-step out to RREC - if 
      			residual never goes below TARGET, earlier versions of 
      			the program set the residual to 180 degrees (it could 
      			then be used in subsequent refinement 
      			as flag with TARGET/THRESH = negative (essentially 
                        a flag for complete and permanent subsequent exclusion 
                        of that particle image)).  This option was most useful 
                        for adding data to an existing structure analysis as 
                        resolution is gradually improved, but it has now been 
      			removed so that the residual is always the real value.

FMAG	T/F	Magnification refinement - one number is calculated for each 
      			new film. The orientation parameters must already be
      			available in the input parameter file and should
      			preferably have been determined using a different 3D
      			reference model than that used for the magnification
			refinement.

FDEF	T/F	Defocus refinement, DF1 and DF2 coupled/uncoupled if FASTIG=F/T

FASTIG	T/F	Astigmatism refinement - if F, then DF1 and DF2 are coupled

FPART	T/F	Defocus refinement for individual particles if FPART=T, otherwise
			defocus change is constrained to be the same for all 
			particles in one image

IEWALD  0/1/2/-1/-2  Ewald correction:
                0 = No correction
                1 = Do correction, simple insertion method (see publication #3
                    above).
                2 = Do correction, reference-based method (see publication #3
                    above).
                -1 = 1 but with inverted handedness (see publication #3 above).
                -2 = 2 but with inverted handedness (see publication #3 above).

FBEAUT	T/F	Apply extra real space symmetry averaging and masking to 
      			beautify final map just prior to output.

FFILT   T/F     Apply single particle Wiener filter to final reconstruction
                (see publication #4 above).

FBFACT  T/F     Determine and apply B-factor to final reconstruction

FMATCH	T/F	Write out matching projections after the refinement (for 
      		diagnostic purposes): this will have spacegroup 0 in MRC header.
      		Card 9 (whether or not FMATCH is set) describes FINPAT2,
      		the matching projections image stack filename.

IFSC    0/1/2/3 Calculation of FSC tabel:
                0 = Internally calculate two reconstructions with odd and even
                    numbered particles and generate FSC table at the end of the
                     run.
                Options 1, 2 and 3 reduce memory usage:
                1 = Only calculate one reconstruction using odd particles.
                2 = Only calculate one reconstruction using even particles.
                3 = Only calculate one reconstruction using all particles.

FDUMP	T/F	If set to T, dumps intermediate files from a 3D reconstruction
                and then terminates run. This feature allows splitting up
                3D reconstruction into several parts (using IFIRST,ILAST). The
		dump files can be read in and merged for a full reconstruction
		using the new program merge_3d (included with the distribution).

IMEM	0/1/2/3 Memory usage: 0 = least memory, 3 = most memory.
                0 - no padding of reference during refinement, no multi-volume
                parallelization during reconstruction (least memory usage)
                1 - padding of reference during refinement, no multi-volume
                parallelization during reconstruction
                2 - no padding of reference during refinement, multi-volume
                parallelization during reconstruction
                3 - padding of reference during refinement, multi-volume
                parallelization during reconstruction (most memory usage)

INTERP	0/1	Interpolation scheme used for 3D reconstruction:
		0 = Nearest neighbor
		1 = Trilinear (more time-consuming)

          ---------------------------------------

RO - Outer radius of reconstruction in Angstroms from centre of particle: 108.0
      Enter the outer radius of the volume to be reconstructed. The program will
      also apply a mask with a cosine edge to the particle image before
      processing (done inside CTFAPPLY using HALFW=6 pixels for cosine bell).

RI - Inner radius of reconstruction in Angstroms from centre of particle: 0.0
      Enter the inner radius of the volume to be reconstructed. This is useful
      for reconstructions of viruses and other particles that might be hollow
      or have a disordered core.

PSIZE - Required pixel size [Angstrom]: 3.00
      Enter pixel size in Angstroms required for output map.  The input 
      particle images are then reinterpolated by using the densitometer step 
      size and the magnification of the micrograph to calculate the relative
      magnification.  RELMAG, normally set to 1.0, can be used to make further 
      manual adjustments.   A warning message is printed if the pixel size in 
      the header of the input 3D map used as the reference differs by more 
      than 1% from the required pixel size for the output map.  The input 3D
      map is not interpolated, but it is assumed the user knows what they are
      doing.  Similarly, a warning is printed if the input particle image 
      pixel size is less than 0.65 or more than 1.50 times PSIZE.

MW - Approximate molecular mass of the partcle, in kDa. This is used to
      calculate the optimal filter for the 3D reconstruction (option FFILT).

WGH - % Amplitude contrast (-1...1): 0.07
      Amplitude contrast for CTF correction.  A negative value will invert the
      sign of the CTF to flip image contrast.  Can use WGH = -1.0 to set
      CTF = 1 to fit images to a model without CTF correction.

XSTD - number of standard deviations above mean for masking of input low-pass 
       filtered 3D model - note this 3D masking does not use RI 
             - if positive, calculates mask with subroutine D3MASK, equiv to 
               solvent flattening with 5-pixel-cosine-bell smoothed mask 
               boundary.  The mask is then used to multiply the input 3D map, 
               which is then used for all parameter refinement and subsequent 
               calculations.
             - if negative, calculates mask with subroutine D2MASK resulting  
               in a sharp binary (0/1) mask boundary for which is used for 
               both parameter refinement and reconstruction, and to mask and 
               output the matching projections.  Each matching particle image 
               is also always masked with a cosine bell edged function of 
               radius RI (see above RI).

PBC - Phase residual / pseudo-B-factor conversion Constant: 1.0
      Automatic weighting is applied to each particle: a pseudo-temperature (B)
      factor is applied to each particle according to its relative phase 
      residual against the reference. The weight is calculated as
                W = exp (-DELTAP/PBC * R^2)
      with DELTAP = relative phase residual (actual phase residual minus BOFF), 
      PBC = conversion constant (5.0 in the example),
      and R^2 the squared resolution in Fourier units (R = 0.0 ... 0.5).
      A large value for PBC (e.g. 100.0) gives equal weighting to each particle

BOFF - average phase residual: 60.0, approximate average phase residual of 
      all particles, used in calculating weights for contributions of different 
      particles to 3D map (see Grigorieff, 1998).

DANG - angular step size for the angular search used in modes IFLAG=3,4
      There is a program default value calculated taking resolution into
      account, but if this input value is non-zero, the program value is
      overridden.

ITMAX - number of cycles of randomised search/refinement used in modes IFLAG=2,4
      There is a program default value (10 cycles), but if this input value is
      non-zero, the program value is overridden.

IPMAX - number of potential matches in a search that should be tested further in
      a subsequent local refinement.

          ---------------------------------------

PMASK - 0/1 mask to exclude parameters from refinement. There are five numbers to
      allow  psi, theta, phi, deltaX, deltaY to be included or excluded from the
      search or refinement. The use of 1 1 1 1 1 causes all 5 parameters to be 
      refined for each particle.  The use of 0 0 0 1 1 would cause only the 
      position of each particle to be optimised.

DMASK - X,Y,Z and radius to discribe an area in the map that should be used for
      classification ("focussed classification"). Not extensively tested - use
      with caution.

          ---------------------------------------

IFIRST,ILAST - First and last particle to be included: 1,5000
      Number of first and last particle to be included in the complete 
      reconstruction, using data from all datasets.

          ---------------------------------------

ASYM - symmetry required Cn,Dn,T,O,I,I1,I2 or N (can be zero)
      		  n  = rotational symmetry required in pointgroup Cn or Dn
      		  N  = number of symmetry matrices to read in.
      		  T  = tetrahedral pointgroup 23
      		  O  = octahedral pointgroup 432
      		  I  = icosahedral 532 symmetry in setting 1 (5fold is on X) 
      		  I1 = also in setting 1 (X) - as used by Imagic
      		  I2 = in setting 2 (Y) - as used by Crowther et al

      For any standard pointgroup symmetries and setting, the symmetry 
        requested is used to calculate the symmetry operators to be applied 
        to the reconstruction at the end of the calculation and to the search
        angles in modes 3 and 4.
        Alternatively the symmetry matrices can be read in explicitly to allow 
        for unusual symmetries or unusual orientations. If the number of 
        symmetry operators, N, is not 0 then the symmetry matrices must follow 
        immediately, e.g.
      	0.5			0.8660254037844	0.0
      	-0.8660254037844	0.5		0.0
      	0.0			0.0		1.0
  	(this matrix carries out a 60deg rotation about the z axis)
        In version 2.01 and earlier, the possibility of applying Cn,Dn,T,O, or I 
        symmetry in standard settings is not available.

          ---------------------------------------

RELMAG,DSTEP,TARGET,THRESH - Relative magnification to apply to data set (0=END) 
      			densitometer step size in microns,
      			target phase residual for parameter refinement,
      			and phase residual threshold:     1.00, 7.0, 15.0, 90.0

      RELMAG: Relative magnification of data set.  The magnification feature 
      	allows for manual variations of magnification between data sets. The
      	micrograph magnification should be in the input parameter file or read 
      	in on card 8b.  Enter 0.0 or a negative value if no further data sets
      	are to be added.   If RELMAG is a negative number (e.g. -100.0), there 
      	is no map output, application of symmetry or statistics.

      DSTEP: Densitometer step size.  Each dataset defined by cards 5 to 10 must
      	have the same values for DSTEP as well as other parameters on this card.

      TARGET: Target phase residual (for resolution between RMAX1 and RMAX2) 
      	during parameter saerch and refinement, below which the search and/or 
      	refinement is terminated.  This is normally set low (e.g. THRESH=10.0)
      	This will give excellent determination of particle orientations.

      THRESH: Phase residual cut-off. Any particles with a higher overall phase
        residual will not be included in the reconstruction when IFLAG=0,1,2,3. 
      	This variable is often used with IFLAG=0 in separate runs to calculate 
      	maps using various values of THRESH to find the optimum value to 
      	produce the best map as judged from the statistics.

CS,AKV - CS [mm], V [kV]: 2.60,200.0
      Microscope parameters for this data set: Cs and voltage.

TX,TY - Beam tilt [mrad] in X, Y direction: 0.0, 0.0

          ---------------------------------------

RREC - Resol. of reconstruction in Angstroms, e.g. 10.0 
      Resolution to which the reconstruction is calculated.
      If several datasets have different values, the data is individually
      limited in the summation to the RREC resolution but symmetry is 
      applied, statistics output and the final map calculated to the 
      maximum resolution requested for any dataset.

RMAX1,RMAX2 - Resol. in refinement in Angstroms, low & high: 200.0,25.0
      Resolution of the data included in the search/refinement.  These
      two parameters are very important.  The successful alignment of
      particles depends critically on the signal-to-noise ratio of the
      cross-correlation or phase residual calculation, and exclusion of
      weak data at high resolution or spurious, very strong artefacts at
      low resolution can make a big difference.  Success can be judged 
      by whether the X,Y coordinates of the particle centres are reasonable.

DFSIG - Defocus uncertainty in Angstroms, e.g. 200.0
      This will restrain the change in defocus when refining defocus values
      for individual particles.

RBFACT - B-factor to apply to particle image projections before orientation 
      determination or refinement.  This allows inclusion of high resolution 
      data but with a reduced (or increased if negative) weight.  WGH and 
      RBFAC can be manipulated in particle parameter refinement as if they 
      were low pass and high pass filters.  WGH and the CTF are used to 
      correct the density in the final map, whereas RBFAC is not.
      NOTE: This option should be used with great care as amplification of 
      noisy high-resolution terms can lead to over-fitting and artificially
      high values in the FSC curve (se publication #2 above). FREALIGN uses an
      automatic weighting scheme and RBFACT should normally be set to 0.0.

          ---------------------------------------

NIN		- number of particle images for this film
ABSMAGPIN	- real magnification for this film (e.g. 45,000 x)
IFILMIN		- film number for this film
DFMID1IN,DFMID2IN,ANGASTIN - defocus and astigmatism for this film
MORE		- use '1' for more cards describing more particle images in
      		  this stack of images, '0' to terminate.

********************************************************************************

File formats for input/output files
==================================

All output files containing 3D maps have spacegroup 1 and 80 bytes of symmetry
information in the header (containing " X,Y,Z   "), whereas the matching image
stack, has spacegroup 0 with no symmetry bytes.  Note that the MRC format 
also defines a mapmode depending on how many bytes of real or complex data is
written out per pixel.

Parameter FILE for input or output:    Twelve entries per particle:
 ##  PSI  THETA  PHI  SHX  SHY  MAG  FILM  DF1  DF2  ANGAST  PRESA
      Particle number
      Phi   - first Eulerian angle, rotates model clockwise about Z (alpha)
      Theta - rotates about new Y after rotation about Z (beta), to produce a 
      	      characteristic view of the particle
      Psi   - third rotation about new Z (gamma), rotates particle without 
      	      changing the viewpoint
      	      this ZYZ rotation is standard Eulerian [a.k.a. alpha, beta, gamma]
      Shift X, Shift Y (in Angstrom)
      Magnification - this is real magnification (e.g. 45,000 x).  If 
      	magnification has been refined, then the output parameter file will 
      	contain the appropriate magnification so that the new value of 
      	RELMAG = 1.0
      Film number - identifies which particles are from same film
      defocus1 in Angstrom, defocus2 in Angstrom, angle for astigmatism
      	The two defocus values (positive for underfocus) give the defocus in two
      	perpendicular directions. They are the same if there is no astigmatism.
      	The angle is measured between the X-axis and the direction of the first
      	defocus measurement (a positive angle means a rotation of the defocus
      	direction to the left).
      Phase residual: This is actually the inverse cosine of a weighted correlation
        coefficient calculated during refinement (see publication #2 above).

Statistical table 1:    Thirteen entries per radius bin in table comparing two halves of data
C  NO. RESOL  RING RAD   FSPR    FSC   QFACT QRAN+SIG    SSNR  RFACT  PFREE PNOISE NONZERO  TOTVOX
C   2  398.1    0.0078   0.00  0.997   0.971  0.030     47.51  0.073  26.22 108.44      17      17
C   3  199.0    0.0156   0.03  0.997   0.559  0.040     19.98  0.080  80.41  75.08      41      41
C   4  132.7    0.0234   0.04  0.998   0.538  0.045     17.86  0.061  77.51  65.14      89      89
C   5   99.5    0.0313   0.66  0.982   0.398  0.057     13.37  0.184  57.01 109.59     129     129
C   6   79.6    0.0391   0.30  0.965   0.426  0.063     13.67  0.290  55.10  93.69     225     225
C   7   66.3    0.0469   3.07  0.908   0.166  0.065      4.93  0.385 104.20  82.42     253     253
C   8   56.9    0.0547  11.57  0.882   0.209  0.069      4.70  0.337  86.66 103.94     393     393
C------------------------------------------------------------------------------------------------
C  17   24.9    0.1250  54.79  0.328   0.122  0.107      1.56  0.671  85.50  80.16    1701    1701
C  18   23.4    0.1328  54.18  0.345   0.100  0.109      1.02  0.687  87.59  87.67    2021    2021
C  19   22.1    0.1406  73.02  0.169   0.083  0.110      0.73  0.649  77.97  87.75    2181    2181
C  20   21.0    0.1484  76.84  0.136   0.088  0.114      0.69  0.661  93.42  83.76    2473    2473
C  21   19.9    0.1563  56.51  0.311   0.132  0.117      1.44  0.712  85.63  86.15    2697    2697
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Statistical table 2:    Eight entries comparing merged data with reference data
C  NO. RESOL  RING RAD   FSPR   FSC   RFACT NONZERO  TOTVOX
C   2  398.1    0.0078   0.14  0.998  0.094      17      17
C   3  199.0    0.0156   1.30  0.969  0.705      41      41
C   4  132.7    0.0234   4.55  0.854  0.575      89      89
C   5   99.5    0.0313   0.68  0.958  0.953     129     129
C   6   79.6    0.0391   1.60  0.937  0.667     225     225
C   7   66.3    0.0469  41.69  0.452  1.213     253     253
C   8   56.9    0.0547   5.74  0.793  1.103     393     393
C----------------------------------------------------------
C  17   24.9    0.1250  38.32  0.520  0.650    1701    1701
C  18   23.4    0.1328  50.69  0.391  0.649    2021    2021
C  19   22.1    0.1406  65.20  0.228  1.052    2181    2181
C  20   21.0    0.1484  58.91  0.306  0.954    2473    2473
C  21   19.9    0.1563  43.31  0.482  0.613    2697    2697
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Abbreviations:

FSPR:       Fourier shell phase residual
FSC:        Fourier shell correlation
APPR:       Average particle phase residual (against reference)
Av. Q:      Average Q-factor
Critical Q: Q-factor for pure noise (as defined in the JMB reference)

SOME IMPORTANT NOTES:
--------------------
1.    Shifts are measured in Angstrom to make it easier to work with binned
      data. This is an important change from earlier versions (<9.00) that
      measured shifts in pixels. Older parameter files should be compatible
      with later versions as Frealign attempts to recognize older parameter
      files and converts pixel shifts into Angstroms.

2.    In earlier (<5.00) versions of the program, the application of symmetry 
      in subroutine APPLYSYMC was carried out using an extra phase shift (PSHFT) 
      which moved the centre of the symmetry axes of the particles by half a 
      pixel to be at (NSAM/2 - 0.5) along all three axes (e.g. [63.5,63.5,63.5] 
      in a 128x128x128 map).  This meant that the output map had the exact centre 
      at the midpoint of 8 pixels.  While this had no effect on the internal 
      consistency of the program itself since each successive round of 
      refinement automatically located the (X,Y) coordinates of each particle 
      image correctly relative to the 3D model, it meant that it was hard to 
      interface the program to others.  Accordingly, the centre for application 
      of the symmetry was moved to be exactly centred on the pixel at NSAM/2+1 
      from version 5.00 onwards.  The program changes were : to set HALFP to 
      be zero in subroutine APPLYSYMC and to move the centre appropriately 
      inside subroutine BEAUTIFY.

3.    The recommended initial parameters for more accurate orientation parameter
      determination are:-
        a high resolution (e.g. 10-15 Angstroms)
        a reasonably but not excessively high value for ITMAX (e.g. 200)
        a reasonably but not excessively high value for IPMAX (e.g. 10)
        a high rather than a small value for DANG (i.e. 200.0 instead of 2.0)
        RBFACT should be set to 0 but different values can be tried to 
          decrease the weight of the noisier data at higher resolution 
          (e.g. use 2000 or 3000) but to keep them in.
      A range of values should always be tried since the optimal values may
        depend on the structure and the quality of the pictures.
