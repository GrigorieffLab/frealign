C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
	SUBROUTINE CARDS15TO18(F3D,FWEIGH,F3D1,F3D2,FPHA,FPOI,NDOC1)
C**************************************************************************
C Cards 15 to 18 are filenames for 3D diagnostic output files (18a is backup)
C
C Card 15         F3D1   - 3D ECONSTRUCTION HALFSET 1 FOR OUTPUT
C Card 16         F3D2   - 3D ECONSTRUCTION HALFSET 2 FOR OUTPUT
C Card 17         FPHA   - 3D PHASE RESIDUAL FILE FOR OUTPUT, Output 3D file with
C                         average phase differences for each voxel (in Fourier
C                         space) between images and reference.
C Card 18         FPOI   - 3D POINT SPREAD FUNCTION FOR OUTPUT, Output 3D file
C                         with point spread function indicating anisotropies
C                         in resolution (in real space).
C Called by Frealign.
C**************************************************************************
        IMPLICIT NONE
        INTEGER NDOC1,SLEN2
	CHARACTER*200 F3D,FWEIGH,F3D1,F3D2,FPHA,FPOI
C**************************************************************************
	  WRITE(*,*)' 3D RECONSTRUCTION HALFSET 1 FOR OUTPUT ?'
          READ(*,7006)F3D1
          WRITE(*,*)' 3D RECONSTRUCTION HALFSET 2 FOR OUTPUT ?'
          READ(*,7006)F3D2
          WRITE(*,*)' 3D PHASE RESIDUAL FILE FOR OUTPUT ?'
          READ(*,7006)FPHA
          WRITE(*,*)' 3D POINT SPREAD FUNCTION FOR OUTPUT ?'
          READ(*,7006)FPOI
7006      FORMAT(A200)
        WRITE(NDOC1+1,6705)F3D(1:SLEN2(F3D)),
     +    FWEIGH(1:SLEN2(FWEIGH)),F3D1(1:SLEN2(F3D1)),
     +    F3D2(1:SLEN2(F3D2)),FPHA(1:SLEN2(FPHA)),FPOI(1:SLEN2(FPOI))
6705    FORMAT(  'C 3D reconstruction file      ',A,
     +          /'C 3D weights file             ',A,
     +          /'C 3D reconstruction halfset 1 ',A,
     +          /'C 3D reconstruction halfset 2 ',A,
     +          /'C 3D ave phase residual file  ',A,
     +          /'C 3D point spread function    ',A/'C')
	RETURN
	END
