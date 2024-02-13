C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
	SUBROUTINE CARD4(IFIRST,ILAST)
C**************************************************************************
C Card 4: global parameters
C         IFIRST,ILAST - FIRST, LAST PARTICLES
C Called by Frealign.
C*************************************************************************
        IMPLICIT NONE
        INTEGER IFIRST,ILAST
C*************************************************************************
      	WRITE(*,*)' FIRST, LAST PARTICLES ?'
        READ(*,*)IFIRST,ILAST
      	WRITE(*,17003)IFIRST,ILAST
        IF (IFIRST.GT.ILAST)
     +    STOP ' ERROR: IFIRST must be smaller'//
     +         ' or equal to ILAST (Card 4)'
        IF (IFIRST.LE.0)
     +    STOP ' ERROR: IFIRST must be greater than 0 (Card 4)'
17003	FORMAT(2I8)
	RETURN
	END
