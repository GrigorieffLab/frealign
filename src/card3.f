C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
	SUBROUTINE CARD3(PMASK,DMASK)
C**************************************************************************
C Card 3: global parameters
C       PMASK - [0/1] parameters to include in refinement [1,1,1,1,1]
C       DMASK - [X, Y, Z, RADIUS] density mask for focused classification
C Called by Frealign.
C*************************************************************************
        IMPLICIT NONE
	INTEGER J,PMASK(5)
        REAL DMASK(4)
        CHARACTER*200 INLINE
C*************************************************************************
C      	WRITE(*,*)' PMASK for parameter refinement (e.g. 1,1,1,1,1)'
      	WRITE(*,*)' PMASK for parameter refinement (e.g. 1,1,1,1,1),',
     +            ' DMASK for focused classification',
     +            ' (e.g. 100.0,80.0,120.0,20.0)?'
        READ(*,'(A)')INLINE
C      	READ(INLINE,*,ERR=99, END=99) (PMASK(J),J=1,5)
      	READ(INLINE,*,ERR=99, END=99) (PMASK(J),J=1,5),(DMASK(J),J=1,4)
C      	WRITE(*,17023) (PMASK(J),J=1,5)
      	WRITE(*,17023) (PMASK(J),J=1,5),(DMASK(J),J=1,4)
17023	FORMAT(5I4,4F8.2)
      	DO 6901 J=1,4
      	  IF(DMASK(J).LE.0.0) STOP ' invalid DMASK values'
C          DMASK(J)=-1.0
6901	CONTINUE
        GOTO 98
99      CONTINUE
      	READ(INLINE,*) (PMASK(J),J=1,5)
      	WRITE(*,17024) (PMASK(J),J=1,5)
17024	FORMAT(5I4)
        WRITE(*,*) 'CARD 3: No DMASK parameters; using default values'
        DMASK(1)=-1.0
98      CONTINUE
      	DO 6900 J=1,5
      	  IF(PMASK(J).NE.0.AND.PMASK(J).NE.1)
     +      STOP ' invalid PMASK values'
6900	CONTINUE
C       PSI/PHI are in different order in parameter file and inside program
        J=PMASK(3)
        PMASK(3)=PMASK(1)
        PMASK(1)=J
C
	RETURN
	END
