C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE CHECKSYM(SYMOP,SYMCOUNT)
C**************************************************************************
C  Checks that repeated application of the symmetry matrix will rotate back
C    to the identity matrix (within error) to make sure the symmetry is 
C    possible - tests up to 50-fold rotational symmetry.
C  Used in CARD5 and GETSYMMAT.
C**************************************************************************
      IMPLICIT NONE
      EXTERNAL MATMUL

      INTEGER I,J,SYMCOUNT
      REAL SYMOP(3,3),TM(3,3),DIFF,DELTA
      PARAMETER (DELTA=0.01)
C**************************************************************************
      DO 10 I=1,3
      	DO 10 J=1,3
          TM(I,J)=0.0
10    CONTINUE
      TM(1,1)=1.0
      TM(2,2)=1.0
      TM(3,3)=1.0

      SYMCOUNT=0
20    CONTINUE
      CALL MATMUL(TM,SYMOP,TM)
      SYMCOUNT=SYMCOUNT+1
      DIFF=0.0
      DO 30 I=1,3
      	DO 30 J=1,3
      	IF (I.EQ.J) THEN
      	  DIFF=DIFF+ABS(TM(I,J)-1.0)
      	ELSE
      	  DIFF=DIFF+ABS(TM(I,J))
      	ENDIF
30    CONTINUE
      IF (DIFF.LT.DELTA) GOTO 99
      IF (SYMCOUNT.GT.50) THEN
      	SYMCOUNT=-1
      	GOTO 99
      ENDIF
      GOTO 20

99    CONTINUE
      RETURN
      END
