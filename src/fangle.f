C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      REAL FUNCTION FANGLE(SIG2,A,AM,SA,STIF)
C**************************************************************************
C  Calculates probability distribution for an angle A
C  Used in CALCFX and PREFINE
C**************************************************************************
      IMPLICIT NONE
      REAL AM,SA,A,SIG2,D,PI,STIF
      PARAMETER (PI=3.1415926535897)
C**************************************************************************
C
      D=ABS(A-AM)
10    CONTINUE
      IF (D.GT.PI) THEN
        D=ABS(D-2.0*PI)
        GOTO 10
      ENDIF      
      FANGLE=SIG2*(-(STIF*D)**2/2.0/SA**2)
C
      RETURN
      END
