C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE CALCSINC(SINCLUT,N)
C**************************************************************************
C  Calculates table of sinc function in 1 deg steps out to N=2000 deg
C**************************************************************************
      IMPLICIT NONE
      REAL PI,SINCLUT(*)
      INTEGER J,N
      PARAMETER  (PI=3.1415926535897)

      DO 100 J=1,N
      	SINCLUT(J)=SIN(FLOAT(J)*PI/180.0)/(FLOAT(J)*PI/180.0)
100   CONTINUE

      RETURN
      END
