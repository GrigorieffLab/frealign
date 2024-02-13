C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C************************************************************************
      REAL FUNCTION ERC(X)
C************************************************************************
C     Calculates complementary error function
C     Used in subroutine OPRESSTATHALF
C************************************************************************
      IMPLICIT NONE
C
      INTEGER I
      REAL X,Y,Z,C(10)
      DATA C/0.17087277,-0.82215223,1.48851587,-1.13520398,
     +       0.27886807,-0.18628806,0.09678418,0.37409196,
     +       1.00002368,-1.26551223/
C************************************************************************
      Y=1.0/(1.0+0.5*ABS(X))
      Z=C(1)
      DO 10 I=2,10
        Z=Y*Z+C(I)
10    CONTINUE
      Z=Y*EXP(-ABS(X)**2+Z)
C
      IF (X.LT.0.0) Z=2.0-Z
C
      ERC=Z
C
      RETURN
      END
