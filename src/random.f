C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      REAL FUNCTION RANDOM(IS)
C**************************************************************************
C     Random number algorithm after L'Ecuyer for uniform distribution
C     Used in LMAIN
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER IS,I1,I2,I3
C**************************************************************************
C
      IF (IS.LT.0) THEN
        IS=-IS*127773
      ENDIF
      I1=IS/127773
      I2=MOD(IS,127773)
      I3=16807*I2-2836*I1
      IF (I3.LE.0) THEN
        IS=I3+2147483647
      ELSE
        IS=I3
      ENDIF
      RANDOM=REAL(IS)/REAL(2147483647)
C
      RETURN
      END
