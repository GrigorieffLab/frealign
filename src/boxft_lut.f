C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      REAL FUNCTION   BOX FT_LUT (ARG,SINCLUT)
C**************************************************************************
C     Computes the Fourier transform of the indicator function of a
C     solid box. The argument ARG should be PI * DSTAR * WIDTH .
C     Used in AINTER03DS and PINSERT.
C     Uses a lookup table for sinc functions created in CALCSINC.
C**************************************************************************
      IMPLICIT   NONE
      INTEGER    I,J
      REAL       ARG(3), SINC, PROD
      REAL       C1, C2
      PARAMETER  ( C1 = 1.0 / 6.0 , C2 = 1.0 / 120.0 )
      REAL       SINCLUT(*), SEMICIRC
      PARAMETER  (SEMICIRC=180.0)
C**************************************************************************
C     If the argument is large enough, use the standard formula.
C     If it is small, use the Taylor series expansion.
C**************************************************************************
      PROD=1.0
      DO 10 I=1,3
      ARG(I)  =  ABS(ARG(I))
      IF ( ARG(I) .LT. 2.0E-02 ) THEN
         SINC  =  1.0
      ELSE
      	 J  = NINT(SEMICIRC*ARG(I))
         SINC  =  SINCLUT (J)
      END IF
      PROD=PROD*SINC
10    CONTINUE
      BOX FT_LUT = PROD

      RETURN
      END
