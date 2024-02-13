C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      REAL FUNCTION REMAP_THETA(THETA)
C**************************************************************************
C  Remaps THETA angle to sample more around THETA=90deg
C  Used in LMAIN
C**************************************************************************
      IMPLICIT NONE
      REAL THETA,PI,D,T
      PARAMETER (PI=3.1415926535897,D=1.0E-3)
C**************************************************************************
C
      T=ABS(THETA-NINT(THETA/2.0/PI)*2.0*PI)
C
      IF (ABS(T).LT.D) THEN
        REMAP_THETA=0.0
      ELSEIF (ABS(ABS(T)-PI).LT.D) THEN
        REMAP_THETA=PI
      ELSE
        REMAP_THETA=-LOG10(1.0/TAN(T)+1.0/SIN(T))/2.0+PI/2.0
      ENDIF
C
      IF (THETA.LT.0.0) REMAP_THETA=-REMAP_THETA
C
      RETURN
      END
