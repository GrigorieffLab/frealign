C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      FUNCTION PDIFF(C1,C2)
C**************************************************************************
C     Phase calculation done if complexx number non-zero
C     used in PINSERT, PRESB, PRES, and SHELTEST
C**************************************************************************
      IMPLICIT NONE
      REAL PDIFF
      COMPLEX C1,C2,C
C**************************************************************************
      C=C1*CONJG(C2)
      IF(C.EQ.(0.,0.)) THEN
        PDIFF=0.
          ELSE
        PDIFF=ABS(ATAN2(AIMAG(C),REAL(C)))
      ENDIF
      RETURN
      END
