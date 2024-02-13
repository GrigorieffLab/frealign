C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------
      SUBROUTINE WINDOW(NSAM,NSAMR,DATA,STD)
C -----------------------------------------------------------
C     WINDOWS PARTICLE FROM NSAM x NSAM TO NSAMR x NSAMR
C     Used in subroutine MATCH, CTFAPPLY, CTFAPPLY_PHASE_ONLY
C -----------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER I,J,M,NSAM,NSAMR,IDR,ID
      REAL DATA(*),STD,MEAN
C -----------------------------------------------------------
      STD=0.0
      MEAN=0.0
      M=(NSAM-NSAMR)/2
      DO 10 J=1,NSAMR
        DO 10 I=1,NSAMR
          IDR=I+(J-1)*NSAMR
          ID=I+M+(J+M-1)*(NSAM+2)
          DATA(IDR)=DATA(ID)
          STD=STD+DATA(IDR)**2
          MEAN=MEAN+DATA(IDR)
10    CONTINUE
      STD=STD/NSAMR/NSAMR
      MEAN=MEAN/NSAMR/NSAMR
      STD=SQRT(STD-MEAN**2)
C
      RETURN
      END
