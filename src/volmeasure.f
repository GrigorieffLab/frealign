C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE VOLMEASURE(NSAM,A3DV,FPART)
C**************************************************************************
C Calculated approximate particle volume (fraction) by thresholding map
C Used by FREALIGN
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,I,J,I1,I2,JMAX,JC,NSAM3,ID,IS,JJ
      PARAMETER (JMAX=10)
      REAL XSTD,A3DV(*),THRESH,FPART
      REAL D1,D2,DATA
      PARAMETER (XSTD=5.0)
C**************************************************************************
C
      JC=NSAM/2+1
      NSAM3=NSAM*NSAM*NSAM
C
C     Threshold filtered map
      THRESH=0.0
      DO 40 J=1,10
        D1=0.0
        D2=0.0
        I1=0
        I2=0
        DO 30 JJ=0,NSAM*NSAM-1
          IS=JJ*(NSAM+2)
          DO 30 I=1,NSAM
            ID=IS+I
            DATA=A3DV(ID)
            IF (DATA.GE.THRESH) THEN
              I2=I2+1
              D2=D2+DATA
            ELSE
              I1=I1+1
              D1=D1+DATA
            ENDIF
30      CONTINUE
        IF (I1.NE.0) D1=D1/I1
        IF (I2.NE.0) D2=D2/I2
        THRESH=(D2-D1)/2.0+D1
40    CONTINUE
      THRESH=(D2-D1)/2.0+D1
C      WRITE(*,*)
C      WRITE(*,*) 'Threshold = ',THRESH
C
C      THRESH=(D2-D1)/2.0+D1-(D2-D1)/XSTD
      THRESH=(D2-D1)/2.0+D1
      I2=0
      DO 31 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 31 I=1,NSAM
          ID=IS+I
          DATA=A3DV(ID)
          IF (DATA.GE.THRESH) THEN
            I2=I2+1
          ENDIF
31    CONTINUE
C
      FPART=REAL(I2)/NSAM3
C
      RETURN
      END
