C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------
      SUBROUTINE PAD(NSAM,IPAD,DATA,OUTD,FFTW_PLANS)
C -----------------------------------------------------------
C     Pad particle by factor IPAD (particle origin at (1,1))
C -----------------------------------------------------------
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE
C
      INTEGER NSAM,I,J,II,JJ,ID,ID4,IPAD,NSAM4
      INTEGER NSAMH
      REAL CIRC,DATA(*),OUTD(*),SCAL,A
      TYPE(C_PTR) FFTW_PLANS(*)
C -----------------------------------------------------------
C
      NSAM4=IPAD*NSAM
      NSAMH=NSAM/2
      SCAL=1.0/NSAM/NSAM
      DO 75 I=1,NSAM4*(NSAM4+2)
        OUTD(I)=0.0
75    CONTINUE
      DO 71 I=1,NSAM*(NSAM+2)
        OUTD(I)=DATA(I)*SCAL
71    CONTINUE

      CALL FFTW_BWD(OUTD,OUTD,FFTW_PLANS(2))
C
      CIRC=0.0
      DO 70 I=1,NSAM
        DO 70 J=1,NSAM
          IF ((I.EQ.NSAM/2+1).OR.(J.EQ.NSAM/2+1)
     +       .OR.(I.EQ.NSAM/2).OR.(J.EQ.NSAM/2)) THEN
            ID=I+(NSAM+2)*(J-1)
            CIRC=CIRC+OUTD(ID)
          ENDIF
70    CONTINUE
      CIRC=CIRC/(4*NSAM-4)
C
      DO 72 J=1,NSAM4
        DO 72 I=1,NSAM4
          ID4=I+(NSAM4+2)*(J-1)
          OUTD(ID4)=OUTD(ID4)+CIRC
72    CONTINUE
      DO 73 J=NSAM,1,-1
        JJ=J
        IF (JJ.GT.NSAMH) JJ=JJ+NSAM4-NSAM
        DO 73 I=NSAM,1,-1
          II=I
          IF (II.GT.NSAMH) II=II+NSAM4-NSAM
          ID=I+(NSAM+2)*(J-1)
          ID4=II+(NSAM4+2)*(JJ-1)
          A=OUTD(ID)-CIRC
          OUTD(ID)=CIRC
          OUTD(ID4)=A
73    CONTINUE
C
      CALL FFTW_FWD(OUTD,OUTD,FFTW_PLANS(9))
C
      RETURN
      END
