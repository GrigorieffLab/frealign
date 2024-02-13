C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE CORRECT3D(NSAM,SINCLUT,A3DV,INTERP,IPAD)
C**************************************************************************
C     Corrects 3D volume for fall-off at edges due to
C     interpolation kernel in PEXTRACT
C     Used in frealign
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER NSAM,I,J,K,NMID,ID,INTERP,IPAD
      REAL SINCLUT(*),A3DV(*),A(3),BOXFT_LUT,SCAL
C
      NMID=NSAM/2+1
      SCAL=1.0/NSAM
      IF (IPAD.NE.0) SCAL=SCAL/IPAD
C
      IF (INTERP.EQ.0) THEN
C
        DO 10 I=1,NSAM
          A(1)=REAL(I-NMID)*SCAL
          DO 10 J=1,NSAM
            A(2)=REAL(J-NMID)*SCAL
            DO 10 K=1,NSAM
              A(3)=REAL(K-NMID)*SCAL
             ID=I+(NSAM+2)*(J-1+NSAM*(K-1))
              A3DV(ID)=A3DV(ID)/BOXFT_LUT(A,SINCLUT)
10      CONTINUE
C
      ELSE
C
        DO 11 I=1,NSAM
          A(1)=REAL(I-NMID)*SCAL
          DO 11 J=1,NSAM
            A(2)=REAL(J-NMID)*SCAL
            DO 11 K=1,NSAM
              A(3)=REAL(K-NMID)*SCAL
             ID=I+(NSAM+2)*(J-1+NSAM*(K-1))
              A3DV(ID)=A3DV(ID)/BOXFT_LUT(A,SINCLUT)**2
11      CONTINUE
C
      ENDIF
C
      RETURN
      END
