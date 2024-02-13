C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PINC(NSET,THRESH,OCC,PRESA,NANG1,NANG)
C**************************************************************************
C     Determines the score threshold to exclude a given
C     percentage of particles
C     Used in card10
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER NSET,NANG1,NANG,I,J
      REAL THRESH(*),PRESA(*),T,P,C,OCC(*),N,OT
C
      OT=0.0
      DO 40 I=NANG1+1,NANG
        OT=OT+OCC(I)
40    CONTINUE
      OT=OT/REAL(NANG-NANG1)
C
      C=1.0/10000.0
      DO 20 J=20000,1,-1
        N=0.0
        T=(REAL(J)-10000.0)*C
        DO 10 I=NANG1+1,NANG
          IF (PRESA(I).GT.T) N=N+OCC(I)
10      CONTINUE
        P=N/REAL(NANG-NANG1)
        IF (P.GE.THRESH(NSET)*OT) GOTO 30
20    CONTINUE
C
30    CONTINUE
      IF (N.EQ.0.0) THEN
        WRITE(*,*)
     + ' ERROR: all particles excluded in data set'
        STOP
     + ' Try changing score threshold/percentage (Card 6)'
      ELSE
        T=T*100.0
        IF (J.EQ.1) THEN
          T=-100.0
          P=1.0
        ENDIF
        WRITE(*,1000)NSET,P*100.0,NINT(N),T
1000    FORMAT(/,' NSET =',I3,': percentage included = ',F7.2,
     +         /,'            particles included  = ',I7,
     +         /,'            threshold set to   ',F10.2,/)
      ENDIF
      THRESH(NSET)=T
C
      RETURN
      END
