C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE FIND_RLIM(N,DATA,NA,THRESH,K)
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER N,NA,I,J,K,N1,N2
      REAL DATA(*),THRESH,S
C**************************************************************************
C
      N1=NA/2
      N2=NA/2
      IF (N2+N1+1.GT.NA) N1=N1-1
      DO 10 K=3,N
        S=0.0
        I=0
        DO 20 J=K-N1,K+N2
          IF ((J.GE.3).AND.(J.LE.N)) THEN
            S=S+DATA(J)
            I=I+1
          ENDIF
20      CONTINUE
        IF (I.NE.0) S=S/I
        IF (S.LT.THRESH) GOTO 99
10    CONTINUE
C
99    CONTINUE
      K=K-1
C
      RETURN
      END
