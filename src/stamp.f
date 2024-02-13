C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE STAMP(A,IA,ID,DATA,N1,N2,IX,IY,CFORM)
C**************************************************************************
      IMPLICIT NONE

      INTEGER I,J,K,L,IDIGIT(3,5,16),IX,IY,IA,N1,N2,ID
      INTEGER IDATA(4*10,5)
      REAL DATA(*),A,MIN,MAX
      CHARACTER CFORM,LINE*20,CDIGIT*16

      DATA CDIGIT/"0123456789+-. °*"/
      DATA IDIGIT/1,1,1,1,0,1,1,0,1,1,0,1,1,1,1,
     +            0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,
     +            1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,
     +            1,1,1,0,0,1,1,1,1,0,0,1,1,1,1,
     +            1,0,1,1,0,1,1,1,1,0,0,1,0,0,1,
     +            1,1,1,1,0,0,1,1,1,0,0,1,1,1,1,
     +            1,1,1,1,0,0,1,1,1,1,0,1,1,1,1,
     +            1,1,1,0,0,1,0,0,1,0,0,1,0,0,1,
     +            1,1,1,1,0,1,1,1,1,1,0,1,1,1,1,
     +            1,1,1,1,0,1,1,1,1,0,0,1,1,1,1,
     +            0,0,0,0,1,0,1,1,1,0,1,0,0,0,0,
     +            0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,
     +            0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
     +            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     +            0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,
     +            1,0,1,1,1,1,1,0,1,1,1,1,1,0,1/
C**************************************************************************
      IF (IA.EQ.0) RETURN
      I=1
      IF (ABS(A).GT.1.0) I=INT(LOG10(ABS(A)))+1
C      IF (A.LT.0.0) I=I+1
      IF (ID.EQ.2) I=I+1
      IF ((I.GT.IA).OR.(IA.GT.10)) THEN
        LINE(1:20)='--------------------'
      ELSE
        IF (ID.EQ.0) THEN
          WRITE(LINE,1000) INT(A)
          LINE(1:IA)=LINE(20-IA+1:20)
        ELSE
C          IF (ID.EQ.2) I=I-1
          WRITE(LINE,1010) A
          LINE(1:IA)=LINE(11-I:11-I+IA-1)
          IF (ID.EQ.2) LINE(IA:IA)=CDIGIT(15:15)
        ENDIF
      ENDIF
1000  FORMAT(I20)
1010  FORMAT(F20.9)

      DO 10 I=1,40
        DO 10 J=1,5
          IDATA(I,J)=0
10    CONTINUE
      DO 20 L=1,IA
        DO 30 K=1,16
          IF (CDIGIT(K:K).EQ.LINE(L:L)) GOTO 40
30      CONTINUE
40      CONTINUE
        DO 50 I=1,5
          DO 50 J=1,3
            IDATA(J+4*L-3,I)=IDIGIT(J,I,K)
50      CONTINUE
20    CONTINUE

      MIN=1.0E30
      MAX=-1.0E30
      DO 60 I=1,N1*N2
        IF(DATA(I).GT.MAX) MAX=DATA(I)
        IF(DATA(I).LT.MIN) MIN=DATA(I)
60    CONTINUE
      DO 70 I=1,4*IA
        IF (IX+I-1.LE.N1) THEN
          DO 80 J=1,5
            K=IX+I-1+N1*(IY+J-1)
            IF (CFORM.EQ."M") THEN
              DATA(K)=IDATA(I,6-J)*(MAX-MIN)+MIN
            ELSE
              DATA(K)=IDATA(I,J)*(MAX-MIN)+MIN
            ENDIF
80        CONTINUE
        ENDIF
70    CONTINUE
      RETURN
      END
