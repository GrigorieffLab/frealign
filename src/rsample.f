C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
      PROGRAM RSAMPLE
C
      IMPLICIT NONE
C
      INTEGER NMAX,N,IRAN,I,J,K,NBOOT,CNT
      PARAMETER (NMAX=2000000)
      INTEGER F(NMAX),SLEN,L,SLEN2,M(NMAX),LGP(NMAX)
      INTEGER ILST(NMAX)
      REAL A(NMAX),B(NMAX),C(NMAX),X(NMAX),PS(NMAX)
      REAL Y(NMAX),MM,D1(NMAX),D2(NMAX),S(NMAX)
      REAL AS(NMAX),P(NMAX),DP(NMAX),RANDOM
      REAL OCC(NMAX),OMAX,PSIZE
      CHARACTER LINE*200,ONAME*200,CN*200,VX*15
      LOGICAL FOLD,FREALIGNX
C**************************************************************************
C
      DATA  VX/'1.02 - 26.03.17'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C
      WRITE(*,1010)VX
1010  FORMAT(/'  RSAMPLE ',A15,/,
     +       /'  Copyright 2013 Howard Hughes Medical Institute.',
     +       /'  All rights reserved.',
     +       /'  Use is subject to Janelia Farm Research Campus',
     +        ' Software Copyright 1.1',
     +       /'  license terms ( http://license.janelia.org/license',
     +        '/jfrc_copyright_1_1.html )'/)
C
      WRITE(*,*)' INPUT PARAMETER FILE ?'
      READ(*,7006)LINE
      WRITE(*,*)LINE(1:SLEN(LINE))
7006  FORMAT(A200)  
      OPEN(77,FILE=LINE,STATUS='OLD')   
C
      WRITE(*,*)' PIXEL SIZE [A] ?'
      READ(*,*)PSIZE
      WRITE(*,*)' # OF BOOTSTRAP VOLUMES ?'
      READ(*,*)NBOOT
C
      N=0
      FOLD=.FALSE.
      FREALIGNX=.FALSE.
7009  CONTINUE
      READ(77,7006,END=7008)LINE
      IF (LINE(1:1).EQ.'C') GOTO 7009
      N=N+1
      IF (N.GT.NMAX) STOP 'INCREASE NMAX'
C
97    CONTINUE
      IF (FOLD) THEN
        READ(LINE,7015)ILST(N),A(N),B(N),
     +     C(N),X(N),Y(N),MM,F(N),
     +     D1(N),D2(N),AS(N),P(N),DP(N)
7015    FORMAT(I7,5F8.2,F8.0,I6,2F9.1,F8.2,F7.2,F8.2)
        M(N)=NINT(MM)
        LGP(N)=5000
        S(N)=1.0
      ELSE
        IF (SLEN(LINE).LE.95) THEN
          GOTO 99
        ELSEIF (SLEN(LINE).LE.128) THEN
          READ(LINE,7012,ERR=99,IOSTAT=CNT)ILST(N),A(N),B(N),
     +       C(N),X(N),Y(N),M(N),F(N),
     +       D1(N),D2(N),AS(N),MM,LGP(N),P(N),DP(N)
7012      FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I13,2F8.2)
          IF (CNT.NE.0) GOTO 99
          S(N)=1.0
          GOTO 98 
        ELSEIF (SLEN(LINE).LE.140) THEN
          IF (FREALIGNX) GOTO 99
          READ(LINE,7011,ERR=99,IOSTAT=CNT)ILST(N),A(N),B(N),
     +       C(N),X(N),Y(N),M(N),F(N),
     +       D1(N),D2(N),AS(N),MM,LGP(N),S(N),P(N),DP(N)
7011      FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I10,F11.4,2F8.2)
          IF (CNT.NE.0) GOTO 99
C          IF ((LGP(N).EQ.0).OR.(S(N).EQ.0.0)) GOTO 99
          GOTO 98
        ELSE
          IF (.NOT.FREALIGNX) THEN
            WRITE(*,*) ' FrealignX parameter file'
            WRITE(*,*)
            FREALIGNX=.TRUE.
          ENDIF
          READ(LINE,7013,ERR=99,IOSTAT=CNT)ILST(N),A(N),B(N),
     +       C(N),X(N),Y(N),M(N),F(N),
     +       D1(N),D2(N),AS(N),PS(N),MM,LGP(N),S(N),P(N),DP(N)
7013      FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,3F8.2,I10,F11.4,2F8.2)
          IF (CNT.NE.0) GOTO 99
C          IF ((LGP(N).EQ.0).OR.(S(N).EQ.0.0)) GOTO 99
          GOTO 98
        ENDIF
99      CONTINUE
        IF (N.EQ.1) THEN
          WRITE(*,*) ' Error reading par file.',
     +               ' Trying pre v9 par file format...'
          FOLD=.TRUE.
          GOTO 97
        ELSE
          WRITE(*,*)  
     +      ' ERROR: Something wrong in parameter file'
          WRITE(*,*) LINE
          STOP
        ENDIF
98      CONTINUE
      ENDIF
C
      IF (FOLD) THEN
        X(N)=X(N)*PSIZE
        Y(N)=Y(N)*PSIZE
      ENDIF
      GOTO 7009
7008  CONTINUE
      CLOSE(77)
      WRITE(*,*)' N = ',N
C
      WRITE(*,*)' OUTPUT PARAMETER FILE (*.par)?'
      READ(*,7006)LINE
      WRITE(*,*)LINE(1:SLEN(LINE))
      L=SLEN(LINE)
      IRAN=-100
      OCC(1)=RANDOM(IRAN)
      DO 10 I=1,NBOOT
        DO 20 J=1,N
          OCC(J)=0.0
20      CONTINUE
        OMAX=0.0
        DO 30 J=1,N
          K=INT(N*RANDOM(IRAN))+1
          OCC(K)=OCC(K)+1.0
          IF (OMAX.LT.OCC(K)) OMAX=OCC(K)
30      CONTINUE
        DO 40 J=1,N
          OCC(J)=OCC(J)/OMAX*100.0
40      CONTINUE
        WRITE(CN,7010)I
7010    FORMAT(I10)
        ONAME=LINE(1:L-4)//CN(SLEN2(CN):10)//LINE(L-3:L)
        WRITE(*,*)ONAME(1:SLEN(ONAME))
        OPEN(77,FILE=ONAME,STATUS='UNKNOWN')
        CLOSE(77,STATUS='DELETE')
        OPEN(77,FILE=ONAME,STATUS='NEW')
        DO 50 J=1,N
          IF (FREALIGNX) THEN
            WRITE(77,7013)ILST(J),A(J),B(J),C(J),
     +          X(J),Y(J),M(J),F(J),D1(J),D2(J),AS(J),
     +          PS(J),OCC(J),LGP(J),S(J),P(J),DP(J)
          ELSE
            WRITE(77,7011)ILST(J),A(J),B(J),C(J),
     +          X(J),Y(J),M(J),F(J),D1(J),D2(J),
     +          AS(J),OCC(J),LGP(J),S(J),P(J),DP(J)
          ENDIF
50      CONTINUE
        CLOSE(77)
10    CONTINUE
C
      WRITE(*,*)
      WRITE(*,*) ' Normal termination of rsample'
C
      END
C**************************************************************************
      REAL FUNCTION RANDOM(IS)
C**************************************************************************
C     Random number algorithm after L'Ecuyer for uniform distribution
C     Used in LMAIN
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER IS,I1,I2,I3
C**************************************************************************
C
      IF (IS.LT.0) THEN
        IS=-IS*127773
      ENDIF
      I1=IS/127773
      I2=MOD(IS,127773)
      I3=16807*I2-2836*I1
      IF (I3.LE.0) THEN
        IS=I3+2147483647
      ELSE
        IS=I3
      ENDIF
      RANDOM=REAL(IS)/REAL(2147483647)
C
      RETURN
      END
C**************************************************************************
      INTEGER FUNCTION SLEN(STRG)
C**************************************************************************
C Determines length of string STRG
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER I
      CHARACTER*200 STRG
C**************************************************************************
C
      DO 10 I=200,1,-1
        IF (STRG(I:I).NE." ") GOTO 20
10    CONTINUE
C
20    CONTINUE
      SLEN=I
C
      RETURN
      END
C**************************************************************************
      INTEGER FUNCTION SLEN2(STRG)
C**************************************************************************
C Determines the first non-blank char in a string STRG
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER I
      CHARACTER*200 STRG
C**************************************************************************
C
      DO 10 I=1,200,1
        IF (STRG(I:I).NE." ") GOTO 20
10    CONTINUE
C
20    CONTINUE
      SLEN2=I
C
      RETURN
      END
