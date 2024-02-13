C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      PROGRAM SET_POLARITY
C**************************************************************************
C Program to determine filament polarity from input parameter file.
C Generates new parameter file with updated PSI angles.
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER MAXP,I,J,K,SLEN2,IERR,MAXH,NH,NF,NP
      INTEGER NFILM,N,NMAX,JMAX,I2,CNT
      INTEGER,ALLOCATABLE :: ILIST(:),FILM(:),LGP(:)
      INTEGER,ALLOCATABLE :: IFILM(:),ABSMAGP(:)
      REAL,ALLOCATABLE :: PSI(:),THETA(:),PHI(:)
      REAL,ALLOCATABLE :: SHX(:),SHY(:),OCC(:)
      REAL,ALLOCATABLE :: DFMID1(:),DFMID2(:),S(:)
      REAL,ALLOCATABLE :: DPRES(:),ANGAST(:),PRES(:)
      REAL PSIVAR,D,PSIAVE,AVPRES,PFLIP,AMAG
      CHARACTER*200,ALLOCATABLE :: LINE(:)
      CHARACTER FPAR*200,VX*15
      LOGICAL EX,FOLD
C**************************************************************************
C
      DATA  VX/'1.02 - 26.05.15'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C
      WRITE(*,1010)VX
1010  FORMAT(/'  SET_POLARITY ',A15,/,
     +       /'  Copyright 2013 Howard Hughes Medical Institute.',
     +       /'  All rights reserved.',
     +       /'  Use is subject to Janelia Farm Research Campus',
     +        ' Software Copyright 1.1',
     +       /'  license terms ( http://license.janelia.org/license',
     +        '/jfrc_copyright_1_1.html )'/)
C
C
C  Read input files
C
      WRITE(*,2001)
2001  FORMAT(' PSI angle variation?')
      READ(*,*)PSIVAR
      WRITE(*,3000)PSIVAR
3000  FORMAT(F12.3,/)
C
      WRITE(*,2002)
2002  FORMAT(' Flip polarity if ave. SCORE below?')
      READ(*,*)PFLIP
      WRITE(*,3000)PFLIP
C
      WRITE(*,4000)
4000  FORMAT(' Input parameter file?')
      READ(*,7007)FPAR
      WRITE(*,8000)FPAR(1:SLEN2(FPAR))
      OPEN(10,FILE=FPAR,STATUS='OLD')
C
      WRITE(*,7000)
7000  FORMAT(/' Output parameter file?')
      READ(*,7007)FPAR
      WRITE(*,8000)FPAR(1:SLEN2(FPAR))
      INQUIRE(FILE=FPAR,EXIST=EX)
      IF (EX) THEN
        OPEN(11,FILE=FPAR,STATUS='OLD')
        CLOSE(11,STATUS='DELETE')
      ENDIF
      OPEN(11,FILE=FPAR,STATUS='NEW')
C
      MAXP=1
      MAXH=1
      I=-1
      NFILM=0
      FOLD=.FALSE.
7006  CONTINUE
      READ(10,8000,END=7005)FPAR
      IF (FPAR(1:1).NE.'C') THEN
        MAXP=MAXP+1
        READ(FPAR,7012)J
7012    FORMAT(55X,I6)
        IF (I.NE.J) THEN
          I=J
          NFILM=NFILM+1
        ENDIF
      ELSE
        MAXH=MAXH+1
      ENDIF
      GOTO 7006
C
7005  CONTINUE
      REWIND(UNIT=10)
C
      ALLOCATE(ILIST(MAXP),FILM(MAXP),PSI(MAXP),
     +       THETA(MAXP),PHI(MAXP),OCC(MAXP),
     +       SHX(MAXP),SHY(MAXP),ABSMAGP(MAXP),
     +       DFMID1(MAXP),DFMID2(MAXP),PRES(MAXP),
     +       DPRES(MAXP),ANGAST(MAXP),IFILM(MAXP),
     +       LINE(2*MAXH),LGP(MAXP),S(MAXP),STAT=IERR)
      IF (IERR.NE.0) THEN
        STOP ' ERROR: Memory allocation failed'
      ENDIF
C
      NH=1
      NF=1
      NP=1
      I=-1
      NFILM=0
      IFILM(1)=1
7009  CONTINUE
      READ(10,7007,END=7008)LINE(NF)
7007  FORMAT(A200)
      IF (LINE(NF)(1:1).EQ.'C') THEN
        NF=NF+1
        GOTO 7009
      ENDIF
      NH=NF-1
97    CONTINUE
      IF (FOLD) THEN
        READ(LINE(NF),7014)ILIST(NP),PSI(NP),
     +      THETA(NP),PHI(NP),SHX(NP),
     +      SHY(NP),AMAG,FILM(NP),
     +      DFMID1(NP),DFMID2(NP),ANGAST(NP),
     +      PRES(NP),DPRES(NP)
7014    FORMAT(I7,5F8.2,F8.0,I6,2F9.1,F8.2,F7.2,F8.2)
            ABSMAGP(NP)=AMAG
            OCC(NP)=100.0
            LGP(NP)=5000
            S(NP)=1.0
      ELSE
        IF (SLEN2(LINE(NF)).LE.95) THEN
          GOTO 99
        ELSEIF (SLEN2(LINE(NF)).LE.128) THEN
          READ(LINE(NF),7015,ERR=99,IOSTAT=CNT)ILIST(NP),
     +      PSI(NP),THETA(NP),PHI(NP),SHX(NP),
     +      SHY(NP),ABSMAGP(NP),FILM(NP),
     +      DFMID1(NP),DFMID2(NP),ANGAST(NP),
     +      OCC(NP),LGP(NP),PRES(NP),DPRES(NP)
7015      FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I13,2F8.2)
          IF (CNT.NE.0) GOTO 99
          S(NP)=1.0
          GOTO 98
        ELSE
          READ(LINE(NF),7011,ERR=99,IOSTAT=CNT)ILIST(NP),
     +      PSI(NP),THETA(NP),PHI(NP),SHX(NP),
     +      SHY(NP),ABSMAGP(NP),FILM(NP),
     +      DFMID1(NP),DFMID2(NP),ANGAST(NP),
     +      OCC(NP),LGP(NP),S(NP),PRES(NP),DPRES(NP)
7011      FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I10,F11.4,2F8.2)
          IF (CNT.NE.0) GOTO 99
          IF ((LGP(NP).EQ.0).OR.(S(NP).EQ.0.0)) GOTO 99
          GOTO 98
        ENDIF
99      CONTINUE
        IF (NP.EQ.1) THEN
          WRITE(*,*) ' Error reading par file.',
     +               ' Trying pre v9 par file format...'
          FOLD=.TRUE.
          GOTO 97
        ELSE
          WRITE(*,*)
     +      ' ERROR: Something wrong in parameter file'
          WRITE(*,*) LINE(NF)
          STOP
        ENDIF
98      CONTINUE
      ENDIF
C
      IF(ILIST(NP).LE.0) THEN
        WRITE(*,*) ' ***WARNING*** blank line in parameter file'
        GO TO 7009
      ENDIF
      IF (I.NE.FILM(NP)) THEN
        I=FILM(NP)
        NFILM=NFILM+1
        IFILM(NFILM)=NP
      ENDIF
      NP=NP+1
      GOTO 7009
C
7008  CONTINUE
      NF=NF-NH-1
      NP=NP-1
      WRITE(*,1000)NH,NF,NP,NFILM
1000  FORMAT(/' No. of header lines read = ',I10,/,
     +        ' No. of footer lines read = ',I10,/,
     +        ' No. of particles read    = ',I10,/,
     +        ' No. of filaments         = ',I10,/)
C
      CLOSE(10)
C
C
C  Determine new PSI angles
C
      DO 10 I=1,NFILM
        JMAX=0
        NMAX=0
        AVPRES=0.0
        IF (I.NE.NFILM) THEN
          I2=IFILM(I+1)-1
        ELSE
          I2=NP
        ENDIF
        DO 20 J=0,359
          N=0
          DO 30 K=IFILM(I),I2
            D=ABS(PSI(K)-J)
            IF (D.GT.180.0) D=ABS(D-360.0)
            IF (D.LE.PSIVAR/2.0) N=N+1
            IF (J.EQ.0) AVPRES=AVPRES+PRES(K)
30        CONTINUE
          IF (N.GT.NMAX) THEN
            NMAX=N
            JMAX=J
          ENDIF
20      CONTINUE
        AVPRES=AVPRES/(I2-IFILM(I)+1)
        WRITE(*,5000) FILM(IFILM(I)),JMAX,AVPRES,
     +          REAL(NMAX)/(I2-IFILM(I)+1)
5000    FORMAT(' Best PSI angle for filament ',I9,
     +    ' = ',I5,' deg,  ave. SCORE = ',F7.2,', PSIscore = ',F6.2)
        DO 40 K=IFILM(I),I2
          D=ABS(PSI(K)-JMAX)
          IF (D.GT.180.0) D=ABS(D-360.0)
          IF (D.GT.PSIVAR/2.0) THEN
            PSI(K)=PSI(K)+180.0
            IF (PSI(K).GE.360.0) PSI(K)=PSI(K)-360.0
          ENDIF
40      CONTINUE
        DO 50 K=IFILM(I),I2
          D=ABS(PSI(K)-JMAX)
          IF (D.GT.180.0) D=ABS(D-360.0)
          IF (D.GT.PSIVAR/2.0) THEN
            PSIAVE=0.0
            N=0
            DO 60 J=-5,5
              IF ((K+J.GE.IFILM(I)).AND.(K+J
     +          .LE.I2)) THEN
                D=ABS(PSI(K+J)-JMAX)
                IF (D.GT.180.0) D=ABS(D-360.0)
                IF (D.LE.PSIVAR/2.0) THEN
                  PSIAVE=PSIAVE+PSI(K+J)
                  N=N+1
                ENDIF
              ENDIF
60          CONTINUE
            IF (N.NE.0) THEN
              PSIAVE=PSIAVE/N
              PSI(K)=PSIAVE
            ELSE
              PSI(K)=JMAX
            ENDIF
          ENDIF
50      CONTINUE
        IF (AVPRES.LT.PFLIP) THEN
          WRITE(*,7013)
7013      FORMAT('   --> Polarity will be flipped')
          DO 70 K=IFILM(I),I2
            PSI(K)=PSI(K)+180.0
            IF (PSI(K).GE.360.0) PSI(K)=PSI(K)-360.0
70        CONTINUE
        ENDIF
10    CONTINUE
C
C  Write output files
C
      DO 80 K=1,NH
        WRITE(11,8000)LINE(K)(1:SLEN2(LINE(K)))
8000    FORMAT(A)
80    CONTINUE
      DO 90 K=1,NP
        WRITE(11,7011)ILIST(K),PSI(K),
     +      THETA(K),PHI(K),SHX(K),
     +      SHY(K),ABSMAGP(K),FILM(K),
     +      DFMID1(K),DFMID2(K),ANGAST(K),
     +      OCC(K),LGP(K),S(K),PRES(K),DPRES(K)
90    CONTINUE
      DO 100 K=NH+1,NH+NF
        WRITE(11,8000)LINE(K)(1:SLEN2(LINE(K)))
100   CONTINUE
C
      CLOSE(11)
C
      WRITE(*,*)
      WRITE(*,*) ' Normal termination of set_polarity'
C
      END
C
C**************************************************************************
      INTEGER FUNCTION SLEN2(STRG)
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
      SLEN2=I
C
      RETURN
      END
