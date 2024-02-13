C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      PROGRAM SELECT_CLASSES
C**************************************************************************
C Program to copy select classes into a new particle image stack
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER N1,N2,N3,SLEN2,MAXP,MAXH,K,I,IERR,LEN,NH
      INTEGER NF,NP,IP,ID,CNT,N,M(4),SETOCC,L
      REAL PSIZE,R(12)
      INTEGER,ALLOCATABLE :: ILIST(:),FILM(:),LGP(:)
      INTEGER,ALLOCATABLE :: ABSMAGP(:),SEL(:),IC(:)
      REAL,ALLOCATABLE :: PSI(:),THETA(:),PHI(:),S(:)
      REAL,ALLOCATABLE :: SHX(:),SHY(:),OCC(:),DATA(:)
      REAL,ALLOCATABLE :: DFMID1(:),DFMID2(:),MAXO(:)
      REAL,ALLOCATABLE :: DPRES(:),ANGAST(:),PRES(:)
      CHARACTER*200,ALLOCATABLE :: LINE(:),PARIN(:)
      CHARACTER*200,ALLOCATABLE :: PAROUT(:)
      CHARACTER*200 STACKIN,STACKOUT
      CHARACTER VX*15,CFORM*1,ASYM*3,C1*1
      CHARACTER CDATE*8,CTIME*10,CZONE*5
      LOGICAL EX
C**************************************************************************
C
      DATA VX/'1.00 - 26.05.15'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C
      WRITE(*,1010)VX
1010  FORMAT(/'  SELECT_CLASSES ',A15,/,
     +       /'  Copyright 2013 Howard Hughes Medical Institute.',
     +       /'  All rights reserved.',
     +       /'  Use is subject to Janelia Farm Research Campus',
     +        ' Software Copyright 1.1',
     +       /'  license terms ( http://license.janelia.org/license',
     +        '/jfrc_copyright_1_1.html )'/)
C
C  Read input file
C
      WRITE(*,*)' Intput stack?'
      READ(*,1020)STACKIN
1020  FORMAT(A200)
      WRITE(*,1000)STACKIN(1:SLEN2(STACKIN))
1000  FORMAT(A)
      WRITE(*,*)' Number of classes to select from?'
      READ(*,*)N
      WRITE(*,*)N
      IF (N.EQ.0) THEN
        WRITE(*,*)
        WRITE(*,*) 'N = 0: Nothing to do. Exiting...'
        GOTO 999
      ENDIF
      ALLOCATE(PARIN(N),PAROUT(N),SEL(N),STAT=IERR)
      IF (IERR.NE.0) THEN
        STOP ' ERROR: Memory allocation failed'
      ENDIF
      WRITE(*,4010)
4010  FORMAT(/'*************************************************')
      DO 20 I=1,N
        WRITE(*,4000)I
4000    FORMAT(/' Input parameter file ',I6)
        READ(*,1020)PARIN(I)
        WRITE(*,1000)PARIN(I)(1:SLEN2(PARIN(I)))
        WRITE(*,*)' Select this class (0 = NO, 1 = YES)?'
        READ(*,*)SEL(I)
        WRITE(*,*)SEL(I)
        IF (SEL(I).EQ.1) THEN
          WRITE(*,4001)I
4001      FORMAT(' Output parameter file ',I6)
          READ(*,1020)PAROUT(I)
          WRITE(*,1000)PAROUT(I)(1:SLEN2(PAROUT(I)))
        ENDIF
20    CONTINUE
      WRITE(*,4010)
      WRITE(*,*)' Output stack?'
      READ(*,1020)STACKOUT
      WRITE(*,1000)STACKOUT(1:SLEN2(STACKOUT))
      WRITE(*,*)' Reset occupancies to 100 (0 = NO, 1 = YES)?'
      READ(*,*)SETOCC
      WRITE(*,*)SETOCC
C
      CALL GUESSF(STACKIN,CFORM,EX)
      IF (EX) THEN
        CALL IOPEN(STACKIN,10,CFORM,N1,N2,N3,'OLD',ASYM,PSIZE,VX)
      ELSE
        WRITE(*,*)' ERROR: File does not exist'
        GOTO 999
      ENDIF
C
      OPEN(11,FILE=PARIN(N),STATUS='OLD')
C
      MAXP=1
      MAXH=1
7006  CONTINUE
      READ(11,8000,END=7005)C1
8000  FORMAT(A)
      IF (C1.NE.'C') THEN
        MAXP=MAXP+1
      ELSE
        MAXH=MAXH+1
      ENDIF
      GOTO 7006
C
7005  CONTINUE
      REWIND(UNIT=11)
C
      ALLOCATE(ILIST(MAXP),FILM(MAXP),PSI(MAXP),
     +       THETA(MAXP),PHI(MAXP),OCC(MAXP),
     +       SHX(MAXP),SHY(MAXP),ABSMAGP(MAXP),
     +       DFMID1(MAXP),DFMID2(MAXP),PRES(MAXP),
     +       DPRES(MAXP),ANGAST(MAXP),IC(MAXP),
     +       LINE(2*MAXH),LGP(MAXP),S(MAXP),
     +       DATA(N1*N2),MAXO(MAXP),STAT=IERR)
      IF (IERR.NE.0) THEN
        WRITE(*,*) ' ERROR: Memory allocation failed'
        CALL ICLOSE(10)
        CLOSE(11)
        GOTO 999
      ENDIF
C
      DO 120 K=N,1,-1
C
      NP=1
      IF (K.NE.N) OPEN(11,FILE=PARIN(K),STATUS='OLD')
7016  CONTINUE
      READ(11,7007,END=7017)LINE(1)
      IF (LINE(1)(1:1).EQ.'C') GOTO 7016
      IF (SLEN2(LINE(1)).LE.95) THEN
        GOTO 199
      ELSEIF (SLEN2(LINE(1)).LE.128) THEN
        READ(LINE(1),7015,ERR=199,IOSTAT=CNT)M(1),
     +    (R(I),I=1,5),M(2),M(3),(R(I),I=6,9),M(4),
     +    R(11),R(12)
        IF (CNT.NE.0) GOTO 199
        R(10)=1.0
        GOTO 198
      ELSE
        READ(LINE(1),7011,ERR=199,IOSTAT=CNT)M(1),
     +    (R(I),I=1,5),M(2),M(3),(R(I),I=6,9),M(4),
     +    (R(I),I=10,12)
        IF (CNT.NE.0) GOTO 199
        GOTO 198
      ENDIF
199   CONTINUE
      IF (NP.EQ.1) THEN
        WRITE(*,*) ' ERROR reading par file.',
     +             ' File format not compatible.'
        CALL ICLOSE(10)
        CLOSE(11)
        GOTO 999
      ELSE
        WRITE(*,*)
     +    ' ERROR: Something wrong in parameter file.'
        WRITE(*,*) LINE(NF)
        CALL ICLOSE(10)
        CLOSE(11)
        GOTO 999
      ENDIF
198   CONTINUE
C
      IF(M(1).LE.0) THEN
        WRITE(*,*) ' ***WARNING*** blank line in parameter file'
        GO TO 7016
      ENDIF
C
      IF (K.EQ.N) THEN
        MAXO(NP)=R(9)
        IC(NP)=K
      ELSEIF (R(9).GE.MAXO(NP)) THEN
        MAXO(NP)=R(9)
        IC(NP)=K
      ENDIF
C
      NP=NP+1
      GOTO 7016
C
7017  CONTINUE
C
      CLOSE(11)
C
120   CONTINUE
C
      DO 70 K=1,N
C
      IF (SEL(K).EQ.1) THEN
C
C      IF (K.NE.N) OPEN(11,FILE=PARIN(K),STATUS='OLD')
      OPEN(11,FILE=PARIN(K),STATUS='OLD')
      NH=1
      NF=1
      NP=1
      L=1
7009  CONTINUE
      READ(11,7007,END=7008)LINE(NF)
7007  FORMAT(A200)
      IF (LINE(NF)(1:1).EQ.'C') THEN
        NF=NF+1
        GOTO 7009
      ENDIF
      NH=NF-1
      IF (SLEN2(LINE(NF)).LE.95) THEN
        GOTO 99
      ELSEIF (SLEN2(LINE(NF)).LE.128) THEN
        READ(LINE(NF),7015,ERR=99,IOSTAT=CNT)M(1),
     +    (R(I),I=1,5),M(2),M(3),(R(I),I=6,9),M(4),
     +    R(11),R(12)
7015    FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I13,2F8.2)
        IF (CNT.NE.0) GOTO 99
        R(10)=1.0
        GOTO 98
      ELSE
        READ(LINE(NF),7011,ERR=99,IOSTAT=CNT)M(1),
     +    (R(I),I=1,5),M(2),M(3),(R(I),I=6,9),M(4),
     +    (R(I),I=10,12)
7011    FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I10,F11.4,2F8.2)
        IF (CNT.NE.0) GOTO 99
        GOTO 98
      ENDIF
99    CONTINUE
      IF (NP.EQ.1) THEN
        WRITE(*,*) ' ERROR reading par file.',
     +             ' File format not compatible.'
        CALL ICLOSE(10)
        CLOSE(11)
        GOTO 999
      ELSE
        WRITE(*,*)
     +    ' ERROR: Something wrong in parameter file.'
        WRITE(*,*) LINE(NF)
        CALL ICLOSE(10)
        CLOSE(11)
        GOTO 999
      ENDIF
98    CONTINUE
C
      IF(M(1).LE.0) THEN
C        WRITE(*,*) ' ***WARNING*** blank line in parameter file'
        GO TO 7009
      ENDIF
C
      IF (SEL(IC(L)).NE.0) THEN
        ILIST(NP)=M(1)
        PSI(NP)=R(1)
        THETA(NP)=R(2)
        PHI(NP)=R(3)
        SHX(NP)=R(4)
        SHY(NP)=R(5)
        ABSMAGP(NP)=M(2)
        FILM(NP)=M(3)
        DFMID1(NP)=R(6)
        DFMID2(NP)=R(7)
        ANGAST(NP)=R(8)
        OCC(NP)=R(9)
        LGP(NP)=M(4)
        S(NP)=R(10)
        PRES(NP)=R(11)
        DPRES(NP)=R(12)
C
        IF (SETOCC.EQ.1) THEN
          IF (IC(L).EQ.K) THEN
            OCC(NP)=100.0
          ELSE
            OCC(NP)=0.0
          ENDIF
        ENDIF
C
        NP=NP+1
C
      ENDIF
C
      L=L+1
      GOTO 7009
C
7008  CONTINUE
C
      CLOSE(11)
C
      NP=NP-1
      NF=NF-NH-1
      L=L-1
      WRITE(*,1030)K,NH,NF,L,NP
1030  FORMAT(/' Parameter file for class      ',I10,/,
     +        ' No. of header lines read    = ',I10,/,
     +        ' No. of footer lines read    = ',I10,/,
     +        ' No. of particles read       = ',I10,/,
     +        ' No. of particles for output = ',I10)
C
      IF (L.NE.N3) THEN
        WRITE(*,*) ' ERROR: Particle number in stack and'//
     +             ' par file differ'
        CALL ICLOSE(10)
        GOTO 999
      ENDIF
      IF (NP.EQ.0) THEN
        WRITE(*,*) ' No particles selected. Stopping...'
        GOTO 999
      ENDIF
C
      INQUIRE(FILE=PAROUT(K),EXIST=EX)
      IF (EX) THEN
        OPEN(13,FILE=PAROUT(K),STATUS='OLD')
        CLOSE(13,STATUS='DELETE')
      ENDIF
      OPEN(13,FILE=PAROUT(K),STATUS='NEW')
C
      DO 80 I=1,NH
        WRITE(13,8000)LINE(I)(1:SLEN2(LINE(I)))
80    CONTINUE
      N=0
      DO 90 I=1,NP
        N=N+1
        WRITE(13,7011)N,PSI(I),
     +    THETA(I),PHI(I),SHX(I),
     +    SHY(I),ABSMAGP(I),FILM(I),
     +    DFMID1(I),DFMID2(I),ANGAST(I),
     +    OCC(I),LGP(I),S(I),PRES(I),DPRES(I)
90    CONTINUE
      DO 100 I=NH+1,NH+NF
        WRITE(13,8000)LINE(I)(1:SLEN2(LINE(I)))
100   CONTINUE
      CLOSE(13)
C
      ENDIF
C
70    CONTINUE
C
      INQUIRE(FILE=STACKOUT,EXIST=EX)
      IF (EX) THEN
        OPEN(12,FILE=STACKOUT,STATUS='OLD')
        CLOSE(12,STATUS='DELETE')
      ENDIF
      ASYM='C1 '
      CALL IOPEN(STACKOUT,12,CFORM,N1,N2,NP,'NEW',ASYM,PSIZE,VX)
C
      N=0
      DO 110 K=1,L
        IF (SEL(IC(K)).NE.0) THEN
          N=N+1
          IP=N2*(K-1)
          DO 19 I=1,N1
            ID=1+N1*(I-1)
            CALL IREAD(10,DATA(ID),I+IP)
19        CONTINUE
          IP=N2*(N-1)
          DO 30 I=1,N1
            ID=1+N1*(I-1)
            CALL IWRITE(12,DATA(ID),I+IP)
30        CONTINUE
        ENDIF
110   CONTINUE
      CALL ICLOSE(10)
      CALL ICLOSE(12)
C
999   CONTINUE
C
      END
