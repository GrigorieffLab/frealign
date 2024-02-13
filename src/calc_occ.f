C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      PROGRAM CALC_OCC
C**************************************************************************
C Program to calculate occupancies from input parameter files based on
C LogP values. Generates new parameter files with updated OCC values.
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER MAXP,MAXC,I,J,K,N,SLEN2,IERR,MAXH,LGP,II
      INTEGER,ALLOCATABLE :: ILIST(:,:),FILM(:,:),NH(:)
      INTEGER,ALLOCATABLE :: NF(:),NP(:),ABSMAGP(:,:)
      REAL,ALLOCATABLE :: PSI(:,:),THETA(:,:),PHI(:,:)
      REAL,ALLOCATABLE :: SHX(:,:),SHY(:,:),SIG(:,:)
      REAL,ALLOCATABLE :: DFMID1(:,:),DFMID2(:,:)
      REAL,ALLOCATABLE :: PSHFT(:,:)
      REAL,ALLOCATABLE :: OCC(:,:),ALGP(:,:),PRES(:,:)
      REAL,ALLOCATABLE :: DPRES(:,:),ANGAST(:,:),A(:)
      REAL AMAX,SUMP,F,ON,PSSNR,SUMA
      CHARACTER*200,ALLOCATABLE :: LINE(:,:)
      CHARACTER*47,ALLOCATABLE :: FC(:)
      CHARACTER*51,ALLOCATABLE :: LC(:)
      CHARACTER FPAR*200,LINEC,VX*15
      LOGICAL EX,LPSSNR,FREALIGNX
C**************************************************************************
C
      DATA  VX/'1.03 - 26.04.17'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C
      WRITE(*,1010)VX
1010  FORMAT(/'  CALC_OCC ',A15,/,
     +       /'  Copyright 2013 Howard Hughes Medical Institute.',
     +       /'  All rights reserved.',
     +       /'  Use is subject to Janelia Farm Research Campus',
     +        ' Software Copyright 1.1',
     +       /'  license terms ( http://license.janelia.org/license',
     +        '/jfrc_copyright_1_1.html )'/)
C
C  Read input files
C
      WRITE(*,2001)
2001  FORMAT(' No. of input files?')
      READ(*,*)N
      WRITE(*,3000)N
3000  FORMAT(I6,/)
      IF (N.EQ.0) THEN
        WRITE(*,*)
        WRITE(*,*) 'N = 0: Nothing to do. Exiting...'
        GOTO 9999
      ENDIF
      MAXC=N
      LPSSNR=.FALSE.
      FREALIGNX=.FALSE.
C
      WRITE(*,2002)
2002  FORMAT(' OCC change multiplier?')
      READ(*,*)F
      WRITE(*,3001)F
3001  FORMAT(F6.2,/)
C
      DO 10 I=1,N
C
        WRITE(*,4000)I
4000    FORMAT('*************************************************',
     +      //,' Input file name for file ',I6)
        READ(*,7007)FPAR
        WRITE(*,8000)FPAR(1:SLEN2(FPAR))
        OPEN(10,FILE=FPAR,STATUS='OLD')
C
        IF (I.EQ.1) THEN
          MAXP=1
          MAXH=1
7006      CONTINUE
          READ(10,8000,END=7005)LINEC
          IF (LINEC.NE.'C') THEN
            MAXP=MAXP+1
          ELSE
            MAXH=MAXH+1
          ENDIF
          GOTO 7006
C
7005      CONTINUE
          REWIND(UNIT=10)
C
          ALLOCATE(ILIST(MAXP,MAXC),FILM(MAXP,MAXC),
     +       NF(MAXC),NP(MAXC),NH(MAXC),PSI(MAXP,MAXC),
     +       THETA(MAXP,MAXC),PHI(MAXP,MAXC),SIG(MAXP,MAXC),
     +       SHX(MAXP,MAXC),SHY(MAXP,MAXC),ABSMAGP(MAXP,MAXC),
     +       DFMID1(MAXP,MAXC),DFMID2(MAXP,MAXC),A(MAXC),
     +       OCC(MAXP,MAXC),ALGP(MAXP,MAXC),PRES(MAXP,MAXC),
     +       DPRES(MAXP,MAXC),ANGAST(MAXP,MAXC),FC(MAXC),
     +       LINE(2*MAXH,MAXC),LC(MAXC),PSHFT(MAXP,MAXC),
     +       STAT=IERR)
          IF (IERR.NE.0) THEN
            STOP ' ERROR: Memory allocation failed'
          ENDIF
        ENDIF
C
        NH(I)=1
        NF(I)=1
        NP(I)=1
        A(I)=0.0
7009    CONTINUE
        READ(10,7007,END=7008)LINE(NF(I),I)
7007    FORMAT(A200)
        IF (LINE(NF(I),I)(1:1).EQ.'C') THEN
          IF (((INDEX(LINE(NF(I),I),'PRES').NE.0).AND.
     +       (INDEX(LINE(NF(I),I),'DPRES').NE.0)).OR.
     +       (INDEX(LINE(NF(I),I),'Phase res. / B factor').NE.0).OR.
     +       (INDEX(LINE(NF(I),I),'Phase residual target').NE.0).OR.
     +       (INDEX(LINE(NF(I),I),'Phase residual threshold').NE.0))
     +    STOP ' ****ERROR**** Old par file. Convert to new format.'
          NF(I)=NF(I)+1
          GOTO 7009
        ELSEIF (LEN_TRIM(LINE(NF(I),I)).LT.120) THEN
          STOP ' ****ERROR**** Old par file. Convert to new format.'
        ELSEIF ((LEN_TRIM(LINE(NF(I),I)).GT.140)
     +           .AND.(.NOT.FREALIGNX)) THEN
          WRITE(*,*) ' FrealignX parameter file'
          WRITE(*,*)
          FREALIGNX=.TRUE.
        ENDIF
        NH(I)=NF(I)-1
        IF (FREALIGNX) THEN
          READ(LINE(NF(I),I),7013)ILIST(NP(I),I),PSI(NP(I),I),
     +      THETA(NP(I),I),PHI(NP(I),I),SHX(NP(I),I),
     +      SHY(NP(I),I),ABSMAGP(NP(I),I),FILM(NP(I),I),
     +      DFMID1(NP(I),I),DFMID2(NP(I),I),ANGAST(NP(I),I),
     +      PSHFT(NP(I),I),OCC(NP(I),I),LGP,SIG(NP(I),I),
     +      PRES(NP(I),I),DPRES(NP(I),I)
7013      FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,3F8.2,I10,F11.4,2F8.2)
        ELSE
          READ(LINE(NF(I),I),7011)ILIST(NP(I),I),PSI(NP(I),I),
     +      THETA(NP(I),I),PHI(NP(I),I),SHX(NP(I),I),
     +      SHY(NP(I),I),ABSMAGP(NP(I),I),FILM(NP(I),I),
     +      DFMID1(NP(I),I),DFMID2(NP(I),I),ANGAST(NP(I),I),
     +      OCC(NP(I),I),LGP,SIG(NP(I),I),PRES(NP(I),I),
     +      DPRES(NP(I),I)
7011      FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I10,F11.4,2F8.2)
          PSHFT(NP(I),I)=0.0
        ENDIF
        ALGP(NP(I),I)=LGP
        IF(ILIST(NP(I),I).LE.0) THEN
          WRITE(*,*) ' ***WARNING*** blank line in parameter file'
          GO TO 7009
        ENDIF
        A(I)=A(I)+OCC(NP(I),I)
        NP(I)=NP(I)+1
        GOTO 7009
C
7008    CONTINUE
        NF(I)=NF(I)-NH(I)-1
        NP(I)=NP(I)-1
        A(I)=A(I)/NP(I)
        WRITE(*,1000)NH(I),NF(I),NP(I),A(I)
1000    FORMAT(' No. of header lines read = ',I10,/,
     +         ' No. of footer lines read = ',I10,/,
     +         ' No. of particles read    = ',I10,/,
     +         ' Average occupancy in %   = ',F10.4,/)
C
        CLOSE(10)
10    CONTINUE
C
C
C  Check par file consistency
C
      DO 20 I=2,N
        IF(NP(I).NE.NP(1))
     +    STOP ' ****ERROR**** No. of particles varies'
20    CONTINUE
C
      DO 30 K=1,NP(1)
        DO 30 I=2,N
          IF(ILIST(K,I).NE.ILIST(K,1)) THEN
            WRITE(*,5000)I,K
5000        FORMAT(' In files 1 and ',I6,', particle ',I10)
            STOP ' ****ERROR**** Particle sequences differ'
          ENDIF
          IF(FILM(K,I).NE.FILM(K,1)) THEN
            WRITE(*,5000)I,K
            STOP ' ****ERROR**** Film no. differs'
          ENDIF
30    CONTINUE
C
C
C  Calculate occupancies
C
      DO 80 K=1,NP(1)
        SUMP=0.0
        DO 120 I=1,N
          SUMP=SUMP+OCC(K,I)
120     CONTINUE
        IF (SUMP.EQ.0.0) THEN
          SUMP=100.0
          DO 140 I=1,N
            OCC(K,I)=100.0/REAL(N)
140       CONTINUE
        ENDIF
        DO 130 I=1,N
          OCC(K,I)=OCC(K,I)/SUMP*100.0
130     CONTINUE
        AMAX=-1.0E30
        SUMP=0.0
        DO 90 I=1,N
          AMAX=MAX(ALGP(K,I),AMAX)
90      CONTINUE
        DO 100 I=1,N
          IF (AMAX-ALGP(K,I).LT.10.0) THEN
            SUMP=SUMP+EXP(ALGP(K,I)-AMAX)*A(I)
          ENDIF
100     CONTINUE
        DO 110 I=1,N
          IF (AMAX-ALGP(K,I).LT.10.0) THEN
            ON=EXP(ALGP(K,I)-AMAX)*A(I)/SUMP*100.0
          ELSE
            ON=0.0
          ENDIF
          OCC(K,I)=F*(ON-OCC(K,I))+OCC(K,I)
110     CONTINUE
        SUMP=0.0
        DO 150 I=1,N
          SUMP=SUMP+SIG(K,I)*OCC(K,I)/100.0
150     CONTINUE
        DO 160 I=1,N
          SIG(K,I)=SUMP
160     CONTINUE
80    CONTINUE
C
C
C  Write output files
C
      WRITE(*,6000)
6000  FORMAT('*************************************************',
     +     /,'*************************************************',/)
      DO 40 I=1,N
C
        WRITE(*,7000)I
7000    FORMAT(' Output file name for file ',I6)
        READ(*,7007)FPAR
        WRITE(*,8000)FPAR(1:SLEN2(FPAR))
        INQUIRE(FILE=FPAR,EXIST=EX)
        IF (EX) THEN
          OPEN(10,FILE=FPAR,STATUS='OLD')
          CLOSE(10,STATUS='DELETE')
        ENDIF
        OPEN(10,FILE=FPAR,STATUS='NEW')
        DO 50 K=1,NH(I)
          WRITE(10,8000)LINE(K,I)(1:SLEN2(LINE(K,I)))
8000      FORMAT(A)
50      CONTINUE
        DO 60 K=1,NP(I)
          IF (FREALIGNX) THEN
            WRITE(10,7013)ILIST(K,I),PSI(K,I),
     +        THETA(K,I),PHI(K,I),SHX(K,I),
     +        SHY(K,I),ABSMAGP(K,I),FILM(K,I),
     +        DFMID1(K,I),DFMID2(K,I),ANGAST(K,I),
     +        PSHFT(K,I),OCC(K,I),NINT(ALGP(K,I)),
     +        SIG(K,I),PRES(K,I),DPRES(K,I)
          ELSE
            WRITE(10,7011)ILIST(K,I),PSI(K,I),
     +        THETA(K,I),PHI(K,I),SHX(K,I),
     +        SHY(K,I),ABSMAGP(K,I),FILM(K,I),
     +        DFMID1(K,I),DFMID2(K,I),ANGAST(K,I),
     +        OCC(K,I),NINT(ALGP(K,I)),SIG(K,I),
     +        PRES(K,I),DPRES(K,I)
          ENDIF
60      CONTINUE
        DO 70 K=NH(I)+1,NH(I)+NF(I)
          IF (.NOT.LPSSNR) THEN
            WRITE(10,8000)LINE(K,I)(1:SLEN2(LINE(K,I)))
            DO 170 J=1,200-8
              IF (LINE(K,I)(J:J+8).EQ.'Part_SSNR') THEN
                LPSSNR=.TRUE.
                GOTO 180
              ENDIF
170         CONTINUE
          ENDIF
          IF (LPSSNR) THEN
200         CONTINUE
            DO 190 J=1,200-6
              IF (LINE(K,I)(J:J+6).EQ.'Average') THEN
                WRITE(10,8000)LINE(K,I)(1:SLEN2(LINE(K,I)))
                LPSSNR=.FALSE.
                GOTO 180
              ENDIF
190         CONTINUE
            SUMP=0.0
            SUMA=0.0
            DO 210 II=1,N
              READ(LINE(K,II),7012)FC(II),PSSNR,LC(II)
7012          FORMAT(A47,F11.4,A51)
              SUMP=SUMP+PSSNR**2*A(II)
              SUMA=SUMA+A(II)
210         CONTINUE
            SUMP=SUMP/SUMA
            WRITE(10,7012)FC(I),SQRT(SUMP),LC(I)
          ENDIF
180       CONTINUE
70      CONTINUE
C
        CLOSE(10)
40    CONTINUE
C
9999  CONTINUE
C
      WRITE(*,*)
      WRITE(*,*) ' Normal termination of calc_occ'
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
