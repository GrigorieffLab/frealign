      PROGRAM BOOTVAR
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE
C
      INTEGER NN1,I,J,ID,IP,LN,SLEN3,SLEN2,K,NP,NX
      INTEGER NSAM,N1,N2,N3,ND,MODE,IFIRST,ILAST
      INTEGER LL,MM,NN,JC,NSAMH,IS,M,N,L,STAR,ST
      INTEGER IERR,HKL,II
      INTEGER,ALLOCATABLE :: R(:)
      REAL PSIZE,ANP,F,DF,P,DP,AM,AS,BM,BS,PI,PERR
      PARAMETER (PI=3.1415926535897)
      REAL,ALLOCATABLE :: DATA(:)
      DOUBLE PRECISION AVE,SIG
      DOUBLE PRECISION,ALLOCATABLE :: SUM(:)
      DOUBLE PRECISION,ALLOCATABLE :: SUM2(:)
      DOUBLE PRECISION,ALLOCATABLE :: SUM3(:)
      CHARACTER FNAME*200,CFORM,VX*15,ASYM*3
      CHARACTER ONAME*200,CN*10,FNAME2*200,Z*10
      CHARACTER FNAME3*200
      LOGICAL EX
      TYPE(C_PTR) FFTW_PLANS(2)
      DATA Z/'0000000000'/
C
      DATA  VX/'1.00 - 28.07.15'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C
      WRITE(*,1010)VX
1010  FORMAT(/'  BOOTVAR ',A15,/,
     +       /'  Copyright 2013 Howard Hughes Medical Institute.',
     +       /'  All rights reserved.',
     +       /'  Use is subject to Janelia Farm Research Campus',
     +        ' Software Copyright 1.1',
     +       /'  license terms ( http://license.janelia.org/license',
     +        '/jfrc_copyright_1_1.html )'/)
C
      WRITE(*,*)'Input seed for 3D maps?'
      READ(*,1000)FNAME
      WRITE(*,1000)FNAME
1000  FORMAT(A200)
      ST=STAR(FNAME)
      LN=SLEN3(FNAME)
      FNAME=FNAME(1:LN-4-ST)//FNAME(LN-3:LN)
      LN=LN-ST
C
      WRITE(*,*)'Input first and last #'
      READ(*,*)IFIRST,ILAST
      WRITE(*,*)IFIRST,ILAST
      WRITE(CN,7010)IFIRST
7010  FORMAT(I10)
      IF (ST.GT.10-SLEN2(CN)+1)
     + CN=CN(1:10-ST)//Z(10-ST+1:SLEN2(CN)-1)
     +    //CN(SLEN2(CN):10)
      ONAME=FNAME(1:LN-4)//CN(SLEN2(CN):10)//FNAME(LN-3:LN)
      WRITE(*,*)'First 3D volume:'
      WRITE(*,*)ONAME(1:SLEN3(ONAME))
C
      WRITE(*,*)'How many particles used for bootstrap?'
      READ(*,*)NP
      WRITE(*,*)NP
      ANP=SQRT(REAL(NP))
C
      WRITE(*,*)'Output sigF, sigP (0=no/1=yes)?'
      READ(*,*)HKL
      WRITE(*,*)HKL
C
      CALL GUESSF(ONAME,CFORM,EX)
      IF (.NOT.EX) THEN
        WRITE(*,1001) FNAME
1001    FORMAT(' File not found ',A200)
        GOTO 999
      ENDIF
C
      IF (HKL.NE.1) THEN
        WRITE(*,*)'Output average map?'
        READ(*,1000)FNAME2
        WRITE(*,1000)FNAME2
        WRITE(*,*)'Output sqrt(variance map)?'
        READ(*,1000)FNAME3
        WRITE(*,1000)FNAME3
      ELSE
        WRITE(*,*)'Output HKL list?'
        READ(*,1000)FNAME2
        WRITE(*,1000)FNAME2
        ANP=1.0
      ENDIF
C
      DO 99 I=IFIRST,ILAST
C
        WRITE(CN,7010)I
        IF (ST.GT.10-SLEN2(CN)+1)
     +    CN=CN(1:10-ST)//Z(10-ST+1:SLEN2(CN)-1)
     +       //CN(SLEN2(CN):10)
        ONAME=FNAME(1:LN-4)//CN(SLEN2(CN):10)//FNAME(LN-3:LN)
        CALL FLUSH(6)
        CALL IOPEN(ONAME,10,CFORM,N1,N2,N3,'OLD',ASYM,
     +             PSIZE,VX)
        IF ((N1.NE.N2).OR.(N1.NE.N3).OR.(N2.NE.N3)) THEN
          WRITE(*,*) 'X,Y,Z dimensions not equal'
          CALL ICLOSE(10)
          GOTO 999
        ENDIF
C
        IF (I.EQ.IFIRST) THEN
          NSAM=N1
          NSAMH=NSAM/2
          JC=NSAMH+1
C
          NN1=1
          IF (HKL.EQ.1) NN1=NSAM*NSAM*NSAM+2*NSAM*NSAM
C
          ALLOCATE(R(NSAM*NSAM*NSAM+2*NSAM*NSAM),
     +       SUM(NSAM*NSAM*NSAM+2*NSAM*NSAM),
     +       DATA(NSAM*NSAM*NSAM+2*NSAM*NSAM),
     +       SUM2(NSAM*NSAM*NSAM+2*NSAM*NSAM),
     +       SUM3(NN1),STAT=IERR)
          IF (IERR.NE.0) THEN
            WRITE(*,*) ' ERROR: Memory allocation failed'
            CALL ICLOSE(10)
            GOTO 999
          ENDIF
C
          CALL FFTW_PLANS_3D(NSAM,DATA,DATA,FFTW_PLANS(1),
     +                       FFTW_PLANS(2))
C
          DO 20 J=1,NSAM*NSAM*(NSAM+2)
            SUM(J)=0.0
            SUM2(J)=0.0
20        CONTINUE
        ELSEIF (N1.NE.NSAM) THEN
          WRITE(*,*)'N = ',I
          WRITE(*,*)'Volume dimensions not all the same'
          GOTO 999
        ENDIF
C
        CALL FLUSH(6)
        DO 10 J=1,NSAM
          IP=NSAM*(J-1)
          DO 10 K=1,NSAM
            ID=1+(NSAM+2)*((K-1)+NSAM*(J-1))
            CALL IREAD(10,DATA(ID),K+IP)
10      CONTINUE
C
        CALL ICLOSE(10)
C
C********************************************
        goto 666
        AVE=0.0
        SIG=0.0
        DO 29 J=0,NSAM*NSAM-1
          IS=J*(NSAM+2)
          DO 29 K=1,NSAM
            ID=IS+K
            AVE=AVE+DATA(ID)
            SIG=SIG+DATA(ID)**2
29      CONTINUE
        AVE=AVE/NSAM/NSAM/NSAM
        SIG=SIG/NSAM/NSAM/NSAM-AVE**2
        IF (SIG.GE.0.0) THEN
          SIG=SQRT(SIG)
        ELSE
          WRITE(*,*)'Error: Negative variance'
          GOTO 999
        ENDIF
C
        DO 28 J=0,NSAM*NSAM-1
          IS=J*(NSAM+2)
          DO 28 K=1,NSAM
            ID=IS+K
cc            DATA(ID)=DATA(ID)-AVE
            DATA(ID)=(DATA(ID)-AVE)/SIG
28      CONTINUE
666     continue
C********************************************
C
        IF (HKL.EQ.1) THEN
          CALL FFTW_FWD(DATA,DATA,FFTW_PLANS(1))
          DO 61 J=1,NSAM*NSAM*JC
              ID=2*J
              F=DATA(ID-1)**2+DATA(ID)**2
              SUM(J)=SUM(J)+SQRT(F)
              SUM2(J)=SUM2(J)+F
61        CONTINUE
          DO 62 J=1,NSAM*NSAM*(NSAM+2)
            SUM3(J)=SUM3(J)+DATA(J)
62        CONTINUE
        ELSE
          DO 60 J=1,NSAM*NSAM*(NSAM+2)
            SUM(J)=SUM(J)+DATA(J)
            SUM2(J)=SUM2(J)+DATA(J)**2
60        CONTINUE
        ENDIF
C
99    CONTINUE
C
      IF (HKL.NE.1) THEN
C
      DO 38 J=1,NSAM*NSAM*(NSAM+2)
        SUM(J)=SUM(J)/(ILAST-IFIRST+1)
        SUM2(J)=SUM2(J)/(ILAST-IFIRST+1)
38    CONTINUE
C
      DO 39 J=1,NSAM*NSAM*(NSAM+2)
        DATA(J)=ANP*SQRT(ABS(SUM2(J)-SUM(J)**2))
39    CONTINUE
C
      CALL IOPEN(FNAME3,10,CFORM,N1,N2,N3,'NEW',ASYM,
     +           PSIZE,VX)
C
      CALL FLUSH(6)
      DO 50 J=1,NSAM
        IP=NSAM*(J-1)
        DO 50 I=1,NSAM
          ID=1+(NSAM+2)*((I-1)+NSAM*(J-1))
          CALL IWRITE(10,DATA(ID),I+IP)
50    CONTINUE
C
      CALL ICLOSE(10)
C
      DO 40 J=1,NSAM*NSAM*(NSAM+2)
        DATA(J)=SUM(J)
40    CONTINUE
C
      CALL IOPEN(FNAME2,10,CFORM,N1,N2,N3,'NEW',ASYM,
     +           PSIZE,VX)
C
      CALL FLUSH(6)
      DO 51 J=1,NSAM
        IP=NSAM*(J-1)
        DO 51 I=1,NSAM
          ID=1+(NSAM+2)*((I-1)+NSAM*(J-1))
          CALL IWRITE(10,DATA(ID),I+IP)
51    CONTINUE
C
      CALL ICLOSE(10)
C
      ELSE
C
      DO 37 J=1,NSAM*NSAM*JC
        SUM(J)=SUM(J)/(ILAST-IFIRST+1)
        SUM2(J)=SUM2(J)/(ILAST-IFIRST+1)
37    CONTINUE
C
      INQUIRE(FILE=FNAME2,EXIST=EX)
      IF (EX) THEN
        OPEN(10,FILE=FNAME2,STATUS='OLD')
        CLOSE(10,STATUS='DELETE')
      ENDIF
      OPEN(10,FILE=FNAME2,STATUS='NEW')
C
      WRITE(10,1030)
1030  FORMAT('     H     K     L               F',
     +       '            SIGF               P',
     +       '            SIGP')
C     +       '            SIGP            RE_F',
C     +       '        SIG_RE_F            IM_F',
C     +       '        SIG_IM_F')
C
      ANP=ANP**2
      DO 88 L=1,JC
        LL=L-1
        DO 88 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          DO 88 N=1,NSAM
            NN=N-1
            IF (NN.GE.JC) NN=NN-NSAM
              IS=L+JC*((M-1)+NSAM*(N-1))
              ID=2*IS
              AM=SUM3(ID-1)
              BM=SUM3(ID)
C              AS=ANP*(ABS(SUM2(ID-1)-SUM(ID-1)**2))
C              BS=ANP*(ABS(SUM2(ID)-SUM(ID)**2))
              F=SUM(IS)
              DF=SUM2(IS)-SUM(IS)**2
C              IF (F.NE.0.0) THEN
C                DF=AS*(AM/F)**2+BS*(BM/F)**2
C                DF=SQRT(DF)*2.0
C              ELSE
C                DF=0.0
C              ENDIF
              IF (AM.EQ.0.0) THEN
                P=0.0
              ELSEIF (AM.GT.0.0) THEN
                P=ATAN(BM/AM)
              ELSE
                P=ATAN(BM/AM)+PI
              ENDIF
              IF (P.GE.0.0) THEN
                P=P/PI*180.0
              ELSE
                P=P/PI*180.0+360.0
              ENDIF
C              DP=AS*(BM/F**2)**2+BS*(AM/F**2)**2
C              DP=SQRT(DP)/PI*180.0
              IF ((LL.EQ.0).AND.(MM.EQ.0)
     +          .AND.(NN.EQ.0)) THEN
                WRITE(10,1020) LL,MM,NN,F,0.0,0.0,
     +                       0.0
C     +                       0.0,0.0,0.0,0.0,0.0
              ELSE
                WRITE(10,1020) LL,MM,NN,F,DF,P,
     +            PERR(F,DF)/PI*180.0
C     +            PERR(F,DF)/PI*180.0,AM,
C     +            SQRT(AS),BM,SQRT(BS)
              ENDIF
1020          FORMAT(3I6,4F16.3)
88    CONTINUE
C
      CLOSE(10)
C
      ENDIF
C
      WRITE(*,*)
      WRITE(*,*)' NORMAL TERMINATION OF BOOTVAR'
C
999   CONTINUE
C
      END
C
C**************************************************************************
      REAL FUNCTION PERR(F,DF)
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER I,N
      PARAMETER (N=90)
      REAL F,DF,PI,R,PM,X,Y,A
      PARAMETER (PI=3.1415926535897)
C
      PM=0.0
      R=PI/N*0.99
      DO 10 I=0,N
        X=F+DF*COS(I*R)
        Y=DF*SIN(I*R)
        IF (ABS(X).LE.ABS(Y)) THEN
          A=PI/2.0-ATAN2(X,Y)
        ELSE
          A=ATAN2(Y,X)
        ENDIF
        PM=PM+ABS(A)
10    CONTINUE
C
      PERR=PM/(N+1)
C
      RETURN
      END
C**************************************************************************
      INTEGER FUNCTION STAR(STRG)
C**************************************************************************
C Determines if STRG contains '*'
C**************************************************************************
      IMPLICIT NONE
C       
      INTEGER I
      CHARACTER*200 STRG
C**************************************************************************
C
      STAR=0
      DO 10 I=1,200,1
        IF (STRG(I:I).EQ."*") STAR=STAR+1
10    CONTINUE
C
      RETURN
      END
C**************************************************************************
      INTEGER FUNCTION SLEN3(STRG)
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
      SLEN3=I
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
