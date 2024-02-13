C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      PROGRAM APPLY_MASK
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE
C
      INTEGER N1,N2,N3,I,J,ID,IP,MODE,NR1,NR2
      INTEGER NSAM,K,N,NN,ITEMP,IM,IS
      INTEGER JC,L,LL,M,MM,IERR,SLEN2,IRAD
      REAL PSIZE,HW,RI2,RI3,RIH,RII,RAD2,EDGE
      REAL MASKR,PI,WEIGHT,SUM,MULT,P,SCAL
      REAL,ALLOCATABLE :: DATA(:),MASK(:),MDATA(:)
      PARAMETER (PI=3.1415926535897)
      CHARACTER FNAME*200,CFORM,VX*15,TITLE*1600,ASYM*3
      TYPE(C_PTR) FFTW_PLANS(2)
C
      DATA  VX/'1.00 - 23.07.15'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C
      WRITE(*,1010)VX
1010  FORMAT(/'  APPLY_MASK ',A15,/,
     +       /'  Copyright 2013 Howard Hughes Medical Institute.',
     +       /'  All rights reserved.',
     +       /'  Use is subject to Janelia Farm Research Campus',
     +        ' Software Copyright 1.1',
     +       /'  license terms ( http://license.janelia.org/license',
     +        '/jfrc_copyright_1_1.html )'/)
C
      WRITE(*,*)'Image format [M,S,I]?'
      READ(*,1100)CFORM
1100  FORMAT(A1)
      WRITE(*,11100)CFORM
11100 FORMAT(3X,A1)
      WRITE(*,*)
      WRITE(*,*)'Input 3D map?'
      READ(*,1000)FNAME
1000  FORMAT(A200)
      WRITE(*,11000) FNAME(1:SLEN2(FNAME))
11000 FORMAT(3X,A)
C
      CALL IOPEN(FNAME,10,CFORM,N1,N2,N3,'OLD',ASYM,
     +           PSIZE,VX)
      IF ((N1.NE.N2).OR.(N1.NE.N3).OR.(N2.NE.N3)) THEN
        WRITE(*,*) 'ERROR: X,Y,Z dimensions not equal'
        GOTO 999
      ENDIF
      IF (MOD(N1,2).NE.0) THEN
        WRITE(*,*) 'ERROR: Dimensions must be even'
        GOTO 999
      ENDIF
      ALLOCATE(DATA(N1*N1*N1+2*N1*N1),
     +  MASK(N1*N1*N1+2*N1*N1),
     +  MDATA(N1*N1*N1+2*N1*N1),STAT=IERR)
      IF (IERR.NE.0) THEN
        WRITE(*,*) ' ERROR: Memory allocation failed'
        GOTO 999
      ENDIF
C
      NSAM=N1
      SCAL=1.0/NSAM/NSAM/NSAM
      JC=NSAM/2+1
C
      CALL FFTW_PLANS_3D(NSAM,DATA,DATA,FFTW_PLANS(1),FFTW_PLANS(2))
C
      DO 10 K=1,N3
        IP=NSAM*(K-1)
        DO 10 I=1,NSAM
          ID=1+(NSAM+2)*((I-1)+NSAM*(K-1))
          CALL IREAD(10,DATA(ID),I+IP)
10    CONTINUE
C
      CALL ICLOSE(10)
C
      WRITE(*,*)'Pixel size from header =',PSIZE
      WRITE(*,*)'Pixel size in A?'
      WRITE(*,*)'(Enter * to use header value)'
      READ(*,1000)FNAME
      IF (FNAME(1:1).NE."*") THEN
        READ(FNAME,*)PSIZE
      ENDIF
      WRITE(*,*)PSIZE
C
      WRITE(*,*)
      WRITE(*,*)'Input 3D mask?'
      READ(*,1000)FNAME
      WRITE(*,11000) FNAME(1:SLEN2(FNAME))
C
      CALL IOPEN(FNAME,10,CFORM,N1,N2,N3,'OLD',ASYM,
     +           P,VX)
      IF ((N1.NE.NSAM).OR.(N1.NE.N2).OR.(N1.NE.N3)) THEN
        WRITE(*,*)
     +    'ERROR: Mask and map X,Y,Z dimensions differ'
        GOTO 999
      ENDIF
C
      DO 13 K=1,NSAM
        IP=NSAM*(K-1)
        DO 13 I=1,NSAM
          ID=1+(NSAM+2)*((I-1)+NSAM*(K-1))
          CALL IREAD(10,MASK(ID),I+IP)
13    CONTINUE
C
      CALL ICLOSE(10)
C
      WRITE(*,*)'Width of cosine edge to add (in pixel)?'
      READ(*,*)IRAD
      WRITE(*,*)IRAD
      IRAD=ABS(IRAD)
C
      IF (IRAD.NE.0) THEN
C
      WRITE(*,*)
      WRITE(*,*)'Smoothing edge...'
C     Threshold mask

      DO 86 J=1,NSAM*NSAM
        IS=(J-1)*(NSAM+2)
        DO 86 I=1,NSAM
          ID=IS+I
          IF (MASK(ID).GT.0.0) THEN
            MASK(ID)=1.0  
          ELSE
            MASK(ID)=0.0
          ENDIF
86    CONTINUE
C
C     Smoothing the edge.....
      NN=(NSAM+2)*NSAM*NSAM
      DO 11 I=1,NSAM
        DO 11 J=1,NSAM
          DO 11 K=1,NSAM
            ID=K+(NSAM+2)*((J-1)+NSAM*(I-1))
            IF (MASK(ID).EQ.1.0) THEN
              DO 12 L=-IRAD,IRAD
                DO 12 M=-IRAD,IRAD
                  DO 12 N=-IRAD,IRAD
                    RAD2=L**2+M**2+N**2
                    RAD2=SQRT(RAD2)
                    EDGE=(1.0+COS(PI*RAD2/IRAD))/2.0
                    ID=K+N+(NSAM+2)*((J+M-1)+NSAM*(I+L-1))
                    IF ((ID.GE.1).AND.(ID.LE.NN)
     +                  .AND.(RAD2.LE.IRAD)) THEN
                      IF (MASK(ID).LT.EDGE) MASK(ID)=EDGE
                    ENDIF
12            CONTINUE
            ENDIF
11    CONTINUE
C
      ENDIF
C
      SUM=0.0
      DO 40 J=1,NSAM*NSAM
        IS=(J-1)*(NSAM+2)
        DO 40 I=1,NSAM
          ID=IS+I
          SUM=SUM+MASK(ID)**2
          MDATA(ID)=DATA(ID)
40    CONTINUE
      WRITE(*,*)
      WRITE(*,1200) INT(SUM)
1200  FORMAT(' Number of voxels inside mask: ',I14)
C
      WRITE(*,*)
      WRITE(*,*)'Weight for density outside mask?'
      READ(*,*)MULT
      WRITE(*,*)MULT
C
      WRITE(*,*)
      WRITE(*,*)
     +  'Low-pass filter outside (0=no, 1=Gauss, 2=cosine edge)?'
      READ(*,*)IM
      WRITE(*,*)IM
C
      IF (IM.NE.0) THEN
C
      WRITE(*,*)
      IF (IM.EQ.1) THEN
        WRITE(*,*)'Gauss filter radius in A?'
      ELSE
        WRITE(*,*)'Cosine edge filter radius in A?'
      ENDIF
      READ(*,*)MASKR
      WRITE(*,*)MASKR
      IF (IM.NE.1) THEN
        WRITE(*,*)'Width of edge in pixels?'
        READ(*,*)HW
        WRITE(*,*)HW
      ENDIF
C
      MASKR=MASKR/PSIZE
      IF (MASKR.NE.0.0) THEN
        MASKR=1.0/MASKR*NSAM
      ELSE
        MASKR=1.0D30
      ENDIF
      RII=MASKR+HW/2.0
      RIH=MASKR-HW/2.0
      IF (RIH.LT.0.0) RIH=0.0
      RI2=RIH**2
      RI3=RII**2
C
      CALL FFTW_FWD(DATA,DATA,FFTW_PLANS(1))
C
      DO 111 L=1,JC
        LL=L-1
        DO 111 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          DO 111 N=1,NSAM
            NN=N-1
            IF (NN.GE.JC) NN=NN-NSAM
            ITEMP=LL**2+MM**2+NN**2
            IF (IM.EQ.1) THEN
              WEIGHT=EXP(-ITEMP/(2.0*MASKR**2))
            ELSEIF ((ITEMP.GT.RI2).AND.(ITEMP.LE.RI3)) THEN
              RAD2=SQRT(REAL(ITEMP))
              WEIGHT=(1.0+COS(PI*(RAD2-RIH)/HW))/2.0
            ELSEIF (ITEMP.GT.RI3) THEN
              WEIGHT=0.0
            ENDIF
              ID=2*(L+JC*((M-1)+NSAM*(N-1)))
              DATA(ID-1)=DATA(ID-1)*WEIGHT*SCAL*MULT
              DATA(ID)=DATA(ID)*WEIGHT*SCAL*MULT
111   CONTINUE
C
      CALL FFTW_BWD(DATA,DATA,FFTW_PLANS(2))
C
      ELSE
C
      DO 50 J=1,NSAM*NSAM
        IS=(J-1)*(NSAM+2)
        DO 50 I=1,NSAM
          ID=IS+I
          DATA(ID)=DATA(ID)*MULT
50    CONTINUE
C
      ENDIF
C
      DO 60 J=1,NSAM*NSAM
        IS=(J-1)*(NSAM+2)
        DO 60 I=1,NSAM
          ID=IS+I
          DATA(ID)=DATA(ID)*(1.0-MASK(ID))+MDATA(ID)*MASK(ID)
60    CONTINUE
C
      WRITE(*,*)
      WRITE(*,*)'Output masked 3D map?'
      READ(*,1000)FNAME
      WRITE(*,11000) FNAME(1:SLEN2(FNAME))
C
      TITLE='APPLY_MASK: Masked output with low-pass filter'
      CALL IOPEN(FNAME,10,CFORM,NSAM,NSAM,NSAM,'NEW',ASYM,
     +           PSIZE,VX)
C
      DO 20 K=1,NSAM
        IP=NSAM*(K-1)
        DO 20 I=1,NSAM
          ID=1+(NSAM+2)*((I-1)+NSAM*(K-1))
          CALL IWRITE(10,DATA(ID),I+IP)
20    CONTINUE
C
      CALL ICLOSE(10)
      WRITE(*,*) ' NORMAL TERMINATION OF APPLY_MASK'
C
999   CONTINUE
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
