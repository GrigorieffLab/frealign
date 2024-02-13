C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C ---------------------------------------------------------------------------------        
      SUBROUTINE FIND_FPART(STD,VTD,JC,NSAM,NSAMH,
     +			A3DF,B3DF,C3DF,D3DF,S3DF,V3DF,
     +                  IREDUN,ASYM,RI,RIC,B3DV,C3DV,
     +                  PSIZE,FPART,FSC,MW,DALT,SSNR,
     +                  FFTW_PLANS)
C ---------------------------------------------------------------------------------        
C       Find fparticle to scale SSNR for Wiener filter
C       Used in Frealign.
C ---------------------------------------------------------------------------------        
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE
C
      INTEGER L,LL,M,MM,N,NN,JC,NSAM,NSAMH,IS,ID,J
      INTEGER NSAM3,IREDUN,ISUM,K,NF,I
      REAL PSHFTR,T,MAVE,MSTD,RI,RIC,FMASK,DF,EP,FRAD
      PARAMETER  (T=1.5,DF=0.2,EP=0.005,FRAD=30.0)
      DOUBLE PRECISION STD,VTD,SUM1,SUM2
      REAL S3DF(*),B3DV(*),C3DV(*),FPART,FPARTMAX
      REAL V3DF(*),CC,CCMAX,FPART_CC,DF1
      REAL PSIZE,FSC(*),SSNR(*),FRAD2,TT,MW,DALT
      COMPLEX A3DF(*),B3DF(*)
      COMPLEX C3DF(*),D3DF(*)
      CHARACTER ASYM*3
      TYPE(C_PTR) FFTW_PLANS(*)
C ---------------------------------------------------------------------------------        
      IF(STD.EQ.0.0.OR.VTD.EQ.0.0) THEN
        WRITE(*,6100) STD,VTD
6100    FORMAT(' Abnormal termination, probably no accept. particles'/
     .         '   either STD or VTD zero, STD,VTD =',2F15.7/
     .         '   at least two particles needed')
        STOP 'STOP 6100'
      ENDIF

      IF (MW.EQ.0.0) THEN

      NSAM3=NSAM*NSAM*NSAM
      FRAD2=2.0*PSIZE/FRAD*NSAM

      DO 60 N=1,NSAM
        NN=N-1
        IF (NN.GE.JC) NN=NN-NSAM
        DO 60 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          DO 60 L=1,JC
            LL=L-1
            ISUM=(LL+MM+NN)
            PSHFTR=1.0
            IF (MOD(ISUM,2).NE.0) PSHFTR=-1.0
              ID=L+JC*((M-1)+NSAM*(N-1))
              B3DF(ID)=A3DF(ID)*PSHFTR/(S3DF(ID)+STD*IREDUN)/NSAM3
              C3DF(ID)=(A3DF(ID)+D3DF(ID))*PSHFTR
     .          /((S3DF(ID)+V3DF(ID))+(STD+VTD)*IREDUN)/NSAM3
              C3DF(ID)=C3DF(ID)*EXP(-(LL**2+MM**2+NN**2)/FRAD2)
60    CONTINUE
C
      CALL FFTW_BWD(B3DF,B3DF,FFTW_PLANS(4))
      CALL FFTW_BWD(C3DF,C3DF,FFTW_PLANS(4))
      IF (ASYM(1:1).EQ.'H') THEN
        CALL MASK3D_C(NSAM,C3DV,RI,RIC,1.0,FMASK,
     +                MAVE,MSTD)
      ELSE
        CALL MASK3D(NSAM,C3DV,RI,RIC,1.0,FMASK,
     +              MAVE,MSTD)
      ENDIF
C
      DO 80 K=0,INT(T/0.1)
        TT=MAVE+(T-K*0.1)*MSTD
        ISUM=0
        DO 90 J=0,NSAM*NSAM-1
          IS=J*(NSAM+2)
          DO 90 I=1,NSAM
            ID=IS+I
            IF (C3DV(ID).GE.TT) ISUM=ISUM+1
90      CONTINUE
        IF (ISUM.GT.0.2*NSAM3*FMASK) GOTO 99
80    CONTINUE
C
99    CONTINUE
C
      ISUM=0
      SUM1=0.0
      SUM2=0.0
      DO 70 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 70 I=1,NSAM
          ID=IS+I
          IF (C3DV(ID).LT.TT) THEN
            B3DV(ID)=0.0
          ELSE
            SUM1=SUM1+B3DV(ID)
            SUM2=SUM2+B3DV(ID)**2
            ISUM=ISUM+1
          ENDIF
70    CONTINUE
      IF (ISUM.NE.0) THEN
        SUM1=SUM1/ISUM
        SUM2=SUM2/ISUM
      ENDIF
      MAVE=SUM1  
      MSTD=SQRT(ABS(SUM2-SUM1**2))
C
      DF1=-DF
      FPART=0.1
      CALL VOLMEASURE(NSAM,C3DV,FPART)
      CCMAX=FPART_CC(JC,NSAM,NSAMH,
     +               B3DF,C3DF,D3DF,S3DF,V3DF,
     +               ASYM,RI,RIC,B3DV,C3DV,
     +               FPART,SSNR,MAVE,MSTD,FSC,
     +               FFTW_PLANS)
10    CONTINUE
      FPART=FPART+DF1
      CC=FPART_CC(JC,NSAM,NSAMH,
     +            B3DF,C3DF,D3DF,S3DF,V3DF,
     +            ASYM,RI,RIC,B3DV,C3DV,
     +            FPART,SSNR,MAVE,MSTD,FSC,
     +            FFTW_PLANS)
      IF (CC.GT.CCMAX) THEN   
        CCMAX=CC
        GOTO 10
      ELSE
        FPART=FPART-DF1
        IF (ABS(DF1).GT.EP) THEN
          DF1=-DF1/3.0
          GOTO 10
        ENDIF
      ENDIF
C
      ELSE
C
      FPART=1000.0*MW/DALT/(PSIZE*NSAM)**3
C
      ENDIF
C
      FPART=ABS(FPART)
      IF (FPART.NE.0.0) THEN
        DO 20 I=1,JC
          SSNR(I)=SSNR(I)/FPART
20      CONTINUE
      ENDIF
C
      RETURN
      END
