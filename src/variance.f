C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE VARIANCE(NSAM,SPEC,SHX,SHY,RI,OUTD,OUTC,
     +   DATD,DATC,IC,VSN,VN,VVSN,VVN,AMAGP,BUF,NA,FFTW_PLANS)
C**************************************************************************
C  Calculates variances of signal+background and background alone
C  Used in LMAIN.
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE

      INTEGER NSAM,JC,L,LL,M,MM,ID,I,NMID,NA,II,NSAMH
      INTEGER IRAD,IC(*)
      REAL PHASE,OUTD(*),FI1,CRAD2,AMAGP,FRACD,FRACO
      REAL DATD(*),PI,BUF(*)
      REAL RI,SHX,SHY,RIP
      DOUBLEPRECISION VSN(*),VN(*),DSUM,DSUM2,BSUM,BSUM2
      DOUBLEPRECISION VVSN(*),VVN(*)
      PARAMETER  (PI=3.1415926535897)
      COMPLEX PSHFT,SPEC(*),OUTC(*),DATC(*)
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
      RIP=RI**2
      JC=NSAM/2+1
      NSAMH=NSAM/2
C
      IF (RIP.LT.NSAMH**2) THEN
C
      NMID=NSAM/2+1
C
      DO 111 L=1,JC
        LL=L-1   
        DO 111 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          PHASE=SHX*LL+SHY*MM
          PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
            ID=L+JC*(M-1)
            DATC(ID)=SPEC(ID)*PSHFT
111   CONTINUE 
C
      CALL FFTW_BWD(DATD,DATD,FFTW_PLANS(2))
C
C     MASK PARTICLE, SHARP EDGE
      DSUM=0.0D0
      DSUM2=0.0D0
      BSUM=0.0D0
      BSUM2=0.0D0
      NA=0
      DO 18 I=1,2*NSAM
        BUF(I)=0.0
18    CONTINUE
      DO 17 I=1,NSAM
        IC(I)=0.0
17    CONTINUE
      DO 21 I=1,NSAM
        FI1=REAL(I-NMID)**2
        DO 21 II=1,NSAM 
          ID=II+(NSAM+2)*(I-1)
          CRAD2=REAL(II-NMID)**2+FI1
          DATD(ID)=DATD(ID)/NSAM/NSAMH
          OUTD(ID)=DATD(ID)
          IF (CRAD2.GT.RIP) THEN
            BSUM=BSUM+OUTD(ID)
            BSUM2=BSUM2+OUTD(ID)**2
          ELSE
            NA=NA+1
            DSUM=DSUM+DATD(ID)
            DSUM2=DSUM2+DATD(ID)**2
          ENDIF
21    CONTINUE
      DSUM=DSUM/NA
      DSUM2=DSUM2/NA
      DSUM2=SQRT(DSUM2-DSUM**2)
      BSUM=BSUM/(NSAM*NSAM-NA)
      BSUM2=BSUM2/(NSAM*NSAM-NA)
      BSUM2=SQRT(BSUM2-BSUM**2)
C
      DO 19 I=1,NSAM
        FI1=REAL(I-NMID)**2
        DO 19 II=1,NSAM
          ID=II+(NSAM+2)*(I-1)
          CRAD2=REAL(II-NMID)**2+FI1
          IF (CRAD2.GT.RIP) THEN
            DATD(ID)=0.0
            OUTD(ID)=(OUTD(ID)-BSUM)/BSUM2
          ELSE
            DATD(ID)=(DATD(ID)-BSUM)/BSUM2
            OUTD(ID)=0.0
          ENDIF
19    CONTINUE
C
      CALL FFTW_FWD(DATD,DATD,FFTW_PLANS(1))
      CALL FFTW_FWD(OUTD,OUTD,FFTW_PLANS(1))
C
C     FRAC will scale the variances to a value which would have
C     been obtained from an unmasked image
      FRACD=REAL(NA)/NSAM/NSAM
      FRACO=REAL(NSAM*NSAM-NA)/NSAM/NSAM
      DO 20 L=1,NSAM/2
        LL=L-1
        DO 20 M=1,NSAM  
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          ID=L+JC*(M-1)
          IRAD=INT(SQRT(REAL(LL**2+MM**2))+0.5)+1
          IF (FRACO.NE.0) THEN
            BUF(IRAD)=BUF(IRAD)+CABS(OUTC(ID))**2/FRACO
          ENDIF
          IF (FRACD.NE.0) THEN
            BUF(NSAM+IRAD)
     +        =BUF(NSAM+IRAD)+CABS(DATC(ID))**2/FRACD
          ENDIF
          IC(IRAD)=IC(IRAD)+1
20    CONTINUE
C
      DO 22 I=1,NSAM
        IF (IC(I).NE.0.0) THEN
          BUF(I)=BUF(I)/IC(I)
          BUF(NSAM+I)=BUF(NSAM+I)/IC(I)
        ENDIF
        VN(I)=VN(I)+BUF(I)
        VVN(I)=VVN(I)+BUF(I)**2
        VSN(I)=VSN(I)+BUF(NSAM+I)
        VVSN(I)=VVSN(I)+BUF(NSAM+I)**2
22    CONTINUE
C
      ELSE
C
      DO 23 I=1,NSAM
        IC(I)=PI*(I-1)
23    CONTINUE
C
      ENDIF
C
      RETURN
      END
