C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE SIGMA(NSAM,IRAD,AMAG,A3DF,DATC,
     +                 PHI,THETA,PSI,SHX,SHY,SINCLUT,IPAD,
     +                 IEWALD,THETATR,CTFF,RBUF,RI2,
     +                 MAXR1,MAXR2,SIGNSET,SIGMA2N,SIGMA2,
     +                 LOGPC,DMASK,ASYM,K,SIGP,PWEIGHTS,
     +                 RCLAS,FFTW_PLANS,SM)
C**************************************************************************
C  Calculates real space correlation coefficient.
C  Uses functions PEXTRACT, FFTW.
C  Used in LMAIN.
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE
C
      INTEGER NSAM,L,M,LL,MM,JC,ID,IPAD,IRAD,IB1,ID2,ITEMP
      INTEGER M2,IEWALD,NSAMH,NSAM2,N1,MAXR1,MAXR2
      INTEGER MAXR12,MAXR22,N2,K,IBIN,MRCLAS
      REAL SHX,SHY,PHASE,SINCLUT(*),DMASK(*),DM(9),SIGNSET
      REAL PHI,THETA,PSI,AMAG,THETATR,RBUF(*),SIGMA2N,S2OLD
      REAL LOGPC,RI2,RI2P,SIGMA2,XM,YM,X2,Y2,RI2F,PI
      REAL SIGMA2OLD,SIGP,PWEIGHTS(*),RCLAS,SM(9)
      DOUBLE PRECISION SUM1,SUM2,AVE1,AVE2,CCPART1,SUM3,AVE3
      DOUBLE PRECISION SUM4,AVE4,SUM5,AVE5,CCPART2,A,S2
      PARAMETER  (PI=3.1415926535897)
      COMPLEX PSHFT,A3DF(*),CTFF(*)
      COMPLEX DATC(*),Z
      CHARACTER ASYM*3
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
      NSAM2=NSAM*NSAM
      IB1=NSAM*(NSAM+2)
      NSAMH=NSAM/2
      JC=NSAM/2+1
      IF (ASYM(1:1).EQ.'H') THEN
        RI2P=(NSAMH-1)**2
      ELSE
        RI2P=RI2/AMAG/AMAG
      ENDIF
      RI2F=(DMASK(4)/AMAG)**2
      MRCLAS=INT(NSAM*RCLAS*ABS(AMAG))
C      MAXR12=1
C      MAXR22=NSAM**2
C  Applying a bandpass filter to the image and reference
C  violates the assumption that the real-space noise is
C  white. However, it appears to help in the classification
C  of noisy data.
C
      MAXR12=MAXR1**2
      MAXR22=MRCLAS**2
C
      CALL PEXTRACT(NSAM,IRAD,RBUF,A3DF,
     +         PHI,THETA,PSI,SINCLUT,IPAD,
     +         THETATR,IEWALD,CTFF,SM)
      DO 222 L=1,JC
        LL=L-1
        DO 222 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          ITEMP=LL**2+MM**2
          IF ((ITEMP.GE.MAXR12).AND.(ITEMP.LT.MAXR22)) THEN
            IBIN=INT(SQRT(FLOAT(ITEMP))+0.5)+1
            PHASE=-SHX*LL-SHY*MM
            PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
              ID=L+JC*(M-1)
              Z=PSHFT*DATC(ID)
C              Z=PSHFT*DATC(ID)*PWEIGHTS(IBIN)
              ID2=2*ID+IB1
              RBUF(ID2-1)=REAL(Z)
              RBUF(ID2)=AIMAG(Z)
          ELSE
              ID=L+JC*(M-1)
              ID2=2*ID+IB1
              RBUF(ID2-1)=0.0
              RBUF(ID2)=0.0
              ID2=2*ID
              RBUF(ID2-1)=0.0
              RBUF(ID2)=0.0
          ENDIF
222   CONTINUE
C      RBUF(IB1+1)=0.0
C      RBUF(IB1+2)=0.0
C      RBUF(1)=0.0
C      RBUF(2)=0.0
      CALL PHASEFLIP(NSAM,RBUF(IB1+1),CTFF)
      CALL PHASEFLIP(NSAM,RBUF,CTFF)
      CALL FFTW_BWD(RBUF(IB1+1),RBUF(IB1+1),FFTW_PLANS(2))
      CALL FFTW_BWD(RBUF,RBUF,FFTW_PLANS(2))
C
      CALL ROTMAT(-PSI,-THETA,-PHI,AMAG,DM)
      CALL MATMUL_T(DM,SM,DM)
      XM=DM(1)*(DMASK(1)-JC)+DM(4)*(DMASK(2)-JC)
     +  +DM(7)*(DMASK(3)-JC)+1.0
      YM=DM(2)*(DMASK(1)-JC)+DM(5)*(DMASK(2)-JC)
     +  +DM(8)*(DMASK(3)-JC)+1.0
C
      SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      SUM4=0.0
      SUM5=0.0
      CCPART1=0.0
      CCPART2=0.0
      AVE1=0.0
      AVE2=0.0
      AVE3=0.0
      AVE4=0.0
      AVE5=0.0
      N1=0
      N2=0
      DO 223 L=1,NSAM
        LL=(L-1)**2
        X2=(L-XM)**2
        IF (L.GE.NSAMH) THEN
          LL=(L-1-NSAM)**2
          X2=(L-XM-NSAM)**2
        ENDIF
        DO 223 M=1,NSAM
          MM=(M-1)**2
          Y2=(M-YM)**2
          IF (M.GE.NSAMH) THEN
            MM=(M-1-NSAM)**2
            Y2=(M-YM-NSAM)**2
          ENDIF
          ID=L+(NSAM+2)*(M-1)
          ID2=IB1+ID
          IF (X2+Y2.LE.RI2F) THEN
C          IF (LL+MM.LE.RI2P) THEN
            CCPART1=CCPART1+RBUF(ID)*RBUF(ID2)
            SUM1=SUM1+RBUF(ID)**2
            SUM2=SUM2+RBUF(ID2)**2
            AVE1=AVE1+RBUF(ID)
            AVE2=AVE2+RBUF(ID2)
            N1=N1+1
          ENDIF
C          IF (X2+Y2.LE.RI2F) THEN
C            CCPART2=CCPART2+RBUF(ID)*RBUF(ID2)
C            SUM3=SUM3+RBUF(ID)**2
C            SUM4=SUM4+RBUF(ID2)**2
C            AVE3=AVE3+RBUF(ID)
C            AVE4=AVE4+RBUF(ID2)
C            N2=N2+1
C          ENDIF
          SUM5=SUM5+RBUF(ID2)**2
          AVE5=AVE5+RBUF(ID2)
223   CONTINUE
C
C      CALL IOPEN("testimage.mrc",66,"M",NSAM,NSAM,1,'NEW',
C     +           "C1 ",2.44,"123456789012345")
C      L=0
C      DO 50 M=1,NSAM
C        L=L+1
C        ID=1+(NSAM+2)*(M-1)
C        CALL IWRITE(66,RBUF(ID),L)
C50    CONTINUE
C      CALL ICLOSE(66)
C
C  SUM1 = sum over reference projection densities
C  SUM2 = sum over particle image densities
      CCPART1=CCPART1-AVE1*AVE2/N1
      SUM1=SUM1-AVE1**2/N1
      SUM2=SUM2-AVE2**2/N1
C      CCPART2=CCPART2-AVE3*AVE4/N2
C      SUM3=SUM3-AVE3**2/N2
C      SUM4=SUM4-AVE4**2/N2
      SUM5=SUM5-AVE5**2/NSAM2
      A=ABS(CCPART1/SUM1)
C
C
C  S2      = sigma^2 of noise in particle image
C            (scale is off due to unnormalized FFTs)
C  SIGMA2  = sigma^2 of noise in particle image
C            (scaled correctly)
C  SIGMA2N = sigma^2 of noise in particle image
C            divided by sqrt of variance in the
C            image and reference. This will make
C            SIGMA2N appropriate for weighting the
C            prior based on parameter distributions
      S2=ABS(SUM2-A**2*SUM1)/N1
ccc      SIGMA2=S2/NSAM2/NSAM2*4.0
      SIGMA2=S2/NSAM2/NSAM2
C
C Let's initially assume that we have mostly noise in our
C images. SUM5 is the total variance of our particle image
C after applying the bandpass filter.
      S2OLD=SUM5/NSAM2
      SIGMA2OLD=S2OLD/NSAM2/NSAM2
ccc      SIGMA2OLD=S2OLD/NSAM2/NSAM2*4.0
C For the current cycle, use old sigma if available
      IF (SIGNSET.NE.0.0) THEN
        IF ((SQRT(SIGMA2OLD)/SIGNSET.LE.3.0).AND.
     +      (SQRT(SIGMA2OLD)/SIGNSET.GE.0.33)) THEN
          S2OLD=SIGNSET**2*NSAM2*NSAM2
ccc          S2OLD=SIGNSET**2*NSAM2*NSAM2/4.0
          SIGMA2OLD=SIGNSET**2
        ELSE
C Something seems wrong with the previously recorded sigma
C because it deviates significantly from the new sigma.
C Let's discard it to keep likelihood calculation and
C weighting of priors approx. correct.
C      WRITE(*,*)' Resetting unrealistic sigma for particle ',K
        ENDIF
      ENDIF
      SIGMA2N=S2OLD/A/SQRT(SUM1*SUM2)
C
C     Now calculate likelihood (priors are added later)
      N1=N1*SIGMA2OLD/SIGP**2
      LOGPC=-N1*SIGMA2/SIGMA2OLD/2.0
     +      -N1*LOG(2.0*PI*SIGMA2OLD)/2.0
C
      RETURN
      END
