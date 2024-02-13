C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PSEARCH(NSAM,IRAD,AMAG,SPEC,A3DF,DANGIN,
     +			 MAXR1,MAXR2,PHI,THETA,PSI,DSHX,DSHY,
     +			 PRESA,SANG,NSANG,CCD,CCC,CBUF,
     +			 IQUADMAX,MASK,RBFACT,SINCLUT,IPAD,
     +                   RBUF,IEWALD,THETATR,CTFF,
     +                   RI2,RI3,RIH,HALFW,TESTPAR,IPMAX,ASYM,
     +                   IRAN,FFTW_PLANS,SM)
C**************************************************************************
C  Searches for 5 parameters describing orientation and position using 
C  cross-correlation method, but also calculates average phase residual.
C  Calls CCP and uses function PRES.
C  Used in LMAIN.
C  MW added IEWALD, THETATR for Ewald sphere correction (used in CCP)
C**************************************************************************
C
      USE ISO_C_BINDING
C
      IMPLICIT NONE

      INTEGER NSAM,MAXR1,MAXR2,IRAD,NSANG,I,I3,IQUAD,J
      INTEGER IQUADMAX,MASK(5),IPAD,IEWALD,IPMAX,K,KK
      INTEGER NSAMH,IC,JC,II,JJ,ITEMP,ID,IRAN
      REAL XPAR(5),SANG(*),CCMAX,PI,CCD(*),PRESA
      REAL PHI,THETA,PSI,DSHX,DSHY,AMAG,CC,RBUF,RANDOM
      REAL PHI1,THETA1,PSI1,RBFACT,SINCLUT(*),WEIGHT
      REAL RI2,RI3,RIH,HALFW,TESTPAR(6,*),SM(9)
      REAL DANGIN,PHIR,PSIR,THETAR,DANG
      PARAMETER  (PI=3.1415926535897)
      COMPLEX SPEC(*),CCC(*)
      COMPLEX A3DF(*),CBUF(*),CTFF(*)
      REAL THETATR
      CHARACTER ASYM*3
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
      CCMAX=-1E30
      DO 30 K=1,IPMAX
        TESTPAR(6,K)=0.0
30    CONTINUE
C
      NSAMH=NSAM/2
      JC=NSAM/2+1
      IC=NSAM*JC
      DO 31 I=0,NSAMH
          II=I+1
          DO 31 J=-NSAMH,NSAMH
            ITEMP=I**2+J**2
            JJ=J+1
            IF (JJ.LT.1) JJ=JJ+NSAM
            ID=II+JC*(JJ-1)
            WEIGHT=EXP(-RBFACT*ITEMP)
            CBUF(ID)=SPEC(ID)*WEIGHT
31    CONTINUE
C
      DANG=DANGIN/180.0*PI
      PHIR=(RANDOM(IRAN)-0.5)*DANG*REAL(MASK(1))
      THETAR=(RANDOM(IRAN)-0.5)*DANG*REAL(MASK(2))
      PSIR=(RANDOM(IRAN)-0.5)*DANG*REAL(MASK(3))
C
      DANG=DANGIN/5.0
      IF (DANG.LT.1.0) DANG=1.0
C
      DO 10 I=1,NSANG
      	I3=3*I
      	XPAR(1)=MASK(1)*(SANG(I3-2)+PHIR)+PHI
        IF (ASYM(1:1).EQ.'H') THEN
      	  XPAR(2)=MASK(2)*(SANG(I3-1)-THETA)+THETA
        ELSE
          XPAR(2)=MASK(2)*(SANG(I3-1)+THETAR)+THETA
        ENDIF
      	XPAR(3)=MASK(3)*(SANG(I3)+PSIR)+PSI
      	IF (XPAR(1).GT.PI) XPAR(1)=XPAR(1)-2.0*PI
      	IF (XPAR(2).GT.PI) XPAR(2)=XPAR(2)-2.0*PI
      	IF (XPAR(3).GT.PI) XPAR(3)=XPAR(3)-2.0*PI
      	XPAR(4)=0.0
      	XPAR(5)=0.0
      	CALL CCP(NSAM,IRAD,AMAG,CBUF,A3DF,RBUF,DANG,
     +		 MAXR1,MAXR2,XPAR(1),XPAR(2),XPAR(3),XPAR(4),XPAR(5),
     +		 CC,CCD,CCC,CBUF(IC+1),IQUAD,MASK,SINCLUT,
     +           IPAD,IEWALD,THETATR,CTFF,RI2,RI3,RIH,HALFW,
     +           IQUADMAX,FFTW_PLANS,SM)
      	IF (XPAR(3).GT.PI) XPAR(3)=XPAR(3)-2.0*PI
        DO 40 K=1,IPMAX
          IF (CC.GT.TESTPAR(6,K)) THEN
            DO 50 KK=IPMAX,K+1,-1
              DO 60 J=1,6
                TESTPAR(J,KK)=TESTPAR(J,KK-1)
60            CONTINUE
50          CONTINUE
            DO 20 J=1,5
              TESTPAR(J,K)=XPAR(J)
20          CONTINUE
            TESTPAR(6,K)=CC
            GOTO 41
          ENDIF
40      CONTINUE
C
41      CONTINUE
      	IF (CC.GT.CCMAX) THEN
      	  CCMAX=CC
      	  PHI1=XPAR(1)
      	  THETA1=XPAR(2)
      	  PSI1=XPAR(3)
      	  DSHX=XPAR(4)
      	  DSHY=XPAR(5)
C		on entry to PSEARCH, DSHX and DSHY are normally zero
      	ENDIF
10    CONTINUE
      PHI=PHI1
      THETA=THETA1
      PSI=PSI1
C      write(*,11) CCMAX,PRESA*180.0/PI
11    format(' Exiting  PSEARCH with CCMAX, PRES',2F8.3)
      RETURN
      END
