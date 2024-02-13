C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PRESB(NSAM,IRAD,AMAG,SPEC,A3DF,
     +		       MAXR1,MAXR2,PHI,THETA,PSI,SHX,SHY,
     +                 RBIN,SINCLUT,IPAD,NS,BUFA,BUFB,BUFAB,
     +                 CS,WL,WGH1,WGH2,DF1,DF2,ANGAST,THET,
     +                 IEWALD,CTFF,RI2,RI3,RIH,HALFW,
     +                 ASYM,FFTW_PLANS,SM)
C**************************************************************************
C Calculates correlation coeff. in resolution bins (RBIN)
C Uses functions AINTERPO3DS and AINTERPO3DBIG.
C Used in LMAIN.
C MW uses function EWALDEX if IEWALD is true
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE

      INTEGER NSAM,JC,I,J,II,JJ,MAXR1,MAXR2,IRAD,NSAMH,ID,MAXR12
      INTEGER MAXR22,ITEMP,IBIN,IPAD,NN1,NS(*),IEWALD
      INTEGER IB1,IB,IC
      REAL DM(9),PHI,THETA,PSI,CPHI,SPHI,CTHE,STHE,CPSI,SPSI
      REAL PHASE,X3,Y3,Z3,SHX,SHY,SUMA,SUMB,CTF,BUFAB(*),SCAL
      REAL WGT,AMAG,RBIN(*),SINCLUT(*),BUFA(*),BUFB(*),SM(9)
      REAL CS,WL,WGH1,WGH2,DF1,DF2,ANGAST,THET,RI2,RI3,RIH,HALFW
      COMPLEX SAMP,SPEC(*),PSHFT,AINPO
      COMPLEX A3DF(*),AINTERPO3DBIG,AINTERPO3DS
      COMPLEX EWALDEX,CTFF(*),CTFR,CTFL
      CHARACTER ASYM*3
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
      JC=NSAM/2+1
      IC=NSAM*JC
      NSAMH=NSAM/2
C      MAXR12=MAXR1**2
C      MAXR22=MAXR2**2
      MAXR12=1
      MAXR22=NSAMH**2
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      CTHE=COS(THETA)
      STHE=SIN(THETA)
      CPSI=COS(PSI)
      SPSI=SIN(PSI)

      DM(1)=(CPHI*CTHE*CPSI-SPHI*SPSI)/ABS(AMAG)
      DM(2)=(SPHI*CTHE*CPSI+CPHI*SPSI)/ABS(AMAG)
      DM(3)=-STHE*CPSI/ABS(AMAG)
      DM(4)=(-CPHI*CTHE*SPSI-SPHI*CPSI)/ABS(AMAG)
      DM(5)=(-SPHI*CTHE*SPSI+CPHI*CPSI)/ABS(AMAG)
      DM(6)=STHE*SPSI/ABS(AMAG)
C MW next lines uncommented for Ewald sphere correction
      DM(7)=STHE*CPHI/ABS(AMAG)
      DM(8)=STHE*SPHI/ABS(AMAG)
      DM(9)=CTHE/ABS(AMAG)
      CALL MATMUL_T(DM,SM,DM)

      DO 10 I=1,NSAM
        NS(I)=0
        BUFA(I)=0.0
        BUFB(I)=0.0
        BUFAB(I)=0.0
10    CONTINUE
      IB1=NSAM
C      IB2=IB1+NSAM*NSAM
      SCAL=1.0/NSAM/NSAM
      DO 11 I=NSAM+1,NSAM*NSAM+3*NSAM
        BUFAB(I)=0.0
11    CONTINUE
C
      DO 30 I=0,JC-1
      	II=I+1
      	DO 30 J=-JC+1,JC-1
      	  ITEMP=I**2+J**2
          IF ((ITEMP.GE.MAXR12).AND.(ITEMP.LT.MAXR22)) THEN
      	    JJ=J+1
      	    IF (JJ.LT.1) JJ=JJ+NSAM
      	      ID=II+JC*(JJ-1)
      	      CTFR=CTFF(ID)
      	      CTFL=CTFF(ID+IC)
            IF (IEWALD.NE.0) THEN
              IF (IEWALD.LT.0) THEN
                CTFR=CONJG(CTFR)
                CTFL=CONJG(CTFL)
              ENDIF
C MW call ewald sphere corrected extraction
              AINPO=EWALDEX(NSAM,IRAD,A3DF,SINCLUT,IPAD,
     +                         I,J,DM,THET,CTFR,CTFL)
C              AINPO=AINPO*(CTFR+CONJG(CTFL))
            ELSE
              X3=DM(1)*I+DM(4)*J
              Y3=DM(2)*I+DM(5)*J
              Z3=DM(3)*I+DM(6)*J
              IF (IRAD.EQ.0) THEN
                AINPO=AINTERPO3DBIG(NSAM,IPAD,A3DF,X3,Y3,Z3)
              ELSE
                AINPO=AINTERPO3DS(IPAD*NSAM,IRAD,A3DF,
     +                              X3,Y3,Z3,SINCLUT,IPAD)
              ENDIF
C              AINPO=AINPO*CABS(CTFR+CONJG(CTFL))**2
              AINPO=AINPO*(CTFR+CONJG(CTFL))
            ENDIF
C      	    IF (II.NE.JC) THEN
              IB=IB1+2*ID
              BUFAB(IB-1)=REAL(AINPO)*SCAL
              BUFAB(IB)=AIMAG(AINPO)*SCAL
C      	    ELSE
C              IB=IB2+2*JJ
C              BUFAB(IB-1)=REAL(AINPO)*SCAL
C              BUFAB(IB)=AIMAG(AINPO)*SCAL
C      	    ENDIF
          ENDIF
30    CONTINUE
      CALL FFTW_BWD(BUFAB(IB1+1),BUFAB(IB1+1),FFTW_PLANS(2))
      IF (ASYM(1:1).EQ.'H')THEN
        CALL MASKCOS2_C(NSAM,BUFAB(IB1+1),RI2,RI3,RIH,
     +                HALFW,AMAG,PSI)
      ELSE
        CALL MASKCOS2(NSAM,BUFAB(IB1+1),RI2,RI3,RIH,
     +                HALFW,AMAG)
      ENDIF
      CALL FFTW_FWD(BUFAB(IB1+1),BUFAB(IB1+1),FFTW_PLANS(1))
C
      DO 32 I=0,JC-1
      	II=I+1
      	DO 32 J=-JC+1,JC-1
      	  ITEMP=I**2+J**2
          IF ((ITEMP.GE.MAXR12).AND.(ITEMP.LT.MAXR22)) THEN
      	    JJ=J+1
      	    IF (JJ.LT.1) JJ=JJ+NSAM
            PHASE=SHX*I+SHY*J
            PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
      	      ID=II+JC*(JJ-1)
      	      SAMP=SPEC(ID)*PSHFT
              IB=IB1+2*ID
              AINPO=CMPLX(BUFAB(IB-1),BUFAB(IB))
      	    IBIN=INT(SQRT(REAL(ITEMP))+0.5)+1
            BUFAB(IBIN)=BUFAB(IBIN)+REAL(SAMP*CONJG(AINPO))
            BUFA(IBIN)=BUFA(IBIN)+CABS(SAMP)**2
            BUFB(IBIN)=BUFB(IBIN)+CABS(AINPO)**2
            NS(IBIN)=NS(IBIN)+1
          ENDIF
32    CONTINUE
C
      DO 31 I=1,NSAM/2
        IF(NS(I).NE.0) THEN
          RBIN(2*I-1)=BUFAB(I)/SQRT(BUFA(I)*BUFB(I))
          RBIN(2*I)=1.0
        ENDIF
31    CONTINUE
C
      RETURN
      END
