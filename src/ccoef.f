C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      REAL FUNCTION CCOEF(NSAM,IRAD,AMAG,A3DF,DATC,
     +                 PHI,THETA,PSI,SHX,SHY,SINCLUT,IPAD,
     +                 IEWALD,THETATR,CTFF,RBUF,FFTW_PLANS,SM)
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
      INTEGER NSAM,L,M,LL,MM,JC,ID,IPAD,IRAD,IB1,ID2
      INTEGER IS,IEWALD
      REAL SHX,SHY,PHASE,SUM1,SUM2,CCPART,SINCLUT(*)
      REAL PHI,THETA,PSI,AMAG,THETATR,RBUF(*),SM(9)
      COMPLEX PSHFT,A3DF(*),CTFF(*)
      COMPLEX DATC(*),Z
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
      IB1=NSAM*(NSAM+2)
      JC=NSAM/2+1
C
      CALL PEXTRACT(NSAM,IRAD,RBUF,A3DF,
     +         PHI,THETA,PSI,SINCLUT,IPAD,
     +         THETATR/AMAG,IEWALD,CTFF,SM)
      DO 222 L=1,JC
        LL=L-1
        DO 222 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          PHASE=-SHX*LL-SHY*MM
          PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
            ID=L+JC*(M-1)
            Z=PSHFT*DATC(ID)
            ID2=2*ID+IB1
            RBUF(ID2-1)=REAL(Z)
            RBUF(ID2)=AIMAG(Z)
222   CONTINUE
      RBUF(IB1+1)=0.0
      RBUF(IB1+2)=0.0
      RBUF(1)=0.0
      RBUF(2)=0.0
      CALL FFTW_BWD(RBUF(IB1+1),RBUF(IB1+1),FFTW_PLANS(2))
      CALL FFTW_BWD(RBUF,RBUF,FFTW_PLANS(2))
C
      SUM1=0.0
      SUM2=0.0
      CCPART=0.0
      DO 223 L=1,NSAM
        LL=(NSAM+2)*L
        DO 223 M=1,NSAM
          ID=M+LL
          IS=IB1+ID
          CCPART=CCPART+RBUF(ID)*RBUF(IS)
          SUM1=SUM1+RBUF(ID)**2
          SUM2=SUM2+RBUF(IS)**2
223   CONTINUE
C
      CCOEF=CCPART/SQRT(SUM1)/SQRT(SUM2)
C
      RETURN
      END
