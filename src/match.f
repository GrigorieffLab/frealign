C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C ------------------------------------------------------------------------	
	SUBROUTINE MATCH(K,NSAM,IRADA,OUTC,OUTD,C3DF,
     +               PHI,THETA,PSI,JC,NSAMH,PHASE,SHX,SHY,
     +		     DSHX,DSHY,XSTD,RI,B3DV,MBUF,ICMP,DATA,
     +		     SPEC,CCPART,PRESA,PRESAOLD,THRESH,
     +		     CFORM,MAXSET,ILIST,IFIRST,NSET,IOPROJ,
     +		     SINCLUT,IPAD,THETATR,IEWALD,
     +               CTFF,PSIZE,FFTW_PLANS,SM,DMASK)
C ------------------------------------------------------------------------	
C         Output of matching projections if requested
C	  CALL MASKENV, FFTW, SHIFT, PEXTRACT, STAMP and IWRITE.
C	  Used in LMAIN.
C ------------------------------------------------------------------------	
C
        USE ISO_C_BINDING
        USE FFTW33
C
        IMPLICIT NONE
        REAL PI,FRAD
	PARAMETER  (PI=3.1415926535897,FRAD=20.0)

	INTEGER NSAM,NSAM1,IRADA,IOPROJ,N1,IDD,ID,IP
        INTEGER ILIST,IFIRST,I,II,J,JJ,ICMP,K,NSET
        INTEGER MAXSET,IPAD,NSAMH,JC

        REAL CCPART,RI,XSTD,DSTD,OSTD
	REAL DATA(*),OUTD(*),PSIZE,DMASK(4)
	REAL B3DV(*),DSHX,DSHY,RAD2,XM,YM,X2,Y2
	REAL PSI,THETA,PHI,MBUF(*),SM(9),DM(9)
	REAL SHX,SHY,PRESA,PRESAOLD,THRESH(*)
	DOUBLE PRECISION SINCLUT(*),PHASE
	COMPLEX OUTC(*),C3DF(*),SPEC(*),CTFF(*)
        CHARACTER*1 CFORM
C MW added parameters for PEXTRACTE (Ewald sphere correction)
C        REAL CS,WL,WGH1,WGH2,DFMID1,DFMID2,ANGAST
        REAL THETATR
        INTEGER IEWALD
        TYPE(C_PTR) FFTW_PLANS(*)
C ------------------------------------------------------------------------	
        write(*,6405) PSI*180./PI,THETA*180./PI,PHI*180./PI
6405    format(' values of PSI,THETA,PHI at FMATCH extraction',3F9.3)

        CALL PEXTRACT(NSAM,IRADA,OUTC,C3DF,PHI,THETA,PSI,
     +	              SINCLUT,IPAD,THETATR,IEWALD,CTFF,SM)

        CALL SHIFT2D(JC,NSAM,NSAMH,0.0,0.0,0.0,0.0,OUTC)

        OUTC(1)=CMPLX(0.0,0.0)
C
        CALL FFTW_BWD(OUTD,OUTD,FFTW_PLANS(2))

        IF (DMASK(4).GT.0.0) THEN
          RAD2=(DMASK(4))**2
          CALL ROTMAT(-PSI,-THETA,-PHI,1.0,DM)
          CALL MATMUL_T(DM,SM,DM)
C          XM=DM(1)*(DMASK(1)-JC)+DM(4)*(DMASK(2)-JC)
C     +      +DM(7)*(DMASK(3)-JC)+1.0
C          YM=DM(2)*(DMASK(1)-JC)+DM(5)*(DMASK(2)-JC)
C     +      +DM(8)*(DMASK(3)-JC)+1.0
          XM=DM(1)*(DMASK(1)-JC)+DM(4)*(DMASK(2)-JC)
     +      +DM(7)*(DMASK(3)-JC)+JC
          YM=DM(2)*(DMASK(1)-JC)+DM(5)*(DMASK(2)-JC)
     +      +DM(8)*(DMASK(3)-JC)+JC
          DO 223 I=1,NSAM
            X2=(I-XM)**2
C            IF (I.GE.NSAMH) THEN
C              X2=(I-XM-NSAM)**2  
C            ENDIF
            DO 223 J=1,NSAM
              Y2=(J-YM)**2
C              IF (J.GE.NSAMH) THEN
C                Y2=(J-YM-NSAM)**2
C              ENDIF
              ID=I+(NSAM+2)*(J-1)
              IF (X2+Y2.GT.RAD2) THEN
                OUTD(ID)=0.0
              ENDIF
223       CONTINUE
        ENDIF
C
        IF (XSTD.LT.0.0) CALL MASKENV(NSAM,RI,OUTD,B3DV,
     +          0.0,0.0,PHI,THETA,PSI,MBUF,0)

C           WRITE OUT PROJECTION
            N1=3.0*RI/ICMP
            IF (N1*ICMP.GT.NSAM) N1=NSAM/ICMP
            CALL SHIFT2D(JC,NSAM,NSAMH,-SHX-DSHX,-SHY-DSHY,
     +                   PI,PI,DATA)
C
C     Low-pass filter....
        DO 88 I=1,JC
          II=I-1
          DO 88 J=1,NSAM
            JJ=J-1
            IF (JJ.GE.JC) JJ=JJ-NSAM
              RAD2=REAL(II**2+JJ**2)/NSAM**2/PSIZE**2*FRAD**2
                ID=I+JC*(J-1)
                SPEC(ID)=SPEC(ID)*EXP(-RAD2)
88      CONTINUE
C
            CALL FFTW_BWD(DATA,DATA,FFTW_PLANS(2))
            CALL WINDOW(NSAM,N1*ICMP,DATA,DSTD)
            CALL WINDOW(NSAM,N1*ICMP,OUTD,OSTD)

            DO 22 I=1,N1
              DO 22 II=1,N1
                DO 22 J=0,ICMP-1
                  DO 22 JJ=0,ICMP-1
                    ID=II+N1*(I-1)
                    IDD=1+ICMP*(II-1)+J+N1*ICMP*(ICMP*(I-1)+JJ)
                    OUTD(ID)=OUTD(IDD)/NSAM/NSAM*2/ICMP**2
     +                       /OSTD*DSTD
22          CONTINUE

            DO 24 I=1,N1
              DO 24 II=1,N1
                DO 24 J=0,ICMP-1
                  DO 24 JJ=0,ICMP-1
                    IDD=1+ICMP*(II-1)+J+N1*ICMP*(ICMP*(I-1)+JJ)
                    ID=II+N1*(I-1+N1)
                    OUTD(ID)=DATA(IDD)/NSAM/NSAM*2/ICMP**2
24          CONTINUE

C            CALL STAMP(REAL(K),6,0,OUTD,N1,2*N1,N1-23,
C     +                 2*N1-5,CFORM)
C            CALL STAMP(PSI/PI*180.0,5,2,OUTD,N1,2*N1,1,
C     +                 N1+1,CFORM)
C            CALL STAMP(THETA/PI*180.0,5,2,OUTD,N1,2*N1,23,
C     +                 N1+1,CFORM)
C            CALL STAMP(PHI/PI*180.0,5,2,OUTD,N1,2*N1,45,
C     +                 N1+1,CFORM)
C            CALL STAMP(PRESA*100.0,7,2,OUTD,N1,2*N1,
C     +                 1,N1,CFORM)
C            CALL STAMP(CCPART,6,1,OUTD,N1,2*N1,
C     +                 1,N1-5,CFORM)
C            CALL STAMP((PRESA-PRESAOLD)*100.0,6,2,OUTD,N1,
C     +                 2*N1,1,N1-5,CFORM)

C            IF (PRESA.GT.ABS(THRESH(NSET))) THEN
C              CALL STAMP(0.0,1,1,OUTD,N1,2*N1,1,2*N1-5,CFORM)
C            ENDIF

C           IP=2*N1*(ABS(ILIST)-1)
            IP=2*N1*(K-IFIRST)
            NSAM1=2*N1

            DO 23 I=1,2*N1
              ID=1+N1*(I-1)
              CALL IWRITE(IOPROJ,OUTD(ID),I+IP)
23          CONTINUE

	RETURN
	END
