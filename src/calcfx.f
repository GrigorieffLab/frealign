C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C******************************************************************************
      SUBROUTINE CALCFX(NX,XPARMAP,RF,XPAR,NPAR,NSAM,MAXR1,MAXR2,
     .		MASK,A3DF,IRAD,
     .		PBUF,SHX,SHY,CS,WL,WGH1,WGH2,THETATR,CTFF,
     . 		AMAG,RIH,HALFW,RI2,RI3,RI4,PHI,THETA,PSI,RMAX1,RMAX2,
     .		XSTD,MBUF,ILST,DATC,IBUF,B3DV,DATD,
     .	        SINCLUT,IVFLAG,RBFACT,OUTD,OUTC,QBUC,IPAD,
     .          RBUF,FSCT,DFMID1,DFMID2,ANGAST,IEWALD,TX,TY,XM,YM,
     .          SX,SY,SIG2,DFSTD,ASYM,THETAM,STHETA,PSIM,SPSI,STIF,
     .          PWEIGHTS,PSSNR,FFTW_PLANS,SM)
C******************************************************************************
C     CALCULATES NEW VALUE FOR F TO INPUT TO SUBROUTIN VA04A
C     There are three modes - in all the score is returned as value of RF
C	IF IVFLAG=0, the score is calculated with 5 new particle parameters
C	IF IVFLAG=1, the score is recalculated with old particle parameters
C		after application of new, but unweighted CTF correction (+/-1.0)
C		for use in defocus refinement
C	IF IVFLAG=2, the score is recalculated with old particle parameters
C		but using the new MAG and unweighted CTF correction (i.e. +/-1.0)
C		for use in magnification refinement

C      Calls CTFAPPLY_PHASE_ONLY and uses functions CC3, CC3M.
C      Used in VA04.
C      MW passed thru THETATR and IEWALD to CC3,CC3M for Ewald sphere corr.     
C******************************************************************************
C
      USE ISO_C_BINDING
C
      IMPLICIT NONE

      INTEGER NSAM,MAXR1,MAXR2,NX,IRAD,IVFLAG,I,K
      INTEGER ILST(*),IBUF,ID,IS,IDD,NCALCFX_CTF
      INTEGER MASK(5),NPAR,IPAR,JPAR,IPAD,IEWALD

      REAL XPAR(*),RF,CC3M,CC3,RFT,RBUF(*),TX,TY
      REAL XPARMAP(*),OUTD(*),SHX(*),SHY(*),SIG2,SPSI
      REAL PI,DMAG,PBUF(*),MBUF(*),FSCT(*),SX,SY,PSIM
      REAL CS,WL,WGH1,WGH2,DFMID1(*),DFMID2(*),DFSTD
      REAL RIH,HALFW,RI2,RI3,RI4,AMAG(*),XM,YM,STHETA
      REAL PHI(*),THETA(*),PSI(*),RMAX1,RMAX2,THETATR
      REAL ANGAST(*),DATD(*),B3DV(*),XSTD,FSH,FD,SM(9)
      REAL PTEMP,RBFACT,SINCLUT(*),DF1,DF2,AST,THETAM
      REAL FANGLE,CC3_C,CC3M_C,STIF,PWEIGHTS(*),PSSNR(*)

      PARAMETER (PI=3.1415926535897)

      COMPLEX OUTC(*),DUMC,CTFF(*)
      COMPLEX A3DF(*),QBUC(*),DATC(*)

      CHARACTER ASYM*3

      TYPE(C_PTR) FFTW_PLANS(*)
C***********************************************************************
      NCALCFX_CTF=0
      IF (IABS(IVFLAG).EQ.1) THEN
C     for defocus and astigmatism, parameters map directly
C	write(*,*)'Entering CALCFX, CTF section'
      	RFT=0.0
      	DO 10 I=1,IBUF
      		K=ILST(I)
      		IDD=1+NSAM/2*NSAM*(I-1)
      		ID=1+NSAM*NSAM*(I-1)
      		IS=1+NSAM*(I-1)
        	MAXR1=INT(NSAM*RMAX1*ABS(AMAG(K)))
        	MAXR2=INT(NSAM*RMAX2*ABS(AMAG(K)))
      		IF(NX.EQ.1) THEN
      			XPARMAP(2)=XPARMAP(1)
      		ENDIF
                DF1=DFMID1(K)+XPARMAP(1)
                DF2=DFMID2(K)+XPARMAP(2)
                AST=ANGAST(K)+XPARMAP(3)
        	CALL CTFAPPLY(NSAM,PBUF(ID),SHX(K),SHY(K),
     +		 CS,WL,WGH1,WGH2,DF1,DF2,AST,
     +		 THETATR,CTFF,OUTD,OUTC,AMAG(K),RIH,HALFW,
     +		 RI2,RI3,RI4,DATD,DATC,B3DV,PHI(K),THETA(K),
     +           PSI(K),MBUF(ID),XSTD,IVFLAG,TX,TY,ASYM,PSSNR,RBUF,
     +           PWEIGHTS,0,FFTW_PLANS)
                IF (ASYM(1:1).EQ.'H') THEN
C
        	RFT=RFT-CC3_C(NSAM,IRAD,AMAG(K),OUTC,A3DF,
     +	         MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),0.0,0.0,IVFLAG,
     +		 QBUC(IDD),RBFACT,SINCLUT,IPAD,RBUF,FSCT,
     +           IEWALD,THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,
     +           FFTW_PLANS,SM)
     +           -FD(SIG2,XPARMAP(1),XPARMAP(2),DFSTD)
C
                ELSE
C
        	RFT=RFT-CC3(NSAM,IRAD,AMAG(K),OUTC,A3DF,
     +	         MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),0.0,0.0,IVFLAG,
     +		 QBUC(IDD),RBFACT,SINCLUT,IPAD,RBUF,FSCT,
     +           IEWALD,THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,
     +           FFTW_PLANS,SM)
     +           -FD(SIG2,XPARMAP(1),XPARMAP(2),DFSTD)
C
                ENDIF
10    	CONTINUE
      	RFT=RFT/REAL(IBUF)
C      	IVFLAG=-1
      	NCALCFX_CTF=NCALCFX_CTF+1
C       print out occasional values for interest
C      	IF(NCALCFX_CTF.LE.5.OR.NCALCFX_CTF.EQ.10.OR.NCALCFX_CTF/50*50
C     .	.EQ.NCALCFX_CTF) print *,XPARMAP(1),XPARMAP(2),
C     .	XPARMAP(3)/pi*180.0,rft/pi*180.0

      ELSEIF (IVFLAG.EQ.2) THEN
C  for magnification, parameters map directly
C	write(*,*)'Entering CALCFX, MAG section',IBUF
      	RFT=0.0
      	DO 20 I=1,IBUF
      		K=ILST(I)
      		IDD=1+NSAM/2*NSAM*(I-1)
      		ID=1+NSAM*NSAM*(I-1)
      		IS=1+NSAM*(I-1)
        	MAXR1=INT(NSAM*RMAX1*ABS(XPARMAP(1)))
        	MAXR2=INT(NSAM*RMAX2*ABS(XPARMAP(1)))
        	CALL CTFAPPLY(NSAM,PBUF(ID),SHX(K),SHY(K),
     +		 CS,WL,WGH1,WGH2,DFMID1(K),DFMID2(K),ANGAST(K),
     +		 THETATR,CTFF,OUTD,OUTC,XPARMAP(1),RIH,
     +		 HALFW,RI2,RI3,RI4,DATD,DATC,B3DV,PHI(K),
     +           THETA(K),PSI(K),MBUF(ID),XSTD,IVFLAG,TX,TY,ASYM,
     +           PSSNR,RBUF,PWEIGHTS,0,FFTW_PLANS)
                IF (ASYM(1:1).EQ.'H') THEN
C
        	PTEMP=-CC3_C(NSAM,IRAD,XPARMAP(1),OUTC,A3DF,
     +	         MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),0.0,0.0,IVFLAG,
     +		 QBUC(IDD),RBFACT,SINCLUT,IPAD,RBUF,FSCT,
     +           IEWALD,THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,
     +           FFTW_PLANS,SM)
C
                ELSE
C
        	PTEMP=-CC3(NSAM,IRAD,XPARMAP(1),OUTC,A3DF,
     +	         MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),0.0,0.0,IVFLAG,
     +		 QBUC(IDD),RBFACT,SINCLUT,IPAD,RBUF,FSCT,
     +           IEWALD,THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,
     +           FFTW_PLANS,SM)
C
                ENDIF
      		RFT=RFT+PTEMP
20    	CONTINUE
      	RFT=RFT/REAL(IBUF)
C	write(*,21) rft*100.0,xparmap(1)
21	format('     exiting CALCFX - rf, relative magnification',7F9.4)

      ELSE
C	write(*,*)'Entering CALCFX, general section'

C  pack variables being refined (XPARMAP) into XPAR using MASK
      	IPAR=0
      	DO 500 JPAR=1,5
      		IF(MASK(JPAR).EQ.1) THEN
      			IPAR=IPAR+1
      			XPAR(JPAR)=XPARMAP(IPAR)
      		ENDIF
500	CONTINUE
      		IF(NPAR.GT.IPAR) XPAR(6)=XPARMAP(NPAR)

                IF (ASYM(1:1).EQ.'H') THEN
C
      	RFT=-CC3M_C(NSAM,IRAD,XPAR(6),OUTC,A3DF,
     +	       MAXR1,MAXR2,XPAR(1),XPAR(2),XPAR(3),XPAR(4),XPAR(5),0,
     +	       DUMC,RBFACT,SINCLUT,IPAD,RBUF,FSCT,IEWALD,
     +         THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,FFTW_PLANS,SM)
     +        -FSH(SIG2,XPAR(4),XPAR(5),XM,YM,SX,SY)
     +        -FANGLE(SIG2,XPAR(2),THETAM,STHETA,STIF)
     +        -FANGLE(SIG2,XPAR(3),PSIM,SPSI,STIF)
C
        ELSE
C
      	RFT=-CC3M(NSAM,IRAD,XPAR(6),OUTC,A3DF,
     +	       MAXR1,MAXR2,XPAR(1),XPAR(2),XPAR(3),XPAR(4),XPAR(5),0,
     +	       DUMC,RBFACT,SINCLUT,IPAD,RBUF,FSCT,IEWALD,
     +         THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,FFTW_PLANS,SM)
     +        -FSH(SIG2,XPAR(4),XPAR(5),XM,YM,SX,SY)
C
        ENDIF
      	IF (ABS(XPAR(1)).GT.PI) RFT=RFT+(PI-ABS(XPAR(1)))**2
      	IF (ABS(XPAR(2)).GT.PI) RFT=RFT+(PI-ABS(XPAR(2)))**2
      	IF (ABS(XPAR(3)).GT.PI) RFT=RFT+(PI-ABS(XPAR(3)))**2
      	IF (ABS(XPAR(4)).GT.1.0) RFT=RFT+(1.0-ABS(XPAR(4)))**2
      	IF (ABS(XPAR(5)).GT.1.0) RFT=RFT+(1.0-ABS(XPAR(5)))**2
C		DMAG=ABS(1.0-XPAR(6))
C		IF (DMAG.GT.0.02) RFT=RFT+(10.0*(DMAG-0.02))**2
      ENDIF

      RF=RFT
C	write(*,25) rf*100.0,xpar,IVFLAG
25	format(' exiting CALCFX - rf,xpar,IVFLAG',7F10.3,I3)

      RETURN
      END
