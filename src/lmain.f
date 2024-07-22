C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -------------------------------------------------------------------
      SUBROUTINE LMAIN(KK,FDEF,FMAG,INPROJ,CFORM,DATA,
     +  RREC,NSAM,SPEC,SHX,SHY,CS,RAWBINS,WL,WGH1,WGH2,
     +  DFMID1,DFMID2,ANGAST,OCC,OUTC,OUTD,AMAGP,RIH,HALFW,HALFWC,
     +  RI2,RI3,RI4,DATD,DATC,B3DV,PHI,THETA,PSI,XSTD,IFLAG,
     +	IRADA,C3DF,RBINS,RBINN,NSET,PRESA,SANG,
     +  NSANG,CCD,CCC,IQUADMAX,NDOC1,NDOC2,PMASK,RBFAC,RMAX1,RMAX2,
     +  IRAN,PSIZE,DSTEP,FILM,NN1,NNPART,MAXSET,FHIST,ILIST,NPB,
     +	ITMAX,DAMAX,ASYMSTORE,TARGET,ICCOUNT,ISWITCH,IPAD,
     +	PALL,IALL,IRAD,RI,A3DF,S3DF,STD,D3DF,ILST,
     +  V3DF,VTD,PBC,BOFF,ASUM,VSUM,PSUM,KSUM,THRESH,
     +	NSAMH,JC,FFLAG,CCAVE,ICMP,IFIRST,IOPROJ,IEWALD,
     +	NSCF,SFILE,NANG,WGH,PSIS,THETAS,PHIS,DSHXS,DSHYS,PRESS,
     +  IEXCL,IBUF,ICOUNT,FMATCH,DFMID1S,DFMID2S,ANGASTS,OCCS,
     +	ALGPS,DPRESS,ABSMAGS,NSTART,FINPAT1,FINPAT2,ASYM,VX,ILAST,
     +	NNSET,IANG,CDATE,CTIME,CZONE,DTVAL,FASTIG,SINCLUT,MBUF,PBUF,
     +  CTFF,RBIN,CCBUF,QBUC,
     +  ISTAT,TNSET,FALL,IC,VSN,VN,VVSN,VVN,CIA,BF,NS,NA,FSCT,TX,TY,
     +  TESTPAR,IPMAX,NSYM,ISYMAX,JSYM,SYMOP,NNSTAT,NKSUM,XM,YM,
     +  SX,SY,IFSC,IMP,NVOL,NVOL1,NVOL2,OALL,DFSTD,ALPHA,RISE,NU,
     +  HSTART,THETAMM,STHETA,PSIMM,SPSI,STIF,NPIX,NS1,ALGP,DALGPS,
     +  DMASK,NKSUM1,TTD1,TTD2,SIG,INTERP,PSSNR,PWEIGHTS,ISW,RECBF,
     +  IRECB,RECB,CTFBF,NN2,NN3,WALL,IWALL,RCLAS,DANGIN,
     +  FFTW_PLANS,NONSYM,ISM,SM,IREDUN)
C ---------------------------------------------------------------------
#ifdef _NAG
        USE F90_UNIX
#endif
C
        USE ISO_C_BINDING
        USE FFTW33
C
        IMPLICIT NONE

        REAL PI
        PARAMETER  (PI=3.1415926535897)

        INTEGER IANG,MAXSET,NN1,NSET,NNPART,IFSC

        CHARACTER*1 CFORM
        CHARACTER*200 SFILE
        CHARACTER*200 FINPAT,FINPAT1(*),FINPAT2(*)
        CHARACTER*3 ASYM
        CHARACTER*15 VX
        CHARACTER CDATE*8,CTIME*10,CZONE*5

        LOGICAL FDEF,FMAG,FHIST,ASYMSTORE,FMATCH,FASTIG
        LOGICAL FALL,CIA(*),CI,NONSYM

        INTEGER I,J,IS,ID,L,LL,M,MM,IBUFB,IEWALD,II,JJ
        INTEGER ITT,IP,MAXRES,IFIRSTS,IPMAX,NNSTAT,NN3
        INTEGER LASTF,IFIRSTF,ISTAT,TNSET(*),NS1(*),NVOL
        INTEGER PMASK(5),NPAR,NPB,IFIRST,IMP,NVOL1,NVOL2
        INTEGER FILM(*),ILST(*),ILIST(*),JC,ISWITCH
        INTEGER KSUM(NN3*NNSTAT+1,NVOL),IWALL,ISM,KK
        INTEGER ICOUNT,NANG,NSCF,HSTART,IRECB,NN2
        INTEGER IVFLAG,NSYM,ISYMAX,JSYM(*),INTERP
        INTEGER DTVAL(8),NNSET(*),IALL,IRAD,NSAMH,NSAM
        INTEGER NSANG,NSAM1,NSTART,IC(*),NS(*),NPIX(*)
        INTEGER NDOC1,NDOC2,N1,N2,N3,INPROJ,IOPROJ,ICMP
        INTEGER IQUADMAX,IRAN,ICCOUNT,ITMAX,ILAST,IRADA
        INTEGER IFLAG,IEXCL,IBUF,K,FFLAG,ISW(*)
        INTEGER MAXR1,MAXR2,IPAD,NA,NU,IREDUN
        INTEGER*8 NKSUM,NKSUM1(*)

        REAL PSIZE,CS(*),WL(*),PBC,PMIN,PMAX,TESTPAR(6,*),LOGPC
        REAL RREC(*),AMAGP(*),DSHXS,DSHYS,MX,MY,SIGX,SIGY,PSX
        REAL PHIS,PSIS,THETAS,DFMID1S,DFMID2S,ANGASTS,ABSMAGS,PRESS
        REAL SINCLUT(*),TX(*),TY(*),PSIZE1,OUTD(*),SIG2,SIG2N,FSH
        REAL PSI(*),THETA(*),PHI(*),PRESA(*),DPSI,DTHETA,DPHI
        REAL DATA(*),DSTEP(*),SYMOP(3,3,*),BF(*),FSCT(*)
        REAL RBIN(*),RBINS(*),RBINN(*),RAWBINS(*),RMAX1(*),RMAX2(*)
        REAL DFMID1(*),DFMID2(*),ANGAST(*),TARGET(*),RBFAC(*)
        REAL B3DV(*),MBUF(*),SHX(*),SHY(*),PBUF(*),DFSTD(*),OCC(*)
        REAL S3DF(NN3,NVOL1),RCLAS(*),SM(9)
        REAL V3DF(NN3,NVOL2),CCOEF,DANGIN
        REAL ASUM(NN3*NNSTAT+1,NVOL),PSIMM(*),WALL
        REAL VSUM(NN3*NNSTAT+1,NVOL),SPSI(*),SIGP,AVEP
        REAL PSUM(NN3*NNSTAT+1,NVOL),THRESH(*),CSN,WLN,OCCS
        REAL PHIOLD,THETAOLD,PSIOLD,DANGAST,PRESAOLD,AMAGOLD
        REAL PHIN,THEN,PSIN,PRESN,AMAGN,PHIM,THETAM,PSIM,AMAGM
        REAL DSHX,DSHY,DSHXN,DSHYN,ALPHA,RISE,ALGPOLD,DALGPS
        REAL BESTPHI,BESTTHETA,BESTPSI,BESTAMAG,BESTP,BESTANGAST
        REAL DDFMID1,DDFMID2,BESTDFMID1,BESTDFMID2,SIGOLD,ALGPS
        REAL DSHXM,DSHYM,BESTDSHX,BESTDSHY,STHETA(*),DMASK(*)
        REAL DATD(*),XM(*),YM(*),SX(*),SY(*),CCD(*),CC3M,OALL
        REAL HALFW,HALFWC,XSTD,SANG(3,*),THETAMM(*),CC3M_C,DPRESS
        REAL CCPART,PALL,BOFF,PHASE,RIH,CCAVE,RANDOM,THETATR
        REAL WGH,WGH1,WGH2,RI,RI2,RI3,RI4,DAMAX,STIF,ALGP(*)
        REAL PRESMIN,RMAX1N,RMAX2N,RBFACT,DSMAX,FANGLE,SIG(*)
        REAL PSSNR(*),PWEIGHTS(*),RECB(10,NVOL),REMAP_THETA
        PARAMETER (DSMAX=4.0)

        DOUBLE PRECISION STD,VTD,VSN(*),VN(*),VVSN(*),VVN(*)
        DOUBLE PRECISION TTD1(*),TTD2(*)

        COMPLEX CTFF(*),DATC(*),OUTC(*)
        COMPLEX SPEC(*),CCBUF(*)
        COMPLEX A3DF(NN3,NVOL1),C3DF(*),D3DF(NN3,NVOL2)
        COMPLEX QBUC(*),CCC(*),PSHFT,DUMC
        COMPLEX RECBF(NN1/2*NN1,NVOL),CTFBF(NN1*NN1,NVOL)

        TYPE(C_PTR) FFTW_PLANS(*)

        SAVE AMAGOLD
c -----------------------------------------------------------------------
C
          CALL FLUSH(6)
C
          IF (NONSYM) THEN
            K=ISM+(KK-1)*IREDUN
          ELSE
            K=KK
            DO 837 J=K,2,-1
              IF (FILM(J).NE.FILM(J-1)) THEN
                IFIRSTF=J
                IF (J.EQ.K) IBUF=0
                GOTO 836
              ENDIF
837         CONTINUE
          ENDIF
          IFIRSTF=1
836       CONTINUE
          MAXR1=INT(NSAM*RMAX1(NSET)*ABS(AMAGP(K)))
          MAXR2=INT(NSAM*RMAX2(NSET)*ABS(AMAGP(K)))
C
          IFIRSTS=ABS(ILIST(1))
          CCPART=0.0
          FINPAT=FINPAT1(NSET)
          RBFACT=RBFAC(NSET)/(4.0*(NSAM*PSIZE*AMAGP(K))**2)
          CI=CIA(NSET)
          IF ((KK.EQ.NSTART).AND.(ISM.EQ.1)) THEN
            CALL IOPEN(FINPAT,INPROJ,CFORM,N1,N2,N3,'OLD',
     .                  ASYM,PSIZE1,VX)
C            IF (ABS(PSIZE-PSIZE1)/PSIZE.GT.0.01)
C     +        WRITE(*,*)
C     +        ' ***WARNING*** Input and file pixel sizes differ!'
            IF (N3.GT.TNSET(NSET)) WRITE(*,2033)N3,TNSET(NSET)
2033        FORMAT(' ***WARNING*** no. of particles in stack',
     +             ' exceeds no. of lines in parameter file:',
     +             ' NSTACK, NPAR = ',2I8)
            IF (N3.LT.TNSET(NSET)) THEN
              WRITE(*,2034)N3,TNSET(NSET)
2034          FORMAT(' ERROR: no. of particles in stack',
     +             ' smaller than no. of lines in parameter file:',
     +             ' NSTACK, NPAR = ',2I8)
              CALL FLUSH(6)
              STOP 'Fatal error'
            ENDIF
            CALL FLUSH(6)
C
            IF (FMATCH) THEN
              N1=INT(3.0*RI/ICMP)
              IF (N1*ICMP.GT.NSAM) N1=NSAM/ICMP
              N2=N1
              N3=NNSET(NSET)
              NSAM1=2*N1
              CALL IOPEN(FINPAT2(NSET),IOPROJ,CFORM,N1,
     +        2*N1,N3,'NEW',ASYM,PSIZE,VX)
            ENDIF
            IF (IFLAG.NE.0) THEN
              WRITE(NDOC1+NSET,6707)
              CALL FLUSH(NDOC1+NSET)
              WRITE(NDOC2+NSET,6709)
              CALL FLUSH(NDOC2+NSET)
6707          FORMAT('C',11X,'PSI   THETA     PHI       SHX       SHY ',
     +      '    MAG  FILM      DF1      DF2  ANGAST     OCC',
     +      '      LogP      SIGMA   SCORE  CHANGE')
6709          FORMAT('C',11X,'PSI   THETA     PHI       SHX       SHY ',
     +      '    MAG  FILM      DF1      DF2  ANGAST     OCC',  
     +      '      LogP      SIGMA   SCORE')
            ENDIF
          ENDIF
          IF (KK.GT.NSTART) THEN
            IF (REAL(ILIST(K))*REAL(ILIST(K-1)).LT.0.0) THEN
              IF (NONSYM) THEN
                WRITE(*,6500)NSET,INT(IEXCL/IREDUN)
              ELSE
                WRITE(*,6500)NSET,IEXCL
              ENDIF
6500          FORMAT(/' PARTICLES BELOW THRESHOLD (NOT INCLUDED)',
     +               ' IN SET ',I3,': ',I10)
              IEXCL=0
              NSET=NSET+1
              IFIRST=KK
              IFIRSTS=ABS(ILIST(K))
              CALL ICLOSE(INPROJ)
              FINPAT=FINPAT1(NSET)
              CALL IOPEN(FINPAT,INPROJ,CFORM,N1,N2,N3,'OLD',ASYM,
     .                  PSIZE1,VX)
C              IF (ABS(PSIZE-PSIZE1)/PSIZE.GT.0.01)
C     +          WRITE(*,*)
C     +          ' ***WARNING*** Input and file pixel sizes differ!'
              CALL FLUSH(6)
              IF (N3.GT.TNSET(NSET)) WRITE(*,2033)N3,TNSET(NSET)
              IF (N3.LT.TNSET(NSET)) THEN
                WRITE(*,2034)N3,TNSET(NSET)
                CALL FLUSH(6)
                STOP 'Fatal error'
              ENDIF
              IF (FMATCH) THEN
                CALL ICLOSE(IOPROJ)
                N1=INT(3.0*RI/ICMP)
                IF (N1*ICMP.GT.NSAM) N1=NSAM/ICMP
                N2=N1
                N3=NNSET(NSET)
                NSAM1=2*N1
                CALL IOPEN(FINPAT2(NSET),IOPROJ,CFORM,N1,
     +          2*N1,N3,'NEW',ASYM,PSIZE,VX)
              ENDIF
              WRITE(NDOC1+NSET,6707)
              CALL FLUSH(NDOC1+NSET)
              WRITE(NDOC2+NSET,6709)
              CALL FLUSH(NDOC2+NSET)
            ENDIF
          ENDIF
C
          DSHX=0.0
          DSHY=0.0
          PHIOLD=PHI(K)
          THETAOLD=THETA(K)
          PSIOLD=PSI(K)
          PRESAOLD=PRESA(K)
          SIGOLD=SIG(K)
          ALGPOLD=ALGP(K)
          IF (.NOT.FDEF) THEN
            DDFMID1=0.0
            DDFMID2=0.0
            DANGAST=0.0
          ENDIF
C
          IF(KK.EQ.NSTART.OR.KK.EQ.NSTART+1.OR.KK.EQ.NSTART+10.OR.
     .      KK.EQ.NSTART+100.OR.KK.EQ.NSTART+1000.OR.
     .      KK.EQ.NSTART+10000.OR.KK.EQ.ILAST) THEN
            CALL DATE_AND_TIME(CDATE,CTIME,CZONE,DTVAL)
            WRITE(*,6708) KK,CTIME(1:2),CTIME(3:4),CTIME(5:6)
6708        FORMAT(' Time before particle',I10,' was  ',A2,':',
     .            A2,':',A2)
            CALL FLUSH(6)
          ENDIF
C
C         FIND FIRST PARTICLE FROM THIS FILM
C
          IVFLAG=0
          IF ((FDEF).AND.(IBUF.EQ.0)) THEN
C                IFIRSTF=KK
                IVFLAG=1
C                DFMID1OLD=DFMID1(K)
C                DFMID2OLD=DFMID2(K)
C                ANGASTOLD=ANGAST(K)
          ELSEIF ((FDEF).AND.(CI)) THEN
                IVFLAG=1
C                DFMID1OLD=DFMID1(K)
C                DFMID2OLD=DFMID2(K)
C                ANGASTOLD=ANGAST(K)
C          ELSEIF (.NOT.(FDEF)) THEN
C                DFMID1OLD=DFMID1(K)
C                DFMID2OLD=DFMID2(K)
C                ANGASTOLD=ANGAST(K)
          ENDIF
          IF ((FMAG).AND.(IBUF.EQ.0)) THEN
C                IFIRSTF=KK
                IVFLAG=2
                AMAGOLD=AMAGP(K)
          ELSEIF (.NOT.(FMAG)) THEN
                AMAGOLD=AMAGP(K)
          ENDIF
C
          THETATR=WL(NSET)/(PSIZE*NSAM)
C
          IF (((FDEF.AND.(.NOT.CI)).OR.FMAG).AND.(IBUF.EQ.0)) THEN 
              IF (NONSYM) THEN
                WRITE(*,*)' ERROR: defocus/mag refinement not',
     +                    ' possible with asymmetric refinement.'
                CALL FLUSH(6)
                STOP 'Fatal error'
              ENDIF
C read ahead to get all particles on this film for magnification or defocus refinement
              DO 840 J=1,NPB
C
                IF (IFIRSTF+J-1.GT.NNPART) GOTO 841
                IF (IFIRSTF+J-1.GT.NANG) GOTO 841
                IF (FILM(IFIRSTF+J-1).NE.FILM(K)) GOTO 841
                LASTF=IFIRSTF+J-1
                IF (PRESA(IFIRSTF+J-1).LE.-1.0) GOTO 840
                IBUF=IBUF+1
                ILST(IBUF)=IFIRSTF+J-1
C                IP=NSAM*(ABS(ILIST(IFIRSTF+J-1))-IFIRSTS)
                IP=NSAM*(ABS(ILIST(IFIRSTF+J-1))-1)
                DO 20 I=1,NSAM
                  ID=1+(NSAM+2)*(I-1)+NSAM*(NSAM+2)*(IBUF-1)
C                  WRITE(*,*)' IREAD1 IP = ',IP
                  CALL IREAD(INPROJ,PBUF(ID),I+IP)
20              CONTINUE
                PMIN=1.0E30
                PMAX=-1.0E30
                AVEP=0.0
                SIGP=0.0
                DO 21 I=1,NSAM
                  ID=(NSAM+2)*(I-1)+NSAM*(NSAM+2)*(IBUF-1)
                  DO 21 II=1,NSAM
                    IF(PBUF(II+ID).GT.PMAX) PMAX=PBUF(II+ID)
                    IF(PBUF(II+ID).LT.PMIN) PMIN=PBUF(II+ID)
                    AVEP=AVEP+PBUF(II+ID)
                    SIGP=SIGP+PBUF(II+ID)**2
21              CONTINUE
                IF(PMAX.EQ.PMIN) THEN
                  WRITE(*,*)' ERROR: blank image, N = ',
     +                      IFIRSTF+J-1
                  CALL FLUSH(6)
                  STOP 'Fatal error'
                ENDIF
                IF((PMAX.LE.0.0).OR.(PMIN.GE.0.0)) THEN
                  WRITE(*,*)' ERROR: particle not normalized,',
     +                      ' N = ',IFIRSTF+J-1
                  CALL FLUSH(6)
                  STOP 'Fatal error'
                ENDIF
                SIGP=SQRT(SIGP/NSAM/NSAM-(AVEP/NSAM/NSAM)**2)
C                IF ((SIG(K)/SIGP.GT.3.0).OR.(SIG(K)/SIGP.LT.0.33))
C     +            THEN
C                  IF (SIG(K).NE.0.0)
C     +      WRITE(*,*)' Resetting unrealistic sigma for particle ',K
C                  SIG(K)=SIGP
C                ENDIF
C
                ID=1+NSAM*(NSAM+2)*(IBUF-1)
                IS=1+NSAM*(IBUF-1)
                CALL FFTW_FWD(PBUF(ID),PBUF(ID),FFTW_PLANS(1))
C               Need to set PBUQ=0 since nonzero terms will be
C               incompatible with shift operation this is the
C               result of using a specialized real fft)
C                DO 839 I=1,NSAM
C                  PBUQ(I)=CMPLX(0.0,0.0)
C839             CONTINUE
C
840             CONTINUE
841         CONTINUE
          ENDIF
C
          IF ((FMAG).AND.(IVFLAG.EQ.2)) THEN
             IF (NONSYM) THEN
              WRITE(*,*)' ERROR: magnification refinement not',
     +                  ' possible with asymmetric refinement.'
              CALL FLUSH(6)
              STOP 'Fatal error'
            ENDIF
C done only at first particle on each new film number
            RMAX1N=RMAX1(NSET)
            RMAX2N=RMAX2(NSET)
            CSN=CS(NSET)
            WLN=WL(NSET)
            WRITE(*,6012) IFIRSTF,LASTF,FILM(K),K
6012        FORMAT(' Refine magnification for particles',I7,' to',I7,
     .             '  on film',I6,'     (K=',I7,')')
            CALL FLUSH(6)
            CALL MAGREFINE(AMAGP,DFMID1,DFMID2,ANGAST,IFIRSTF,LASTF,
     +			NPAR,NSAM,MAXR1,MAXR2,PMASK,C3DF,IRADA,PBUF,
     +			SHX,SHY,IPAD,CSN,WLN,WGH1,WGH2,THETATR,
     +                  CTFF,AMAGP(IFIRSTF),RIH,HALFWC,RI2,RI3,RI4,
     +			PHI,THETA,PSI,RMAX1N,RMAX2N,XSTD,MBUF,
     +			ILST,DATC,IBUF,B3DV,DATD,SINCLUT,IVFLAG,
     +			RBFACT,OUTD,OUTC,QBUC,BF,FSCT,FILM,ILAST,
     +                  IEWALD,TX(NSET),TY(NSET),ASYM,PWEIGHTS,PSSNR,
     +                  FFTW_PLANS,SM)

            WRITE(*,6013) (1.0/AMAGP(IFIRSTF))*DSTEP(NSET)*10000./PSIZE
6013        FORMAT('                      absolute magnification =',
     .      F9.1)
            CALL FLUSH(6)
            IF (FDEF) THEN
                IVFLAG=1
            ELSE
                IVFLAG=0
            ENDIF
          ENDIF
C
          IF ((FDEF.AND.(.NOT.CI)).AND.(IVFLAG.EQ.1)) THEN
             IF (NONSYM) THEN
              WRITE(*,*)' ERROR: defocus refinement not',
     +                  ' possible with asymmetric refinement.'
              CALL FLUSH(6)
              STOP 'Fatal error'
            ENDIF
C done only at first particle on each new film number
            RMAX1N=RMAX1(NSET)
            RMAX2N=RMAX2(NSET)
            CSN=CS(NSET)
            WLN=WL(NSET)
            WRITE(*,6011) IFIRSTF,LASTF,FILM(K),K
6011        FORMAT(' Refine defocus for particles',I7,' to',I7,
     .             '  on film',I6,'     ( K=',I7,')')
            CALL FLUSH(6)
            CALL CTFREFINE(DFMID1,DFMID2,ANGAST,IFIRSTF,LASTF,FASTIG,
     +		NPAR,NSAM,MAXR1,MAXR2,PMASK,DDFMID1,DDFMID2,DANGAST,
     +	        C3DF,IRADA,PBUF,SHX,SHY,IPAD,
     +		CSN,WLN,WGH1,WGH2,THETATR,CTFF,AMAGP,RIH,
     +		HALFWC,RI2,RI3,RI4,PHI,THETA,PSI,
     +          RMAX1N,RMAX2N,XSTD,MBUF,ILST,DATC,IBUF,B3DV,
     +	        DATD,SINCLUT,IVFLAG,RBFACT,OUTD,OUTC,QBUC,CI,
     +	 	BF,FSCT,FILM,ILAST,IEWALD,TX(NSET),TY(NSET),0.0,
     +          DFSTD(NSET),ASYM,PWEIGHTS,PSSNR,FFTW_PLANS,SM)
C
            IVFLAG=0
          ENDIF
C
C criterion for inclusion when THRESH < 0 (not use anymore)
C        IF (((THRESH(NSET).LT.0.0).AND.(PRESA(K).GT.-1.0))
C     +     .OR.(THRESH(NSET).GE.0.0)) THEN

C          IF (((.NOT.FDEF).OR.CI).AND.(.NOT.FMAG)) THEN

C then need to read in particle images because they were not
C already read in for defocus or magnification refinement
              IF (ISM.EQ.1) THEN
C                IP=NSAM*(ABS(ILIST(K))-IFIRSTS)
                IP=NSAM*(ABS(ILIST(K))-1)
                DO 19 I=1,NSAM
                   ID=1+(NSAM+2)*(I-1)
C                   WRITE(*,*)' IREAD2 IP = ',IP
                   CALL IREAD(INPROJ,DATA(ID),I+IP)
19              CONTINUE
                PMIN=1.0E30
                PMAX=-1.0E30
                AVEP=0.0
                SIGP=0.0
                DO 22 I=1,NSAM
                  ID=(NSAM+2)*(I-1)
                  DO 22 II=1,NSAM
                    IF(DATA(II+ID).GT.PMAX) PMAX=DATA(II+ID)
                    IF(DATA(II+ID).LT.PMIN) PMIN=DATA(II+ID)
                    AVEP=AVEP+DATA(II+ID)
                    SIGP=SIGP+DATA(II+ID)**2
22              CONTINUE
                IF(PMAX.EQ.PMIN) THEN
                  WRITE(*,*)' ERROR: blank image, N = ',KK
                  CALL FLUSH(6)
                  STOP 'Fatal error'
                ENDIF
                IF((PMAX.LE.0.0).OR.(PMIN.GE.0.0)) THEN
                  WRITE(*,*)' ERROR: particle not normalized,',
     +                      ' N = ',KK
                  CALL FLUSH(6)
                  STOP 'Fatal error'
                ENDIF
                SIGP=SQRT(SIGP/NSAM/NSAM-(AVEP/NSAM/NSAM)**2)
C                IF ((SIG(K)/SIGP.GT.3.0).OR.(SIG(K)/SIGP.LT.0.33))
C     +            THEN
C                  IF (SIG(K).NE.0.0)
C     +      WRITE(*,*)' Resetting unrealistic sigma for particle ',K
C                  SIG(K)=SIGP
C                ENDIF
                CALL FFTW_FWD(DATA,DATA,FFTW_PLANS(1))
              ENDIF
C
          IF ((FDEF.AND.CI).AND.(IVFLAG.EQ.1)) THEN
            IF (NONSYM) THEN
              WRITE(*,*)' ERROR: defocus refinement not',
     +                  ' possible with asymmetric refinement.'
              CALL FLUSH(6)
              STOP 'Fatal error'
            ENDIF
            RMAX1N=RMAX1(NSET)
            RMAX2N=RMAX2(NSET)
            CSN=CS(NSET)
            WLN=WL(NSET)
            IBUFB=IBUF
            IBUF=1
            ILST(IBUF)=K
C
            CALL CTFAPPLY(NSAM,SPEC,SHX(K),SHY(K),CS(NSET),
     +           WL(NSET),WGH1,WGH2,DFMID1(K),DFMID2(K),ANGAST(K),
     +           THETATR,CTFF,OUTD,OUTC,AMAGP(K),RIH,HALFW,
     +           RI2,RI3,RI4,DATD,DATC,B3DV,PHI(K),THETA(K),
     +           PSI(K),MBUF,XSTD,0,TX(NSET),TY(NSET),ASYM,
     +           PSSNR,BF,PWEIGHTS,1,FFTW_PLANS)
C
            CALL SIGMA(NSAM,IRADA,AMAGP(K),C3DF,DATC,
     +             PHI(K),THETA(K),PSI(K),DSHX,DSHY,SINCLUT,
     +             IPAD,IEWALD,THETATR,CTFF,BF,RI4,MAXR1,MAXR2,
     +             SIG(K),SIG2N,SIG2,LOGPC,DMASK,ASYM,K,SIGP,
     +             PWEIGHTS,RCLAS(NSET),FFTW_PLANS,SM)
            ALGP(K)=LOGPC+FSH(SIG2N,DSHX,DSHY,XM(NSET)-SHX(K),
     +                  YM(NSET)-SHY(K),SX(NSET),SY(NSET))
C
            CALL CTFREFINE(DFMID1,DFMID2,ANGAST,IFIRSTF,
     +        LASTF,FASTIG,NPAR,NSAM,MAXR1,MAXR2,PMASK,DDFMID1,
     +        DDFMID2,DANGAST,C3DF,IRADA,SPEC,
     +	      SHX,SHY,IPAD,CSN,WLN,WGH1,WGH2,THETATR,CTFF,
     +	      AMAGP,RIH,HALFWC,RI2,RI3,RI4,PHI,THETA,PSI,
     +        RMAX1N,RMAX2N,XSTD,MBUF,ILST,DATC,IBUF,B3DV,
     +        DATD,SINCLUT,IVFLAG,RBFACT,OUTD,OUTC,QBUC,CI,
     +        BF,FSCT,FILM,ILAST,IEWALD,TX(NSET),TY(NSET),
     +        SIG2N,DFSTD(NSET),ASYM,PWEIGHTS,PSSNR,FFTW_PLANS,SM)
            ANGAST(K)=ANGAST(K)-NINT(ANGAST(K)/PI)*pi
            WRITE(*,6014) K,DFMID1(K),DFMID2(K),ANGAST(K)/PI*180.0
6014        FORMAT(' Refined defocus for particle',I7,' = ',3F9.1)
            CALL FLUSH(6)
            IVFLAG=0
            IBUF=IBUFB
          ENDIF
C
          CALL CTFAPPLY(NSAM,SPEC,SHX(K),SHY(K),CS(NSET),
     +           WL(NSET),WGH1,WGH2,DFMID1(K),DFMID2(K),ANGAST(K),
     +           THETATR,CTFF,OUTD,OUTC,AMAGP(K),RIH,HALFW,
     +           RI2,RI3,RI4,DATD,DATC,B3DV,PHI(K),THETA(K),
     +           PSI(K),MBUF,XSTD,0,TX(NSET),TY(NSET),ASYM,
     +           PSSNR,BF,PWEIGHTS,1,FFTW_PLANS)
C
          IF (IFLAG.NE.0)
     +      CALL SIGMA(NSAM,IRADA,AMAGP(K),C3DF,DATC,
     +             PHI(K),THETA(K),PSI(K),DSHX,DSHY,SINCLUT,
     +             IPAD,IEWALD,THETATR,CTFF,BF,RI4,MAXR1,MAXR2,
     +             SIG(K),SIG2N,SIG2,LOGPC,DMASK,ASYM,K,SIGP,
     +             PWEIGHTS,RCLAS(NSET),FFTW_PLANS,SM)
C
          IF ((FDEF).OR.(FMAG).OR.(IFLAG.EQ.5)) THEN
            IF (ASYM(1:1).EQ.'H') THEN
C
            PRESA(K)=CC3M_C(NSAM,IRADA,AMAGP(K),OUTC,C3DF,
     +             MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,0,
     +             DUMC,RBFACT,SINCLUT,IPAD,BF,FSCT,IEWALD,
     +             THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,
     +             FFTW_PLANS,SM)
     +             +FSH(SIG2N,DSHX,DSHY,XM(NSET)-SHX(K),
     +                  YM(NSET)-SHY(K),SX(NSET),SY(NSET))
     +             +FANGLE(SIG2N,THETA(K),THETAMM(K),STHETA(NSET),
     +                     STIF)
     +             +FANGLE(SIG2N,PSI(K),PSIMM(K),SPSI(NSET),STIF)
C
            ELSE
C
            PRESA(K)=CC3M(NSAM,IRADA,AMAGP(K),OUTC,C3DF,
     +             MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,0,
     +             DUMC,RBFACT,SINCLUT,IPAD,BF,FSCT,IEWALD,
     +             THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,
     +             FFTW_PLANS,SM)
     +             +FSH(SIG2N,DSHX,DSHY,XM(NSET)-SHX(K),
     +                  YM(NSET)-SHY(K),SX(NSET),SY(NSET))
C
            ENDIF
C
          ENDIF
C
          IF (IFLAG.GE.3) THEN
C
            MAXR1=INT(NSAM*RMAX1(NSET)*ABS(AMAGP(K)))
            MAXR2=INT(NSAM*RMAX2(NSET)*ABS(AMAGP(K)))
            MAXRES=INT(NSAM*RREC(NSET)*ABS(AMAGP(K)))
            RBFACT=RBFAC(NSET)/(4.0*(NSAM*PSIZE*AMAGP(K))**2)
            BESTP=-1.0
            IF (IFLAG.EQ.5) BESTP=PRESA(K)
            JJ=MIN0(NSANG,IPMAX)
            IF (IFLAG.EQ.5) JJ=1
            TESTPAR(1,JJ+1)=PHI(K)
            TESTPAR(2,JJ+1)=THETA(K)
            TESTPAR(3,JJ+1)=PSI(K)
            TESTPAR(4,JJ+1)=0.0
            TESTPAR(5,JJ+1)=0.0
            DO 230 I=1,ITMAX
9       format(/' Entering PSEARCH particle',I6,
     .          ',   NSAM,IRADA,AMAGP,NSANG,I',2I5,F8.4,2I6)
C       write(*,9) K,NSAM,IRADA,AMAGP(K),NSANG,I
            CALL PSEARCH(NSAM,IRADA,AMAGP(K),OUTC,C3DF,DANGIN,
     +             MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,
     +             DSHY,PRESA(K),SANG,NSANG,CCD,CCC,CCBUF,IQUADMAX,
     +             PMASK,RBFACT,SINCLUT,IPAD,BF,IEWALD,THETATR,
     +             CTFF,RI2,RI3,RIH,HALFW,TESTPAR,JJ,ASYM,
     +             IRAN,FFTW_PLANS,SM)
C       write(*,8)-1,PHI(K)*180.0/PI,THETA(K)*180.0/PI,PSI(K)*180.0/PI,
C     .         DSHX*NSAM/PI/2.0,DSHY*NSAM/PI/2.0,PRESA(K)*100.0
C
            IF (IFLAG.EQ.5) THEN
              IF (ASYM(1:1).EQ.'H') THEN
C
                PRESA(K)=CC3M_C(NSAM,IRADA,AMAGP(K),OUTC,C3DF,
     +             MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,0,
     +             DUMC,RBFACT,SINCLUT,IPAD,BF,FSCT,IEWALD,
     +             THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,
     +             FFTW_PLANS,SM)
     +             +FSH(SIG2N,DSHX,DSHY,XM(NSET)-SHX(K),
     +                  YM(NSET)-SHY(K),SX(NSET),SY(NSET))
     +             +FANGLE(SIG2N,THETA(K),THETAMM(K),STHETA(NSET),
     +                     STIF)
     +             +FANGLE(SIG2N,PSI(K),PSIMM(K),SPSI(NSET),STIF)
C
              ELSE
C
                PRESA(K)=CC3M(NSAM,IRADA,AMAGP(K),OUTC,C3DF,
     +             MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,0,
     +             DUMC,RBFACT,SINCLUT,IPAD,BF,FSCT,IEWALD,
     +             THETATR,CTFF,RI2,RI3,RIH,HALFW,PWEIGHTS,
     +             FFTW_PLANS,SM)
     +             +FSH(SIG2N,DSHX,DSHY,XM(NSET)-SHX(K),
     +                  YM(NSET)-SHY(K),SX(NSET),SY(NSET))
C
              ENDIF
              IF (BESTP.GE.PRESA(K)) THEN
                PRESA(K)=BESTP
                PHI(K)=TESTPAR(1,JJ)
                THETA(K)=TESTPAR(2,JJ)
                PSI(K)=TESTPAR(3,JJ)
                DSHX=TESTPAR(4,JJ)
                DSHY=TESTPAR(5,JJ)
              ENDIF
              GOTO 250
            ENDIF
C
            PRESA(K)=-1.0
C
            JJ=MIN0(NSANG,IPMAX)
            IF (I.EQ.1) JJ=JJ+1
            DO 231 J=1,JJ
C       write(*,*) TESTPAR(1,J)*180.0/PI,TESTPAR(2,J)*180.0/PI,
C     .            TESTPAR(3,J)*180.0/PI,
C     .            TESTPAR(4,J)*NSAM/PI/2.0,
C     .            TESTPAR(5,J)*NSAM/PI/2.0,TESTPAR(6,J)
              CALL PREFINE(NSAM,IRADA,AMAGP(K),
     +           MAXR1,MAXR2,TESTPAR(1,J),TESTPAR(2,J),TESTPAR(3,J),
     +           TESTPAR(4,J),TESTPAR(5,J),PRESN,PMASK,RBFACT,NPAR,
     +		 C3DF,PBUF,SHX(K),SHY(K),IPAD,
     +		 CSN,WLN,WGH1,WGH2,THETATR,CTFF,AMAGP(K),RIH,
     +           HALFWC,RI2,RI3,RI4,RMAX1N,RMAX2N,XSTD,MBUF,
     +		 ILST,DATC,IBUF,B3DV,DATD,SINCLUT,IVFLAG,
     +           OUTD,OUTC,QBUC,BF,FSCT,IEWALD,XM(NSET),YM(NSET),
     +           SX(NSET),SY(NSET),SIG2N,ASYM,THETAMM(K),STHETA(NSET),
     +           PSIMM(K),SPSI(NSET),STIF,PWEIGHTS,FFTW_PLANS,SM)
              IF (PRESN.GT.PRESA(K)) THEN
                PRESA(K)=PRESN
                PHI(K)=TESTPAR(1,J)
                THETA(K)=TESTPAR(2,J)
                PSI(K)=TESTPAR(3,J)
                DSHX=TESTPAR(4,J)
                DSHY=TESTPAR(5,J)
              ENDIF
231         CONTINUE
C       write(*,8) 0,PHI(K)*180.0/PI,THETA(K)*180.0/PI,PSI(K)*180.0/PI,
C     .         DSHX*NSAM/PI/2.0,DSHY*NSAM/PI/2.0,PRESA(K)*100.0
8       format('  parameters at',I3,' are',5F9.3,F7.2)
C
            IF (IFLAG.EQ.3) GOTO 250
C
                IF(FHIST) THEN
                  CALL SIGMA(NSAM,IRADA,AMAGP(K),C3DF,DATC,
     +              PHI(K),THETA(K),PSI(K),DSHX,DSHY,SINCLUT,
     +              IPAD,IEWALD,THETATR,CTFF,BF,RI4,MAXR1,MAXR2,
     +              SIG(K),SIG2N,SIG2,LOGPC,DMASK,ASYM,K,SIGP,
     +              PWEIGHTS,RCLAS(NSET),FFTW_PLANS,SM)
                  ALGP(K)=LOGPC+FSH(SIG2N,DSHX,DSHY,XM(NSET)-SHX(K),
     +                  YM(NSET)-SHY(K),SX(NSET),SY(NSET))
                  PSX=PSIZE*ABS(AMAGP(K))
                  WRITE(NDOC1+NSET,7011)ABS(ILIST(K)),
     +            PSI(K)/PI*180.0,THETA(K)/PI*180.0,PHI(K)/PI*180.0,
     +            (SHX(K)+DSHX)*PSX*NSAM/PI/2.0,
     +            (SHY(K)+DSHY)*PSX*NSAM/PI/2.0,
     +            NINT(DSTEP(NSET)*10000.0/PSX),
     +            FILM(K),DFMID1(K),DFMID2(K),
     +            ANGAST(K)/PI*180.0,OCC(K)*100.0,
     +            NINT(ALGP(K)),SQRT(SIG2),PRESA(K)*100.0
                  CALL FLUSH(NDOC1+NSET)
                ENDIF
C
C               here for IFLAG=4, store best parameters
            IF (PRESA(K).GT.BESTP) THEN
                BESTP=PRESA(K)
                BESTPHI=PHI(K)
                BESTTHETA=THETA(K)
                BESTPSI=PSI(K)
                BESTDSHX=DSHX
                BESTDSHY=DSHY
                BESTAMAG=AMAGP(K)
                BESTDFMID1=DFMID1(K)
                BESTDFMID2=DFMID2(K)
                BESTANGAST=ANGAST(K)
            ENDIF
C criterion for ending search/refine depends on TARGET
            IF (PRESA(K).GE.ABS(TARGET(NSET))) GOTO 240
            PHI(K)=PHI(K)+RANDOM(IRAN)*DAMAX*REAL(PMASK(1))
            PSI(K)=PSI(K)+RANDOM(IRAN)*DAMAX*REAL(PMASK(3))
            THETA(K)=THETA(K)
     +        +(REMAP_THETA(RANDOM(IRAN)*DAMAX+THETA(K))-THETA(K))
     +         *REAL(PMASK(2))
C
230         CONTINUE
C               if IFLAG=4, do not set to -1.0 if score < THRESH, but restore
C               best params and best score
            PRESA(K)=BESTP
            PHI(K)=BESTPHI
            THETA(K)=BESTTHETA
            PSI(K)=BESTPSI
            DSHX=BESTDSHX
            DSHY=BESTDSHY
            AMAGP(K)=BESTAMAG
            DFMID1(K)=BESTDFMID1
            DFMID2(K)=BESTDFMID2
            ANGAST(K)=BESTANGAST
            GOTO 250
C
300         CONTINUE
            CALL PREFINE(NSAM,IRADA,AMAGP(K),
     +           MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,
     +           PRESA(K),PMASK,RBFACT,NPAR,
     +		 C3DF,PBUF,SHX(K),SHY(K),IPAD,
     +		 CSN,WLN,WGH1,WGH2,THETATR,CTFF,AMAGP(K),RIH,
     +           HALFWC,RI2,RI3,RI4,RMAX1N,RMAX2N,XSTD,MBUF,
     +		 ILST,DATC,IBUF,B3DV,DATD,SINCLUT,IVFLAG,
     +           OUTD,OUTC,QBUC,BF,FSCT,IEWALD,XM(NSET),YM(NSET),
     +           SX(NSET),SY(NSET),SIG2N,ASYM,THETAMM(K),STHETA(NSET),
     +           PSIMM(K),SPSI(NSET),STIF,PWEIGHTS,FFTW_PLANS,SM)
c       write(*,8) 1,PHI(K)*180./PI,THETA(K)*180./PI,PSI(K)*180./PI,
c     .         DSHX*NSAM/PI/2.0,DSHY*NSAM/PI/2.0,PRESA(K)*100.0
240         CONTINUE
C           extend resolution out to limit of reconstruction when IFLAG.eq.4
            MAXR2=MAXR2+1
            IF (MAXR2.LE.MAXRES) GOTO 300
C
250         CONTINUE
C
            IF (IFLAG.EQ.4) WRITE(*,6201)KK,PSIZE/RMAX1(NSET),
     +                          PSIZE/RMAX2(NSET),PRESA(K)*100.0
6201        FORMAT(' Best score for particle',I6,
     +             ' at Rmin/Rmax',2F6.1,':',F13.3)
                CALL FLUSH(6)
C accumulate overall score in resolution bins
                CALL PRESB(NSAM,IRADA,AMAGP(K),OUTC,C3DF,
     +            MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,RBIN,
     +            SINCLUT,IPAD,NPIX,BF(1),BF(NSAM+1),BF(2*NSAM+1),
     +            CS(NSET),WL(NSET),WGH1,WGH2,DFMID1(K),DFMID2(K),
     +            ANGAST(K),THETATR/AMAGP(K),IEWALD,CTFF,
     +            RI2,RI3,RIH,HALFW,ASYM,FFTW_PLANS,SM)
                DO 6221 I=1,NSAM
                        RAWBINS(I)=RAWBINS(I)+RBIN(I)
6221            CONTINUE
C
          ENDIF
C
	  IF ((IFLAG.EQ.1).OR.(IFLAG.EQ.2)) THEN

            FFLAG=0
            ITT=1
            IF (IFLAG.EQ.2) ITT=ITMAX
            PRESMIN=-1000.0
C
            DO 210 I=1,ITT
            MAXR1=INT(NSAM*RMAX1(NSET)*ABS(AMAGP(K)))
            MAXR2=INT(NSAM*RMAX2(NSET)*ABS(AMAGP(K)))
            RBFACT=RBFAC(NSET)/(4.0*(NSAM*PSIZE*AMAGP(K))**2)
200         CONTINUE
C       write(*,8)-2,PHI(K)*180.0/PI,THETA(K)*180.0/PI,PSI(K)*180.0/PI,
C     .         DSHX,DSHY,PRESA(K)*100.0
            CALL PREFINE(NSAM,IRADA,AMAGP(K),
     +           MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,
     +           PRESA(K),PMASK,RBFACT,NPAR,
     +	         C3DF,PBUF,SHX(K),SHY(K),IPAD,
     +		 CSN,WLN,WGH1,WGH2,THETATR,CTFF,AMAGP(K),RIH,
     +           HALFWC,RI2,RI3,RI4,RMAX1N,RMAX2N,XSTD,MBUF,
     +		 ILST,DATC,IBUF,B3DV,DATD,SINCLUT,IVFLAG,
     +           OUTD,OUTC,QBUC,BF,FSCT,IEWALD,XM(NSET),YM(NSET),
     +           SX(NSET),SY(NSET),SIG2N,ASYM,THETAMM(K),STHETA(NSET),
     +           PSIMM(K),SPSI(NSET),STIF,PWEIGHTS,FFTW_PLANS,SM)
C       write(*,8) 2,PHI(K)*180.0/PI,THETA(K)*180.0/PI,PSI(K)*180.0/PI,
C     .         DSHX,DSHY,PRESA(K)*100.0
C
            IF (PRESMIN.LT.PRESA(K)) THEN
              PRESMIN=PRESA(K)
              PHIM=PHI(K)
              THETAM=THETA(K)
              PSIM=PSI(K)
              DSHXM=DSHX
              DSHYM=DSHY
              AMAGM=AMAGP(K)
              IF (I.NE.1) FFLAG=1
            ELSE
              PRESA(K)=PRESMIN
              PHI(K)=PHIM
              THETA(K)=THETAM
              PSI(K)=PSIM
              DSHX=DSHXM
              DSHY=DSHYM
              AMAGP(K)=AMAGM
            ENDIF
            PHI(K)=PHI(K)+(RANDOM(IRAN)-0.5)*DAMAX*REAL(PMASK(1))
            THETA(K)=THETA(K)
     +        +(REMAP_THETA((RANDOM(IRAN)-0.5)*DAMAX+THETA(K))
     +        -THETA(K))*REAL(PMASK(2))
            PSI(K)=PSI(K)+(RANDOM(IRAN)-0.5)*DAMAX*REAL(PMASK(3))
            DSHX=DSHX+(RANDOM(IRAN)-0.5)*DSMAX*REAL(PMASK(4))
            DSHY=DSHY+(RANDOM(IRAN)-0.5)*DSMAX*REAL(PMASK(5))
C            AMAGP(K)=AMAGOLD
210         CONTINUE
C
            PRESA(K)=PRESMIN 
            PHI(K)=PHIM 
            THETA(K)=THETAM 
            PSI(K)=PSIM 
            DSHX=DSHXM 
            DSHY=DSHYM 
            AMAGP(K)=AMAGM 
C accumulate overall score in resolution bins
                CALL PRESB(NSAM,IRADA,AMAGP(K),OUTC,C3DF,
     +            MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,RBIN,
     +            SINCLUT,IPAD,NPIX,BF(1),BF(NSAM+1),BF(2*NSAM+1),
     +            CS(NSET),WL(NSET),WGH1,WGH2,DFMID1(K),DFMID2(K),
     +            ANGAST(K),THETATR/AMAGP(K),IEWALD,CTFF,
     +            RI2,RI3,RIH,HALFW,ASYM,FFTW_PLANS,SM)
                DO 221 I=1,NSAM
                        RAWBINS(I)=RAWBINS(I)+RBIN(I)
221             CONTINUE
C
C           Calculate cross-correlation
C
            CCPART=CCOEF(NSAM,IRADA,AMAGP(K),C3DF,DATC,
     +                PHI(K),THETA(K),PSI(K),DSHX,DSHY,SINCLUT,
     +                IPAD,IEWALD,THETATR,CTFF,BF,FFTW_PLANS,SM)
            WRITE(*,6602) KK,CCPART,PRESA(K)*100.0
6602        FORMAT(' CC for particle ',I6,' is',F12.8,
     .          '  score is',F7.2)
            CALL FLUSH(6)
            CCAVE=CCAVE+CCPART
            ICCOUNT=ICCOUNT+1
C
          ENDIF
C
C criterion for insertion depends on THRESH
          IF (((IFLAG.LE.3).AND.(PRESA(K).GE.ABS(THRESH(NSET)))).OR.
     +       ((IFLAG.EQ.4).AND.(PRESA(K).GT.-1.0))) THEN

            MAXR1=INT(NSAM*RMAX1(NSET)*ABS(AMAGP(K)))
            MAXR2=INT(NSAM*RMAX2(NSET)*ABS(AMAGP(K)))
            RBFACT=RBFAC(NSET)/(4.0*(NSAM*PSIZE*AMAGP(K))**2)

C         Put alternating particles into A3 or D3     

            IF ((IFLAG.EQ.0).AND.(ISTAT.EQ.1)) THEN
              CALL PRESB(NSAM,IRADA,AMAGP(K),OUTC,C3DF,
     +          MAXR1,MAXR2,PHI(K),THETA(K),PSI(K),DSHX,DSHY,RBIN,
     +          SINCLUT,IPAD,NPIX,BF(1),BF(NSAM+1),BF(2*NSAM+1),
     +          CS(NSET),WL(NSET),WGH1,WGH2,DFMID1(K),DFMID2(K),
     +          ANGAST(K),THETATR/AMAGP(K),IEWALD,CTFF,
     +          RI2,RI3,RIH,HALFW,ASYM,FFTW_PLANS,SM)
              DO 6222 I=1,NSAM
                RAWBINS(I)=RAWBINS(I)+RBIN(I)
6222          CONTINUE
            ENDIF

            IF (FALL) THEN
              IF (ISTAT.EQ.1)
     +        CALL VARIANCE(NSAM,SPEC,SHX(K),SHY(K),RI,
     +          BF(1),BF(1),BF(2*NSAM*NSAM+1),
     +          BF(2*NSAM*NSAM+1),IC,VSN,VN,
     +          VVSN,VVN,AMAGP(K),BF(4*NSAM*NSAM+1),NA)
              IF (IFLAG.EQ.0) THEN
                SIG2N=SIG(K)**2
                IF (SIG2N.EQ.0.0) SIG2N=1.0
              ENDIF
              CALL A3D3(NSET,ISWITCH,NSAM,IRAD,RI,RREC,K,ILAST,
     +          AMAGP(K),DATC,CTFF,ASUM,VSUM,PSUM,C3DF,
     +          KSUM,PHI(K),THETA(K),PSI(K),DSHX,DSHY,PRESA(K),PBC,
     +          BOFF,A3DF,S3DF,STD,D3DF,V3DF,VTD,
     +          NN1,MAXSET,SINCLUT,IRADA,IPAD,THETATR,IEWALD,
     +          NSYM,ISYMAX,JSYM,SYMOP,NNSTAT,NKSUM,IFSC,IMP,NVOL,
     +          NVOL1,NVOL2,OCC(K),SIG2N,PSIZE,ASYM,ALPHA,RISE,NU,
     +          HSTART,NKSUM1,TTD1,TTD2,NS,NS1,INTERP,RECBF,
     +          RECB,ISW,IRECB,CTFBF,NN2,NN3,SM,NONSYM,ISM)
            ENDIF
            PALL=PALL+PRESA(K)*OCC(K)
            OALL=OALL+OCC(K)
            IF (OCC(K).NE.0.0) THEN
              WALL=WALL+OCC(K)/SIG2N
              IWALL=IWALL+1
            ENDIF
            IALL=IALL+1
          ELSE
            IF (FALL) WRITE(*,6600)K
6600        FORMAT(' **** PARTICLE',I10,' BELOW THRESHOLD,',
     +             ' NOT INCLUDED')
            CALL FLUSH(6)
            IEXCL=IEXCL+1
          ENDIF

C        ELSE
C          WRITE(*,6600)K
C          CALL FLUSH(6)
C          IEXCL=IEXCL+1
C        ENDIF

C      	  IF (KK.EQ.LASTF) IBUF=0
C
      	  IF (IFLAG.GT.0.OR.FDEF.OR.FMAG) THEN
C
            CALL SIGMA(NSAM,IRADA,AMAGP(K),C3DF,DATC,
     +             PHI(K),THETA(K),PSI(K),DSHX,DSHY,SINCLUT,
     +             IPAD,IEWALD,THETATR,CTFF,BF,RI4,MAXR1,MAXR2,
     +             SIG(K),SIG2N,SIG2,LOGPC,DMASK,ASYM,K,SIGP,
     +             PWEIGHTS,RCLAS(NSET),FFTW_PLANS,SM)
            ALGP(K)=LOGPC+FSH(SIG2N,DSHX,DSHY,XM(NSET)-SHX(K),
     +                  YM(NSET)-SHY(K),SX(NSET),SY(NSET))
            ALGPS=ALGPS+ALGP(K)*100.0*OCC(K)
            PRESS=PRESS+PRESA(K)*100.0*OCC(K)
      	    IF (PSI(K).LT.0.0) PSI(K)=PSI(K)+2.0*PI
      	    IF (THETA(K).LT.0.0) THETA(K)=THETA(K)+2.0*PI
      	    IF (PHI(K).LT.0.0) PHI(K)=PHI(K)+2.0*PI
            ANGAST(K)=ANGAST(K)-NINT(ANGAST(K)/PI)*pi
            PSX=PSIZE*ABS(AMAGP(K))
            IF (PRESA(K).GT.2.0) PRESA(K)=2.0
            IF (PRESA(K).LT.0.0) PRESA(K)=0.0
       	    WRITE(NDOC1+NSET,7011)ABS(ILIST(K)),PSI(K)/PI*180.0,
     +		THETA(K)/PI*180.0,PHI(K)/PI*180.0,
     +		(SHX(K)+DSHX)*PSX*NSAM/PI/2.0,
     +          (SHY(K)+DSHY)*PSX*NSAM/PI/2.0,
     +          NINT(DSTEP(NSET)*10000.0/PSX),FILM(K),DFMID1(K),
     +          DFMID2(K),ANGAST(K)/PI*180.0,OCC(K)*100.0,
     +          NINT(ALGP(K)),SQRT(SIG2),PRESA(K)*100.0,
     +          (PRESA(K)-PRESAOLD)*100.0
	    CALL FLUSH(NDOC1+NSET)
7011        FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I10,
     +             F11.4,2F8.2)
            DPSI=PSI(K)-PSIOLD
     +             -INT((INT((PSI(K)-PSIOLD)/PI)+1)/2.0)*2.0*PI
            DTHETA=THETA(K)-THETAOLD
     +             -INT((INT((THETA(K)-THETAOLD)/PI)+1)/2.0)*2.0*PI
            DPHI=PHI(K)-PHIOLD
     +             -INT((INT((PHI(K)-PHIOLD)/PI)+1)/2.0)*2.0*PI
            DANGAST=DANGAST
     +             -INT((INT(DANGAST/PI)+1)/2.0)*2.0*PI
            WRITE(NDOC2+NSET,7011)ABS(ILIST(K)),
     +        DPSI/PI*180.0,DTHETA/PI*180.0,DPHI/PI*180.0,
     +        DSHX*PSX*NSAM/PI/2.0,DSHY*PSX*NSAM/PI/2.0,
     +        NINT((1.0/ABS(AMAGP(K))-1.0/ABS(AMAGOLD))
     +          *DSTEP(NSET)*10000.0/PSIZE),
     +        FILM(K),DDFMID1,DDFMID2,DANGAST/PI*180.0,
     +        0.0,NINT(ALGP(K)-ALGPOLD),SQRT(SIG2)-SIGOLD,
     +        (PRESA(K)-PRESAOLD)*100.0
	    CALL FLUSH(NDOC2+NSET)
      	    PSIS=PSIS+(DPSI/PI*180.0)**2
      	    THETAS=THETAS+(DTHETA/PI*180.0)**2
      	    PHIS=PHIS+(DPHI/PI*180.0)**2
      	    DSHXS=DSHXS+(DSHX*NSAM/PI/2.0)**2
      	    DSHYS=DSHYS+(DSHY*NSAM/PI/2.0)**2
      	    ABSMAGS=ABSMAGS+(ABS(AMAGP(K))-ABS(AMAGOLD))**2
     .			*(DSTEP(NSET)*10000./PSIZE)**2
      	    DFMID1S=DFMID1S+DDFMID1**2
      	    DFMID2S=DFMID2S+DDFMID2**2
      	    ANGASTS=ANGASTS+(DANGAST/PI*180.0)**2
            OCCS=OCCS+OCC(K)*100.0
            DALGPS=DALGPS+(ALGP(K)-ALGPOLD)*100.0*OCC(K)
            DPRESS=DPRESS+(PRESA(K)-PRESAOLD)*100.0*OCC(K)
      	    ICOUNT=ICOUNT+1
      	    IF (FFLAG.EQ.1) WRITE(*,6400)KK,(PSI(K)-PSIOLD)/PI*180.0,
     +		(THETA(K)-THETAOLD)/PI*180.0,(PHI(K)-PHIOLD)/PI*180.0,
     +		DSHX*NSAM/PI/2.0,DSHY*NSAM/PI/2.0
6400  	    FORMAT(' **** NEW LOCAL MINIMUM FOR PROJECTION',I6,/
     +      ' dPSI, dTHETA, dPHI, dX, dY = ',3F8.2,3X,2F6.2)
      	    CALL FLUSH(6)
C
          ENDIF
C
C         Output of matching projections if requested and MASKENV
          IF (FMATCH) THEN
	    CALL MATCH(K,NSAM,IRADA,OUTC,OUTD,C3DF,
     +               PHI(K),THETA(K),PSI(K),JC,NSAMH,PHASE,SHX(K),
     +               SHY(K),DSHX,DSHY,XSTD,RI,B3DV,MBUF,ICMP,DATA,
     +               SPEC,CCPART,PRESA(K),PRESAOLD,THRESH,
     +               CFORM,MAXSET,ILIST(K),IFIRST,NSET,IOPROJ,
     +		     SINCLUT,IPAD,THETATR,IEWALD,CTFF,PSIZE,
     +               FFTW_PLANS,SM,DMASK)
      	  ENDIF
C
9999	RETURN
	END
