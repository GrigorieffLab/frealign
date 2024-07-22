C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C------------------------------------------------------------------------------
	PROGRAM FREALIGN
C------------------------------------------------------------------------------
C
#ifdef _NAG
        USE F90_UNIX
#endif
C
        USE ISO_C_BINDING
        USE FFTW33
C
	PARAMETER  (NSCF=50)
      	PARAMETER  (MAXSET=20,NDOC1=30,NDOC2=NDOC1+MAXSET)
CNIKO       MAXSET=100 does not work on some machines....  set back to 20
      	PARAMETER  (PI=3.1415926535897)
        PARAMETER  (HALFW=6.0)
        PARAMETER  (ITMAXA=10)
        PARAMETER  (DAMAX=360.0*PI/180.0)
      	PARAMETER  (ISYMAX=100)
      	PARAMETER  (IPMAXX=200)
      	PARAMETER  (MSIZE=128)
        PARAMETER  (DALT=0.81)

        CHARACTER*200 SFILE,FINPAT1(MAXSET),FINPAT2(MAXSET)
        CHARACTER*200 FINPAR,F3D,FWEIGH,F3D1,F3D2,FPHA,FPOI
        CHARACTER*200 FOUTPAR
        CHARACTER CFORM*1,ASYMTEST*16,ASYM*3,VX*15
        CHARACTER CDATE*8,CTIME*10,CZONE*5,NCPUS*10

        INTEGER NSANG,NNSET(MAXSET),JSYM(ISYMAX),TNSET(MAXSET)
        INTEGER PMASK(5),DTVAL(8),IPAD,ISTAT,IEWALD,IRAD,IRADA
        INTEGER IREDUN,IERR,IBLOW,HSTART,NNSTAT,NVOL1,NVOL2,NVOL
        INTEGER*8 NKSUM
        INTEGER OMP_GET_NUM_PROCS,ISYM(ISYMAX)

        INTEGER,ALLOCATABLE :: NPIX(:),NPIXT(:),NS(:),IC(:)
        INTEGER,ALLOCATABLE :: KSUM(:,:),ILST(:),NS1(:)
        INTEGER,ALLOCATABLE :: FILM(:),ILIST(:),ISW(:)
        INTEGER*8,ALLOCATABLE :: NKSUM1(:)

        LOGICAL FMAG,FDEF,FPART,FMATCH,ASYMSTORE,FHIST,FALL,FBOOST
        LOGICAL FASTIG,FBEAUT,CI(MAXSET),FCREF,FBFACT,FSTAT,FDUMP
        LOGICAL NONSYM

        REAL TX(MAXSET),TY(MAXSET),DFSTD(MAXSET),DMASK(4),SM(9)
        REAL RELMAG(MAXSET),CS(MAXSET),WL(MAXSET),DSTEP(MAXSET)
        REAL RREC(MAXSET),RMAX1(MAXSET),RMAX2(MAXSET),RBFAC(MAXSET)
        REAL TARGET(MAXSET),THRESH(MAXSET),STHETA(MAXSET)
        REAL SYMOP(3,3,ISYMAX),SINCLUT(2000),SPSI(MAXSET)
        REAL CMIN(MAXSET),CMAX(MAXSET),XM(MAXSET),YM(MAXSET)
        REAL SX(MAXSET),SY(MAXSET),QM(MAXSET),RCLAS(MAXSET)

        REAL,ALLOCATABLE :: RBIN(:),SANG(:,:)
        REAL,ALLOCATABLE :: MBUF(:),PBUF(:)
        REAL,ALLOCATABLE :: S3DF(:,:),BF(:)
        REAL,ALLOCATABLE :: V3DF(:,:),SSNR(:)
        REAL,ALLOCATABLE :: QCP(:),PWEIGHTS(:)
        REAL,ALLOCATABLE :: ASUM(:,:),PSIM(:)
        REAL,ALLOCATABLE :: VSUM(:,:),THETAM(:)
        REAL,ALLOCATABLE :: PSUM(:,:),OCC(:),SIG(:)
        REAL,ALLOCATABLE :: TESTPAR(:,:)
        REAL,ALLOCATABLE :: DATA(:),OUTD(:)
        REAL,ALLOCATABLE :: DATD(:),CCD(:)
        REAL,ALLOCATABLE :: RBINS(:),RBINN(:),RAWBINS(:)
        REAL,ALLOCATABLE :: PR(:),FSC(:),ACC(:),PSSNR(:)
        REAL,ALLOCATABLE :: QF(:),PF(:),QC(:),RBUF(:)
        REAL,ALLOCATABLE :: PSI(:),THETA(:),PHI(:)
        REAL,ALLOCATABLE :: AMAGP(:),PRESA(:),ANGAST(:)
        REAL,ALLOCATABLE :: DFMID1(:),DFMID2(:),SHX(:),SHY(:)
        REAL,ALLOCATABLE :: FSCT(:),ALGP(:),P(:)
        REAL,ALLOCATABLE :: REF3DV(:),RECB(:,:)

        DOUBLEPRECISION STD,VTD

        DOUBLEPRECISION,ALLOCATABLE :: VSN(:),VN(:),TTD1(:)
        DOUBLEPRECISION,ALLOCATABLE :: VVSN(:),VVN(:),TTD2(:)

        COMPLEX,ALLOCATABLE :: A3DF(:,:),B3DF(:),C3DF(:)
        COMPLEX,ALLOCATABLE :: D3DF(:,:),RECBF(:,:)
        COMPLEX,ALLOCATABLE :: CTFBF(:,:),CBUF(:),QBUC(:)
        COMPLEX,ALLOCATABLE :: CCBUF(:),CTFF(:)

        TYPE(C_PTR) FFTW_PLANS(10)
C  FFTW_PLANS(1),FFTW_PLANS(2)  = NSAM x NSAM           transforms, fwd,bwd
C  FFTW_PLANS(3),FFTW_PLANS(4)  = NSAM x NSAM x NSAM    transforms, fwd,bwd
C  FFTW_PLANS(5),FFTW_PLANS(6)  = NNBIG x NNBIG x NNBIG transforms, fwd,bwd (SUBROUTINE CARDS13AND14)
C  FFTW_PLANS(7),FFTW_PLANS(8)  = NSAMR x NSAMR         transforms, fwd,bwd (SUBROUTINE CCP)
C  FFTW_PLANS(9),FFTW_PLANS(10) = NSAM4 x NSAM4         transforms, fwd,bwd (SUBROUTINE PAD)

	DATA  IPOINT/3/
	DATA  I3D1/1/
	DATA  I3D2/2/
	DATA  IPSTAT/10/
	DATA  INSTAT/9/
	DATA  INPIC/8/
	DATA  INPROJ/4/
	DATA  IOPROJ/7/
      	DATA  FALL/.TRUE./
      	DATA  FHIST/.FALSE./
      	DATA  ASYMTEST/' CDTOI0123456789'/
      	DATA  VX/'9.12 - 19.07.24'/
C       15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C------------------------------------------------------------------------------
        IMP=0
#ifdef _OPENMP
        CALL GETENV('OMP_NUM_THREADS',NCPUS)
        READ(NCPUS,*,ERR=111,END=111)IMP
111     CONTINUE
        IF (IMP.LE.0) THEN
          CALL GETENV('NCPUS',NCPUS)
          READ(NCPUS,*,ERR=112,END=112)IMP
112       CONTINUE
        ENDIF
        IF (IMP.LE.0) THEN
          IMP=OMP_GET_NUM_PROCS()
        ENDIF
#endif
        IF (IMP.LE.0) IMP=1
C
	CALL CALCSINC(SINCLUT,2000)
	CALL CARD1(VX,CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,
     +		FMATCH,FSTAT,FBEAUT,FCREF,FBFACT,IFSC,IMP,
     +          INTERP,IMEM,FDUMP,FBOOST)
        IBLOW=1
        IF ((IMEM.EQ.1).OR.(IMEM.EQ.3)) IBLOW=4
        NNSTAT=0
        IF (FSTAT) NNSTAT=1
        CALL FLUSH(6)
	CALL CARD2(RI,RIC,PSIZE,MW,WGH,XSTD,PBC,BOFF,DANGIN,IPMAX,
     +             ITMAX,ITMAXA,IPMAXX)
        CALL FLUSH(6)
	CALL CARD3(PMASK,DMASK)
        IF (DMASK(1).GE.0.0) THEN
          DMASK(1)=DMASK(1)/PSIZE
          DMASK(2)=DMASK(2)/PSIZE
          DMASK(3)=DMASK(3)/PSIZE
          DMASK(4)=DMASK(4)/PSIZE
        ENDIF
        CALL FLUSH(6)
	CALL CARD4(IFIRST,ILAST)
        CALL FLUSH(6)
        CALL CARD5(ASYMTEST,ASYM,JSYM,IREDUN,SYMOP,ASYMSTORE,
     +        NSYM,NASYM,ISYMAX,ALPHA,RISE,NU,HSTART,STIF,NONSYM)
        IF (NONSYM) THEN
          IF (FMATCH) THEN
            WRITE(*,*) ' ERROR: cannot generate matching',
     +        ' projections with asymmetric refinement.'
            CALL FLUSH(6)
            STOP 'Fatal error'
          ENDIF
          IF (FMAG) THEN
            WRITE(*,*)' ERROR: magnification refinement not',
     +                ' possible with asymmetric refinement.'
            CALL FLUSH(6)
            STOP 'Fatal error'
          ENDIF
          IF (FDEF) THEN
            WRITE(*,*)' ERROR: defocus refinement not',
     +                ' possible with asymmetric refinement.'
            CALL FLUSH(6)
            STOP 'Fatal error'
          ENDIF
          IF (FBEAUT) THEN
            WRITE(*,*)' ERROR: the beautify option is',
     +                ' incompatible with asymmetric refinement.'
            CALL FLUSH(6)
            STOP 'Fatal error'
          ENDIF
          NNPART=ILAST*IREDUN
        ELSE
          NNPART=ILAST
        ENDIF
        CALL FLUSH(6)
C-----------------------------------------------------------------------------------
      	HALFWC=HALFW
      	WGH1=SQRT(1.0-WGH**2)
      	WGH2=WGH
      	WC2=0.1
      	IRAD=1
      	AMEANPHAS=0.0
      	NPHAS=0
      	RLIM=0.0
	NANG=0			! NANG - number of angles (projections)
      	NSET=1
      	NSETST=0
      	NSS=1
      	RMM=0.0
        NPB=1
        BFACT=0.0
        BMASK=0.0
        IRECB=0
c----------------------------------------------------------------------------------
1999    CALL CARD6(FALL,NSET,MAXSET,RELMAG,DSTEP,TARGET,THRESH,CS,
     .	           AKV,WL,TX,TY)
        CALL FLUSH(6)
        IF (RELMAG(NSET).EQ.0.0.OR.(.NOT.FALL)) GOTO 1998
C
        CALL  CARD7(NSET,MAXSET,RREC,RMAX1,RMAX2,RCLAS,DFSTD,RBFAC)
        CALL FLUSH(6)
        CI(NSET)=FPART
C
        CALL CARDS8AND9(NSET,MAXSET,FINPAT1,FINPAT2,INPROJ,
     .                  CFORM,ASYM,PSIZE,VX,NSAM,ICMP,RLIM,
     .                  RMM,RMAX1,RMAX2,RREC,ABS(RI),MSIZE,RCLAS)
#ifdef _OPENMP
        IF ((NSAM/8.LT.IMP).AND.(IMEM.LT.2)) THEN
          IMP=NSAM/8
          WRITE(*,6801) IMP
6801      FORMAT(/' Number of parallel processes too large.',
     .            ' Resetting IMP = ',I5/)
        ENDIF
        CALL OMP_SET_NUM_THREADS(IMP)
        CALL FLUSH(6)
#endif
C
C++++++++ ALLOCATE MEMORY START ++++++++++++++++++++++++++++++++++++++++++++++
C
        IF (NSET.EQ.1) THEN
C
        NN1=NSAM
        NVOL=1
        NVOL1=1
        NVOL2=1
        IF (IMEM.GE.2) THEN
          NVOL=IMP
          IF ((FALL).AND.(NVOL.GT.1)) WRITE(*,1100) NVOL
1100      FORMAT(/,' Parallelized reconstruction using NVOL =',I4,/)
          NVOL1=IMP
          IF (IFSC.LE.0) THEN
            NVOL1=INT(NVOL1/2.0+0.5)
            IF (NVOL1.EQ.0) NVOL1=1
            NVOL2=NVOL1
          ENDIF
        ENDIF
C
        ALLOCATE(NPIX(NN1),NPIXT(NN1/2+1),NS(NN1),
     +           NS1(NN1/2*IMP+IMP),FILM(NNPART+IREDUN),
     +           ILIST(NNPART+IREDUN),STAT=IERR)
        IF (IERR.NE.0)
     +    STOP ' ERROR: Memory allocation failed for INTEGER arrays'
C
        ALLOCATE(RBIN(NN1*2),RBUF(NN1),BF(25*NN1*NN1+50*NN1),
     +           ALGP(NNPART+IREDUN),PSIM(NNPART+IREDUN),
     +           TESTPAR(6,IPMAXX+1),THETAM(NNPART+IREDUN),
     +           SIG(NNPART+IREDUN),OCC(NNPART+IREDUN),
     +           DATA(NN1*NN1+2*NN1),OUTD(2*NN1*NN1+4*NN1),
     +           DATD(NN1*NN1+2*NN1),CCD(NN1*NN1+2*NN1),
     +           RBINS(NN1*2),RBINN(NN1*2),RAWBINS(NN1*2),
     +           PSI(NNPART+IREDUN),THETA(NNPART+IREDUN),
     +           PHI(NNPART+IREDUN),ANGAST(NNPART+IREDUN),
     +           AMAGP(NNPART+IREDUN),PRESA(NNPART+IREDUN),
     +           DFMID1(NNPART+IREDUN),DFMID2(NNPART+IREDUN),
     +           SHX(NNPART+IREDUN),SHY(NNPART+IREDUN),FSCT(NN1/2+1),
     +           PSSNR(NN1/2+1),PWEIGHTS(NN1/2+1),ISW(NVOL),
     +           RECBF(NN1/2*NN1+NN1,NVOL),CTFBF(NN1*NN1+2*NN1,NVOL),
     +           RECB(10,NVOL),P(IREDUN),STAT=IERR)
        IF (IERR.NE.0) THEN
          WRITE(*,*) ' ERROR: Memory allocation failed for REAL arrays'
          STOP ' Try setting FSTAT=F'
        ENDIF
C
        NNBIG=NN1*IBLOW
        ALLOCATE(C3DF(NNBIG/2*NNBIG*NNBIG+NNBIG*NNBIG),
     +           CBUF(NN1/2+1),CCBUF(NN1*NN1+2*NN1),
     +           CTFF(2*NN1*NN1+2*NN1),STAT=IERR)
        IF (IERR.NE.0) THEN
          WRITE(*,*) ' ERROR: Memory allocation failed for',
     +               ' COMPLEX arrays'
          STOP ' Try setting IMEM=0'
        ENDIF
C
        CALL FFTW_PLANS_2D(NN1,DATA,DATA,FFTW_PLANS(1),FFTW_PLANS(2))
        NSAMRH=NN1/4
        IF (2*NSAMRH.NE.NN1/2) NSAMRH=(NN1/2+1)/2
        NSAMR=2*NSAMRH
        CALL FFTW_PLANS_2D(NSAMR,OUTD,OUTD,FFTW_PLANS(7),FFTW_PLANS(8))
        CALL FFTW_PLANS_2D(NN1*4,BF,BF,FFTW_PLANS(9),FFTW_PLANS(10))
C
        IPAD=INT(NNBIG/NSAM)
        IF (IPAD.GT.4) IPAD=4
        IF (IPAD.EQ.3) IPAD=4
        IF (IPAD.EQ.2) IPAD=1
        CALL FFTW_PLANS_3D(NN1*IPAD,C3DF,C3DF,FFTW_PLANS(5),
     +                     FFTW_PLANS(6))
C
C       determine the need for interpolation with radius IRADA
C       based on particle radius RI. To be safe, pretend that 3*RI
C       is particle diameter. If volume dimension IPAD*NSAM is 4 times
C       that size, IRADA=0 (one point interpolation), if IPAD*NSAM is
C       2 times the particle diameter, IRADA=1 (3 x 3 x 3 points 
C       interpolation), if IPAD*NSAM still smaller, than use IRADA=2
C       (5 x 5 x 5 points interpolation).
        I=INT(IPAD*NSAM/(3*ABS(RI)))
        IF (I.LE.1) IRADA=1
        IF (I.GE.2) IRADA=1
        IF (I.GE.3) IRADA=0
        IF (IPAD.NE.1)
     .    WRITE (*,*) 'Padding reference volume, IPAD =',IPAD
        WRITE (*,*) 'Interpolation radius for reference ',
     .    'IRADA =',IRADA
        WRITE(*,*)
        CALL FLUSH(6)

        IF (NSET.EQ.1) THEN
          DO 81 I=1,NSAM/2
            FSCT(I)=0.0
81        CONTINUE
        ENDIF
C
        ENDIF
C
        CALL CARD10(NANG,IFLAG,FINPAR,NNPART,ILIST,PSI,THETA,PHI,
     +              SHX,SHY,FILM,DFMID1,DFMID2,ANGAST,OCC,PRESA,
     +              AMEANPHAS,NPHAS,NSAM,AMAGP,NSET,PSIZE,MAXSET,
     +              RELMAG,NSS,IFIRST,ILAST,NNSET,NSETST,DSTEP,TNSET,
     +              FSCT,THRESH,CMIN,CMAX,XM,YM,SX,SY,NPB,
     +              ASYM,THETAM,STHETA,PSIM,SPSI,ALPHA,RISE,QM,ALGP,
     +              SIG,ABS(RI),PSSNR,FBOOST,IREDUN,NONSYM,IRAN,P)
        IF (NPB*NN1*NN1.GT.1024*1024*256) NPB=1024*1024*256/NN1/NN1
        IF (FMAG.OR.(FDEF.AND.(.NOT.FPART))) THEN
C
          ALLOCATE(ILST(NPB),MBUF(NPB*NN1*NN1+2*NPB*NN1),
     +             PBUF(NPB*NN1*NN1+2*NPB*NN1),
     +             QBUC(NPB*NN1/2*NN1+NPB*NN1),STAT=IERR)
          IF (IERR.NE.0)
     +    STOP ' ERROR: Memory allocation failed for BUFFER arrays'
        ELSEIF (FDEF) THEN
          ALLOCATE(ILST(1),MBUF(NN1*NN1+2*NPB*NN1),
     +             PBUF(NN1*NN1+2*NPB*NN1),
     +             QBUC(NN1/2*NN1+NN1),STAT=IERR)
          IF (IERR.NE.0)
     +    STOP ' ERROR: Memory allocation failed for BUFFER arrays'
        ELSEIF (XSTD.LT.0.0) THEN
          ALLOCATE(ILST(1),MBUF(NN1*NN1+2*NN1),STAT=IERR)
          IF (IERR.NE.0)
     +    STOP ' ERROR: Memory allocation failed for BUFFER arrays'
        ENDIF
        CALL FLUSH(6)
        IF (DMASK(1).LT.0.0) THEN
          DMASK(1)=NSAM/2
          DMASK(2)=NSAM/2
          DMASK(3)=NSAM/2
          DMASK(4)=ABS(RI)
        ENDIF
	CALL CARDS11AND12(NDOC1,NDOC2,NSET,CDATE,CTIME,CZONE,DTVAL,VX,
     .		    CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FMATCH,
     .		    FDUMP,FBEAUT,FCREF,FBFACT,RI,RIC,PSIZE,WGH,XSTD,PBC,
     .              BOFF,ASYM,NSYM,SYMOP,IPMAX,ITMAX,DANGIN,IFIRST,
     .	            ILAST,NSS,MAXSET,RELMAG,DSTEP,TARGET,THRESH,CS,AKV,
     .              RREC,RMAX1,RMAX2,RBFAC,FINPAT1,FINPAT2,FINPAR,
     .              ISYMAX,TX,TY,IFSC,DFSTD,IMEM,ALPHA,RISE,NU,HSTART,
     .              STIF,PMASK,DMASK,MW,INTERP,RCLAS,FOUTPAR)
        DMASK(1)=DMASK(1)+1.0
        DMASK(2)=DMASK(2)+1.0
        DMASK(3)=DMASK(3)+1.0
        CALL FLUSH(6)
C
      	GOTO 1999
1998	CONTINUE
C   End of input data cards - only map files are read in after this
C
        IF (FALL) THEN
          NN2=NN1*NN1
          NN3=NN1/2*NN1*(NN1+2)
C
C  Arrays used for FFTs: KSUM
C
          ALLOCATE(KSUM(NN3*NNSTAT+1,NVOL),IC(NN1),STAT=IERR)
          IF (IERR.NE.0)
     +      STOP ' ERROR: Memory allocation failed for INTEGER arrays'
        ELSE
          NN2=1
          NN3=1
          ALLOCATE(KSUM(NN3*NNSTAT+1,NVOL),STAT=IERR)
          IF (IERR.NE.0)
     +      STOP ' ERROR: Memory allocation failed for INTEGER arrays'
        ENDIF
C
        IF (FALL) THEN
C
          ALLOCATE(S3DF(NN3,NVOL1),
     +             SSNR(NN1/2*NN1*NN1*NNSTAT+NN1*NN1*NNSTAT+1),
     +             QCP(NN1/2*NN1*NN1*NNSTAT+NN1*NN1*NNSTAT+1),
     +             NKSUM1(IMP),
     +             ASUM(NN3*NNSTAT+1,NVOL),TTD1(IMP),
     +             VSUM(NN3*NNSTAT+1,NVOL),
     +             PSUM(NN3*NNSTAT+1,NVOL),
     +             REF3DV(NN1*NN1*NN1*NNSTAT+2*NN1*NN1*NNSTAT+2),
     +             STAT=IERR)
          IF (IERR.NE.0) THEN
            WRITE(*,*) ' ERROR: Memory allocation failed',
     +                 ' for REAL arrays'
            STOP ' Try setting FSTAT=F and/or IMEM=0'
          ENDIF
          IF (IFSC.LE.0) THEN
C
            ALLOCATE(V3DF(NN3,NVOL2),TTD2(IMP),STAT=IERR)
            IF (IERR.NE.0) THEN
              WRITE(*,*) ' ERROR: Memory allocation failed',
     +                 ' for REAL arrays'
              STOP ' Try setting IFSC=3'
            ENDIF
          ENDIF
          IF ((IFSC.LE.0).OR.(NNSTAT.NE.0)) THEN
            ALLOCATE(PR(NN1/2+1),FSC(NN1/2+1),ACC(NN1/2+1),
     +        QF(NN1/2+1),PF(NN1/2+1),QC(NN1/2+1),STAT=IERR)
            IF (IERR.NE.0) THEN
              WRITE(*,*) ' ERROR: Memory allocation failed',
     +                 ' for REAL arrays'
              STOP ' Try setting IFSC=3 or NNSTAT=0'
            ENDIF
          ENDIF
        ELSE
          ALLOCATE(S3DF(NN3,NVOL1),
     +             ASUM(NN3*NNSTAT+1,NVOL),
     +             VSUM(NN3*NNSTAT+1,NVOL),
     +             PSUM(NN3*NNSTAT+1,NVOL),
     +             STAT=IERR)
          IF (IERR.NE.0) THEN
            WRITE(*,*) ' ERROR: Memory allocation failed',
     +                 ' for REAL arrays'
            STOP
          ENDIF
          IF (IFSC.LE.0) THEN
            ALLOCATE(V3DF(NN3,NVOL2),STAT=IERR)
            IF (IERR.NE.0) THEN
              WRITE(*,*) ' ERROR: Memory allocation failed',
     +                 ' for REAL arrays'
              STOP ' Try setting IFSC=3'
            ENDIF
          ENDIF
        ENDIF
C
        IF (FALL) THEN
          ALLOCATE(VSN(NN1),VN(NN1),VVSN(NN1),VVN(NN1),STAT=IERR)
          IF (IERR.NE.0)
     +    STOP ' ERROR: Memory allocation failed for DPRECISION arrays'
        ENDIF
C
        IF (FALL) THEN
C
          ALLOCATE(A3DF(NN3,NVOL1),STAT=IERR)
          IF (IERR.NE.0)
     +    STOP ' ERROR: Memory allocation failed for COMPLEX arrays'
          IF (IFSC.LE.0) THEN
C
            ALLOCATE(B3DF(NN1/2*NN1*NN1+NN1*NN1),
     +               D3DF(NN3,NVOL2),STAT=IERR)
            IF (IERR.NE.0) THEN
              WRITE(*,*) ' ERROR: Memory allocation failed for',
     +                   ' COMPLEX arrays'
              STOP ' Try setting IFSC=3'
            ENDIF
          ELSEIF (FBEAUT) THEN
            ALLOCATE(B3DF(NN1/2*NN1*NN1+NN1*NN1),STAT=IERR)
            IF (IERR.NE.0) THEN
              WRITE(*,*) ' ERROR: Memory allocation failed for',
     +                   ' COMPLEX arrays used in BEAUTIFY'
              STOP ' Try setting FBEAUT=F'
            ENDIF
          ENDIF
        ELSEIF ((IPAD.NE.1).OR.(XSTD.NE.0.0)) THEN
            ALLOCATE(B3DF(NN1/2*NN1*NN1+NN1*NN1),STAT=IERR)
            IF (IERR.NE.0) THEN
              WRITE(*,*) ' ERROR: Memory allocation failed for',
     +                   ' COMPLEX arrays'
              STOP ' Try setting IMEM=0'
            ENDIF
        ENDIF
        IF (.NOT.FALL) THEN
          ALLOCATE(A3DF(NN3,NVOL1),STAT=IERR)
          IF (IERR.NE.0)
     +    STOP ' ERROR: Memory allocation failed for COMPLEX arrays'
          IF (IFSC.LE.0) THEN
            ALLOCATE(D3DF(NN3,NVOL2),STAT=IERR)
            IF (IERR.NE.0) THEN
              WRITE(*,*) ' ERROR: Memory allocation failed for',
     +                   ' COMPLEX arrays'
              STOP
            ENDIF
          ENDIF
        ENDIF
C
        CALL FFTW_PLANS_3D(NN1,A3DF,A3DF,FFTW_PLANS(3),FFTW_PLANS(4))
C
C----------------------------------------------------------------------------------
C++++++++ ALLOCATE MEMORY END ++++++++++++++++++++++++++++++++++++++++++++++++
C
        IF ((FSCT(1).EQ.0.0).AND.(IFLAG.NE.0)) THEN
           WRITE(*,*)
           WRITE(*,*) ' *** No Part_SSNR table - refinement ',
     .                           'will be less accurate'
           WRITE(*,*)
        ENDIF
      	IFLAG=IABS(IFLAG)
      	IF (NANG.LT.ILAST) ILAST=NANG
      	NSET=NSET-1
      	WRITE(*,*) 'End of input datasets, number of datasets =',NSET
      	IF (BOFF.EQ.0.0) THEN
          IF (NPHAS.NE.0) AMEANPHAS=AMEANPHAS/NPHAS
      	  BOFF=AMEANPHAS
      	  WRITE(*,6800)BOFF
6800	  FORMAT(' Average score used = ',F8.4)
      	ENDIF
        CALL FLUSH(6)
c----------------------------------------------------------------------------------
        IF (IFLAG.GE.3) THEN
          WRITE(*,*)
	  CALL SEARCHANG(RMM,ABS(RI),DANGIN,NSANG,IQUADMAX,ASYM,
     .	  NASYM,SANG,NDOC1,NSET,ASYMSTORE,PMASK,.FALSE.)
          ALLOCATE(SANG(3,NSANG),STAT=IERR)
          IF (IERR.NE.0)
     +      STOP ' ERROR: Memory allocation failed for SANG array'
	  CALL SEARCHANG(RMM,ABS(RI),DANGIN,NSANG,IQUADMAX,ASYM,
     .	  NASYM,SANG,NDOC1,NSET,ASYMSTORE,PMASK,.TRUE.)
          WRITE(*,*)
	ENDIF
C  end of loop for search angle generation.
c----------------------------------------------------------------------------------
      	RII=ABS(RI)+HALFW/2.0
      	RIH=ABS(RI)-HALFW/2.0
      	IF (RIH.LT.0.0) RIH=0.0
      	RI2=RIH**2
      	RI3=RII**2
      	RI4=RI**2
      	CALL DATE_AND_TIME(CDATE,CTIME,CZONE,DTVAL)
      	IRAN=-DTVAL(8)
        IF (ASYM(1:1).EQ.'H') THEN
441       CONTINUE
          IF (ALPHA.GT.180.0) THEN
            ALPHA=ALPHA-360.0
            GOTO 441
          ENDIF
          IF (ALPHA.LE.-180.0) THEN
            ALPHA=ALPHA+360.0
            GOTO 441
          ENDIF
          ALPHA=ALPHA/180.0*PI
          RISE=-RISE/PSIZE/NSAM*PI*2.0
        ENDIF

      	JC=NSAM/2+1
      	NSAMH=NSAM/2
      	STD=0.0
      	VTD=0.0

        IF (FALL) THEN
          DO 107 L=1,NVOL   
            NKSUM1(L)=0
            TTD1(L)=0.0D0
            IF (IFSC.LE.0) TTD2(L)=0.0D0
107       CONTINUE
      	  DO 103 K=1,NSAM
      	    DO 103 J=1,NSAM
      	      DO 103 I=1,JC
                ID=I+JC*((J-1)+NSAM*(K-1))
                DO 103 L=1,NVOL1
                  A3DF(ID,L)=CMPLX(0.0,0.0)
                  S3DF(ID,L)=0.0
103	  CONTINUE
          IF (IFSC.LE.0) THEN
      	    DO 102 K=1,NSAM
      	      DO 102 J=1,NSAM
      	        DO 102 I=1,JC
      	          ID=I+JC*((J-1)+NSAM*(K-1))
                  DO 102 L=1,NVOL2
      	            D3DF(ID,L)=CMPLX(0.0,0.0)
      	            V3DF(ID,L)=0.0
102	    CONTINUE
          ENDIF
          IF (NNSTAT.NE.0) THEN
      	    DO 101 K=1,NSAM
      	      DO 101 J=1,NSAM
      	        DO 101 I=1,JC
      	          ID=I+JC*((J-1)+NSAM*(K-1))
                  DO 104 L=1,NVOL
      	            KSUM(ID,L)=0
      	            PSUM(ID,L)=0.0
      	            VSUM(ID,L)=0.0
      	            ASUM(ID,L)=0.0
104               CONTINUE
      	          QCP(ID)=0.0
      	          SSNR(ID)=0.0
101	    CONTINUE
          ENDIF
          DO 80 K=1,NSAM/2
            VSN(K)=0.0D0
            VN(K)=0.0D0
            VVSN(K)=0.0D0
            VVN(K)=0.0D0
            IC(K)=0
80        CONTINUE
        ENDIF
c----------------------------------------------------------------------------------
	CALL CARDS13AND14(F3D,FWEIGH,FMATCH,FDEF,FMAG,INPIC,
     .              INSTAT,IFLAG,ISTAT,NSAM,NSAMH,NN1,JC,
     .              PSIZE,UTD,STDD,B3DF,C3DF,CFORM,ASYM,VX,
     .              ABS(RI),HALFW,XSTD,IPAD,IEWALD,FALL,FDUMP,
     .              FFTW_PLANS)
        CALL FLUSH(6)
c----------------------------------------------------------------------------------
      	NKSUM=0
      	NSTART=IFIRST
        ISWITCH=0
        IF (2*(NSTART/2).EQ.NSTART) ISWITCH=1
      	IF (IFSC.EQ.2) ISWITCH=1
        IF (IFSC.LT.0) IFSC=0
      	FFLAG=0
      	PSIS=0.0
      	THETAS=0.0
      	PHIS=0.0
      	DSHXS=0.0
      	DSHYS=0.0
      	DFMID1S=0.0
      	DFMID2S=0.0
      	ANGASTS=0.0
        OCCS=0.0
        ALGPS=0.0
        DALGPS=0.0
        DPRESS=0.0
      	ABSMAGS=0.0
      	PRESS=0.0
      	NSET=NSETST
        IEXCL=0
      	IBUF=0
      	IALL=0
      	IWALL=0
      	PALL=0.0
      	OALL=0.0
      	WALL=0.0
      	CCAVE=0.0
      	ICOUNT=0
      	ICCOUNT=0
      	DO 253 I=1,NSAM*2
      	  RBINS(I)=0.0
      	  RBINN(I)=0.0
      	  RAWBINS(I)=0.0
253	CONTINUE
c----------------------------------------------------------------------------------
	CALL CARDS15TO18(F3D,FWEIGH,F3D1,F3D2,FPHA,FPOI,NDOC1)
        CALL FLUSH(6)
c----------------------------------------------------------------------------------
C	MAIN LOOP OVER PARTICLES.....

      	DO 11 K=NSTART,ILAST
C
          DO 940 I=1,NSYM
            ISYM(I)=1
940       CONTINUE
          ISM=0
C
900       CONTINUE
          ISM=ISM+1
C
          DO 910 I=2,8
            SM(I)=0.0
910       CONTINUE
          SM(1)=1.0
          SM(5)=1.0
          SM(9)=1.0
          IF (NONSYM) THEN
            DO 920 I=1,NSYM
              IF (ISYM(I).NE.JSYM(I)) THEN
                DO 930 JJ=1,ISYM(I)
                  CALL MATMUL_T(SM,SYMOP(1,1,I),SM)
930             CONTINUE
              ENDIF
920         CONTINUE
          ENDIF
C
	  CALL LMAIN(K,FDEF,FMAG,INPROJ,CFORM,DATA,RREC,
     +      NSAM,DATA,SHX,SHY,CS,RAWBINS,WL,WGH1,WGH2,
     +      DFMID1,DFMID2,ANGAST,OCC,OUTD,OUTD,AMAGP,RIH,HALFW,
     +      HALFWC,RI2,RI3,RI4,DATD,DATD,B3DF,PHI,THETA,PSI,XSTD,
     +      IFLAG,IRADA,C3DF,RBINS,RBINN,NSET,PRESA,
     +      SANG,NSANG,CCD,CCD,IQUADMAX,NDOC1,NDOC2,PMASK,RBFAC,
     +      RMAX1,RMAX2,IRAN,PSIZE,DSTEP,FILM,NN1,NNPART,MAXSET,
     +      FHIST,ILIST,NPB,ITMAX,DAMAX,ASYMSTORE,TARGET,ICCOUNT,
     +      ISWITCH,IPAD,PALL,IALL,IRAD,ABS(RI),A3DF,S3DF,
     +      STD,D3DF,ILST,V3DF,VTD,PBC,BOFF,ASUM,
     +      VSUM,PSUM,KSUM,THRESH,NSAMH,JC,FFLAG,CCAVE,ICMP,IFIRST,
     +      IOPROJ,IEWALD,NSCF,SFILE,NANG,WGH,PSIS,THETAS,PHIS,
     +      DSHXS,DSHYS,PRESS,IEXCL,IBUF,ICOUNT,FMATCH,DFMID1S,
     +      DFMID2S,ANGASTS,OCCS,ALGPS,DPRESS,ABSMAGS,NSTART,
     +      FINPAT1,FINPAT2,ASYM,VX,ILAST,NNSET,IANG,CDATE,
     +      CTIME,CZONE,DTVAL,FASTIG,SINCLUT,MBUF,PBUF,CTFF,
     +      RBIN,CCBUF,QBUC,ISTAT,
     +      TNSET,FALL,IC,VSN,VN,VVSN,VVN,CI,BF,NS,NA,FSCT,
     +      TX,TY,TESTPAR,IPMAX,NSYM,ISYMAX,JSYM,SYMOP,NNSTAT,NKSUM,
     +      XM,YM,SX,SY,IFSC,IMP,NVOL,NVOL1,NVOL2,OALL,DFSTD,ALPHA,
     +      RISE,NU,HSTART,THETAM,STHETA,PSIM,SPSI,STIF,NPIX,NS1,
     +      ALGP,DALGPS,DMASK,NKSUM1,TTD1,TTD2,SIG,INTERP,PSSNR,
     +      PWEIGHTS,ISW,RECBF,IRECB,RECB,CTFBF,NN2,NN3,WALL,IWALL,
     +      RCLAS,DANGIN,FFTW_PLANS,NONSYM,ISM,SM,IREDUN)
C
          IF ((IREDUN.GT.1).AND.NONSYM) THEN
            I=1
990         CONTINUE
            ISYM(I)=ISYM(I)+1
            IF (ISYM(I).GT.JSYM(I)) THEN
              ISYM(I)=1
              I=I+1
              IF (I.LE.NSYM) GOTO 990
            ELSE
              I=1
              GOTO 900
            ENDIF
          ENDIF
C
11	CONTINUE

C END MAIN LOOP

	CALL ICLOSE(INPROJ)
        IF (FMATCH) CALL ICLOSE(IOPROJ)
C
        IF (FALL)
     +  CALL COMBINE_ARRAYS(NSAM,ASUM,VSUM,PSUM,KSUM,
     +           A3DF,S3DF,STD,D3DF,V3DF,
     +           VTD,NN1,NNSTAT,NKSUM,IFSC,IMP,NVOL,
     +           NVOL1,NVOL2,NKSUM1,TTD1,TTD2,NN2,NN3)
C
      	WRITE(*,6500)NSET,IEXCL
6500          FORMAT(/' PARTICLES BELOW THRESHOLD (NOT INCLUDED)',
     +               ' IN SET ',I3,': ',I10)
        CALL FLUSH(6)
C
        IF (OCCS.NE.0.0) THEN
          PRESS=PRESS/OCCS
          DPRESS=DPRESS/OCCS
          ALGPS=ALGPS/OCCS
          DALGPS=DALGPS/OCCS
        ENDIF
C
      	IF (ICOUNT.NE.0.0) THEN
          WRITE(NDOC2+NSET,7012)SQRT(PSIS/ICOUNT),SQRT(THETAS/ICOUNT),
     +          SQRT(PHIS/ICOUNT),SQRT(DSHXS/ICOUNT),SQRT(DSHYS/ICOUNT),
     +          NINT(SQRT(ABSMAGS/ICOUNT)),SQRT(DFMID1S/ICOUNT),
     +          SQRT(DFMID2S/ICOUNT),SQRT(ANGASTS/ICOUNT),0.0,
     +          NINT(DALGPS),DPRESS
            CALL FLUSH(NDOC2+NSET)
7012      FORMAT('C',6X,3F8.2,2F10.2,I8,6X,2F9.1,2F8.2,I10,
     +           11X,F8.2)

      	ENDIF
C
        IF (OALL.NE.0.0) PALL=PALL/OALL*100.0
        IF (IALL.NE.0) OALL=OALL/IALL*100.0
        IF (IWALL.NE.0) WALL=WALL/IWALL
        WRITE(*,16401)IALL,PALL,OALL
16401   FORMAT(I10,' PARTICLES FOR OVERALL SCORE CALCULATION',
     +        /6X,' SCORE (between resolution limits) =',F11.6,
     +        /6X,' AVERAGE OCCUPANCY                 =',F11.6)
        CALL FLUSH(6)
        IF (.NOT.FDUMP) WRITE(NDOC1+1,16402)IALL,PALL,OALL
16402   FORMAT('C  Total particles included, overall score, ',
     +         'average occupancy',I12,2F11.6)
C
      	IF (ICCOUNT.NE.0) THEN
      	  CCAVE=CCAVE/ICCOUNT
      	  WRITE (*,6403)ICCOUNT,CCAVE
6403	  FORMAT(/I6,' PARTICLES REFINED, AVERAGE CC',F12.8/)
          CALL FLUSH(6)
      	  WRITE(NDOC1+1,6404)ICCOUNT,CCAVE
6404	  FORMAT('C No. of particles refined, ave. CC',I6,F12.8)
      	ENDIF
C
      	IF(.NOT.FALL) THEN
          WRITE(*,*) ' Normal termination, no output data-files'
          GOTO 9999
        ENDIF
C
        IF (FDUMP) THEN
          CALL DUMP(STD,VTD,NN1,NN2,NN3,IREDUN,FCREF,RAWBINS,
     +          WC2,RLIM,RI,RIC,HALFW,ALPHA,RISE,HSTART,
     +          PSIZE,MW,DALT,CFORM,ASYM,VX,F3D,FWEIGH,
     +          F3D1,F3D2,FPHA,FPOI,FOUTPAR,QCP,ASUM,KSUM,
     +          PSUM,VSUM,A3DF,D3DF,S3DF,V3DF,
     +          IPAD,NNSTAT,NKSUM,NSYM,JSYM,SYMOP,
     +          INTERP,IFSC,INPIC,FBEAUT,FBFACT,IALL,NA,WALL,
     +          IWALL,VSN,VN,VVSN,VVN,NDOC1,NDOC2,NSET,IC,
     +          IEXCL,PALL,OALL,ICCOUNT,CCAVE)
          WRITE(*,*) ' Normal termination, intermediate files dumped'
          GOTO 9999
        ENDIF
C
C      	STD=STD/NSAM/NSAM/NSAM		! invalid for low resolution
      	STD=STD/NKSUM
      	STD=STD*WC2
C      	VTD=VTD/NSAM/NSAM/NSAM		! invalid for low resolution
      	VTD=VTD/NKSUM
      	VTD=VTD*WC2
C
C S3DF : stats file 1
C V3DF : stats file 2
C A3DF : reconstruction 1
C D3DF : reconstruction 2
C ASUM : sum of amplitudes
C VSUM : sum of amplitudes**2
C PSUM : sum of phase residuals
C KSUM : number of terms contributing to voxel
C
	CALL STORESHIFT(STD,VTD,JC,NSAM,NSAMH,IREDUN,NN1,
     +                  QCP,ASUM,KSUM,PSUM,VSUM,SSNR,A3DF,B3DF,
     +                  C3DF,D3DF,S3DF,V3DF,REF3DV,IPAD,NNSTAT,
     +                  IFSC)
C	TEST RESOLUTION - ISTP is step size in pixels for resolution statistics
      	ISTP=1
        IF (ASYM(1:1).EQ.'H') THEN
          I=NINT(2.0*PI/ALPHA)
          IZ1=NSAMH-NINT(ABS(I*RISE*NSAM/PI/4.0))
          IF (IZ1.LT.1) IZ1=1
          IZ2=NSAMH+NINT(ABS(I*RISE*NSAM/PI/4.0))-1
          IF (IZ2.GT.NSAM) IZ2=NSAM
        ENDIF
        IF (IFSC.EQ.0) THEN
!$OMP PARALLEL SECTIONS
!$OMP SECTION
          CALL FFTW_BWD(B3DF,B3DF,FFTW_PLANS(4))
          IF (ASYM(1:1).EQ.'H') THEN
            CALL CORRECT3D_C(NSAM,SINCLUT,B3DF,INTERP,0) 
          ELSE
            CALL CORRECT3D(NSAM,SINCLUT,B3DF,INTERP,0)
          ENDIF
          IF (ASYM(1:1).EQ.'H') THEN
            CALL MASK3D_C(NSAM,B3DF,ABS(RI),RIC,HALFW,
     +                    AFMASK,MAVE,MSTD)
          ELSE
            CALL MASK3D(NSAM,B3DF,ABS(RI),RIC,HALFW,
     +                  AFMASK,MAVE,MSTD)
          ENDIF
!$OMP SECTION
          CALL FFTW_BWD(C3DF,C3DF,FFTW_PLANS(4))
          IF (ASYM(1:1).EQ.'H') THEN
            CALL CORRECT3D_C(NSAM,SINCLUT,C3DF,INTERP,0)
          ELSE
            CALL CORRECT3D(NSAM,SINCLUT,C3DF,INTERP,0)
          ENDIF
          IF (ASYM(1:1).EQ.'H') THEN
            CALL MASK3D_C(NSAM,C3DF,ABS(RI),RIC,HALFW,
     +                    DUMMY,DUMMY,DUMMY)
          ELSE
            CALL MASK3D(NSAM,C3DF,ABS(RI),RIC,HALFW,
     +                    DUMMY,DUMMY,DUMMY)
          ENDIF
!$OMP END PARALLEL SECTIONS
          CALL OPMAPS2(F3D1,F3D2,I3D1,I3D2,
     .             JC,NSAM,NSAMH,PSIZE,B3DF,C3DF,
     .             CFORM,ASYM,VX,RBUF)
C
!$OMP PARALLEL SECTIONS
!$OMP SECTION
          CALL FFTW_FWD(B3DF,B3DF,FFTW_PLANS(3))
!$OMP SECTION
          CALL FFTW_FWD(C3DF,C3DF,FFTW_PLANS(3))
!$OMP END PARALLEL SECTIONS

C	PSUM contains average phase errors - not as informative as SSNR
C	-> PSUM replaced by SSNR. Change back if PSUM needed. Output of
C	shell-averaged values is still PF.
      	  CALL SHELTEST(NSAM,ISTP,B3DF,C3DF,ASUM,SSNR,KSUM,
     +		  QCP,NSTP,PR,FSC,ACC,QF,QC,PF,NS,NPIX,NPIXT,CBUF,
     +		  IREDUN,NNSTAT,PSSNR,S3DF,V3DF,AFMASK)
C
          CALL FIND_FPART(STD,VTD,JC,NSAM,NSAMH,
     +                    A3DF,B3DF,C3DF,D3DF,S3DF,V3DF,
     +                    IREDUN,ASYM,ABS(RI),RIC,B3DF,C3DF,
     +                    PSIZE,AFPART,FSC,MW,DALT,PSSNR,FFTW_PLANS)
C
          WRITE(*,9995) AFMASK,AFPART
9995      FORMAT(/' Fraction_mask, Fraction_particle =',2F10.4/)
          WRITE(*,9994) (PSIZE*NSAM)**3*AFPART*DALT/1000.0
9994      FORMAT(/' This correpsonds to a particle molecular mass of',
     +           F10.4,' kDa'/)
        ELSE
          WRITE(*,9996)IFSC
9996      FORMAT(/' *** IFSC set to',I2,' to save memory:',
     +            ' FSC statistics not done ***'/)
          CALL FLUSH(6)
        ENDIF

        IF (.NOT.FCREF) AFPART=0.0
	CALL SHIFTVOL(NSAM,STD,VTD,JC,RLIM,NSAMH,IREDUN,
     +        A3DF,D3DF,S3DF,V3DF,ASUM,QCP,NNSTAT,
     +        IFSC,PSSNR,AFPART,FSC,AFMASK)
C        IF ((FCREF).AND.(IFSC.EQ.0)) THEN
C !$OMP PARALLEL SECTIONS
C !$OMP SECTION
C          CALL APPLYCREF(NSAM,A3DF,FSC)
C !$OMP SECTION
C          CALL APPLYCREF(NSAM,D3DF,FSC) 
C !$OMP END PARALLEL SECTIONS
C        ENDIF
!$OMP PARALLEL SECTIONS
!$OMP SECTION
        CALL FFTW_BWD(A3DF,A3DF,FFTW_PLANS(4))
        IF (ASYM(1:1).EQ.'H') THEN
          CALL CORRECT3D_C(NSAM,SINCLUT,A3DF,INTERP,0)
        ELSE
          CALL CORRECT3D(NSAM,SINCLUT,A3DF,INTERP,0)
        ENDIF
!$OMP SECTION
        IF (IFSC.EQ.0) THEN
          CALL FFTW_BWD(D3DF,D3DF,FFTW_PLANS(4))
          IF (ASYM(1:1).EQ.'H') THEN
            CALL CORRECT3D_C(NSAM,SINCLUT,D3DF,INTERP,0)
          ELSE
            CALL CORRECT3D(NSAM,SINCLUT,D3DF,INTERP,0)
          ENDIF
        ENDIF
!$OMP END PARALLEL SECTIONS

        IF(FBEAUT) THEN
          IF (ASYM(1:1).EQ.'H') THEN
            CALL HEXTEND(NSAM,ALPHA,RISE*NSAM/PI/2.0,
     +                   A3DF,B3DF,ABS(RI),HSTART,IZ1,IZ2,BF)
          ELSE
            CALL BEAUTIFY(NSAM,NSYM,SYMOP,JSYM,A3DF,B3DF,
     +                           ABS(RI),HALFW,IMP)
          ENDIF
        ENDIF
        IF (RI.GT.0.0) THEN
          IF (ASYM(1:1).EQ.'H') THEN
            CALL MASK3D_C(NSAM,A3DF,RI,RIC,HALFW,AFMASK,
     +                    MAVE,MSTD)
          ELSE
            CALL MASK3D(NSAM,A3DF,RI,RIC,HALFW,AFMASK,
     +                  MAVE,MSTD)
          ENDIF
        ENDIF

        IF (FBFACT) CALL BFACTORSUB(NSAM,A3DF,PSIZE,FSC,
     +              IFSC,AFMASK,AFPART,BFACT,BMASK,NS,BF,FFTW_PLANS)

        IF (IFSC.EQ.0) THEN
	  CALL OPRESSTATHALF(NDOC1,VSN,VN,RAWBINS,PF,PR,ISTP,NSTP,NSAM,
     .                     RLIM,PSIZE,FSC,QF,QC,NPIX,NPIXT,CDATE,
     .                     CTIME,CZONE,DTVAL,VX,NN1,IC,VVSN,VVN,IALL,
     .                     NA,PSSNR,AFMASK,AFPART,BFACT,BMASK,
     .                     WALL)
        ENDIF
	CALL OPMAPS(FPOI,FPHA,INSTAT,IPSTAT,
     .              IPOINT,INPIC,PSUM,
     .              NN1,JC,NSAM,NSAMH,PSIZE,A3DF,D3DF,S3DF,
     .              CFORM,ASYM,VX,RBUF,NNSTAT,IFSC)
C
        DO 59 I=1,NSET
          CLOSE(NDOC1+NSET)
          CLOSE(NDOC2+NSET)
59      CONTINUE

        IF (NNSTAT.NE.0) THEN

C	NOW CALCULATE BRIEF STATISTICS BETWEEN STORED REFERENCE TRANSFORM 
C         REF3D AND THE COMPLETE MERGED NEW DATA A3DF VERSUS RESOLUTION
        CALL FFTW_FWD(A3DF,A3DF,FFTW_PLANS(3))
      	CALL SHELTEST(NSAM,ISTP,A3DF,REF3DV,ASUM,SSNR,KSUM,
     +		  QCP,NSTP,PR,FSC,ACC,QF,QC,PF,NS,NPIX,NPIXT,CBUF,
     +            IREDUN,NNSTAT,PSSNR,S3DF,V3DF,AFMASK)
	CALL OPRESSTATMAPS(NSTP,ISTP,NSAM,RLIM,PSIZE,PR,
     +			FSC,ACC,NPIX,NPIXT,NN1)

        ELSE
          WRITE(*,9997)
9997      FORMAT(/' *** FSTAT set to F to save memory:',
     +            ' some statistics not done ***'/)
          CALL FLUSH(6)
        ENDIF

      	CALL DATE_AND_TIME(CDATE,CTIME,CZONE,DTVAL)
      	WRITE(*,9998) CDATE(7:8),CDATE(5:6),CDATE(1:4),
     .		CTIME(1:2),CTIME(3:4),CTIME(5:8)
9998	FORMAT(' Final date/time ',A2,'-',A2,'-',A4,'/',A2,':',A2,':',A4)
        CALL FLUSH(6)

      	WRITE(*,*) ' Normal termination of frealign'
C
9999    CONTINUE
C
	END
