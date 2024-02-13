C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      PROGRAM MERGE_3D
C**************************************************************************
C Program to merge data from parallelized 3D reconstructions
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE
C
      INTEGER IREDUN,ISTAT,INTERP,NN1,NN2,NN3,NA
      INTEGER IPAD,NNSTAT,IFSC,IERR,SLEN2,I,N,NSAM,ISTP
      INTEGER IZ1,IZ2,NSTP,HSTART,IMP,NSYM,ISYMAX,IALL
      PARAMETER (ISYMAX=100)
      INTEGER JSYM(ISYMAX),DTVAL(8),NDOC1,NSET
      INTEGER OMP_GET_NUM_PROCS,INSTAT,IPSTAT,IPOINT
      INTEGER INPIC,I3D1,I3D2,NSAMH,JC,IALLM,IREC
      INTEGER IWALL,IWALLM,IEXCL,IEXCLM,ICCOUNT,ICCOUNTM
      INTEGER*8 NKSUM,NKSUMM
      INTEGER,ALLOCATABLE :: KSUM(:),NPIX(:),NPIXT(:)
      INTEGER,ALLOCATABLE :: NS(:),IC(:),KSUMM(:)
      DOUBLEPRECISION STD,VTD,STDM,VTDM
      REAL WC2,RLIM,RI,RIC,HALFW,ALPHA,RISE,AFMASK,MW
      REAL SINCLUT(2000),MAVE,MSTD,DUMMY,PSIZE,DALT,PI
      REAL AFPART,SYMOP(3,3,ISYMAX),BFACT,BMASK,WALL
      REAL WALLM,PALL,OALL,CCAVE,PALLM,OALLM,CCAVEM
      REAL,ALLOCATABLE :: ASUM(:),VSUM(:),PSUM(:),BF(:)
      REAL,ALLOCATABLE :: S3DF(:),QCP(:),PSSNR(:)
      REAL,ALLOCATABLE :: V3DF(:),SSNR(:)
      REAL,ALLOCATABLE :: REF3DV(:),RBUF(:),PF(:),QC(:)
      REAL,ALLOCATABLE :: PR(:),FSC(:),ACC(:),QF(:)
      REAL,ALLOCATABLE :: RAWBINS(:),RAWBINSM(:)
      REAL,ALLOCATABLE :: ASUMM(:),VSUMM(:),PSUMM(:)
      REAL,ALLOCATABLE :: S3DFM(:),V3DFM(:)
      PARAMETER (PI=3.1415926535897)
      DOUBLEPRECISION,ALLOCATABLE :: VSN(:),VN(:)
      DOUBLEPRECISION,ALLOCATABLE :: VVSN(:),VVN(:)
      DOUBLEPRECISION,ALLOCATABLE :: VSNM(:),VNM(:)
      DOUBLEPRECISION,ALLOCATABLE :: VVSNM(:),VVNM(:)
      COMPLEX,ALLOCATABLE :: A3DF(:),B3DF(:)
      COMPLEX,ALLOCATABLE :: D3DF(:),C3DF(:)
      COMPLEX,ALLOCATABLE :: CBUF(:)
      COMPLEX,ALLOCATABLE :: A3DFM(:),D3DFM(:)
      CHARACTER*200 F3D,FWEIGH,F3D1,F3D2,FPHA,FPOI
      CHARACTER*200 F3DM,FWEIGHM,F3D1M,F3D2M,FPHAM,FPOIM
      CHARACTER*200 FOUTPAR,FOUTPARM
      CHARACTER VX*15,CFORM*1,ASYM*3,VXX*15
      CHARACTER CDATE*8,CTIME*10,CZONE*5,NCPUS*10
      CHARACTER*200,ALLOCATABLE :: FDUMP(:)
      LOGICAL FCREF,FBEAUT,FBFACT
      TYPE(C_PTR) FFTW_PLANS(10)
C**************************************************************************
C
      DATA VXX/'1.01 - 17.10.15'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C
      IMP=0
#ifdef _OPENMP
      CALL GETENV('OMP_NUM_THREADS',NCPUS)
      READ(NCPUS,*,ERR=111,END=111)IMP
111   CONTINUE
      IF (IMP.LE.0) THEN
        CALL GETENV('NCPUS',NCPUS)
        READ(NCPUS,*,ERR=112,END=112)IMP
112     CONTINUE
      ENDIF
      IF (IMP.LE.0) THEN
        IMP=OMP_GET_NUM_PROCS()
      ENDIF
#endif  
      IF (IMP.LE.0) IMP=1
C
      IF (IMP.GT.1) THEN
C
        WRITE(*,1011) VXX,IMP
1011    FORMAT(/'  MERGE_3D ',A15,/,
     +         /'  Copyright 2013 Howard Hughes Medical Institute.',
     +         /'  All rights reserved.',
     +         /'  Use is subject to Janelia Farm Research Campus',
     +          ' Software Copyright 1.1',
     +         /'  license terms ( http://license.janelia.org/license',
     +          '/jfrc_copyright_1_1.html )'/,
     *         /' Parallel processing: NCPUS =   ',I8,/)
C
      ELSE
C
        WRITE(*,1010)VXX
1010    FORMAT(/'  MERGE_3D ',A15,/,
     +         /'  Copyright 2013 Howard Hughes Medical Institute.',
     +         /'  All rights reserved.',
     +         /'  Use is subject to Janelia Farm Research Campus',
     +          ' Software Copyright 1.1',
     +         /'  license terms ( http://license.janelia.org/license',
     +          '/jfrc_copyright_1_1.html )'/)
C
      ENDIF
C
C  Read input files
C
      WRITE(*,2001)
2001  FORMAT(' No. of input files?')
      READ(*,*)N
      WRITE(*,3000)N
3000  FORMAT(I6)
      IF (N.EQ.0) THEN
        WRITE(*,*)
        WRITE(*,*) 'N = 0: Nothing to do. Exiting...'
        GOTO 9999
      ENDIF
      ALLOCATE(FDUMP(N),STAT=IERR)
      IF (IERR.NE.0) THEN
        STOP ' ERROR: Memory allocation failed'
      ENDIF
      WRITE(*,4010)
4010  FORMAT(/'*************************************************'/)
      DO 20 I=1,N
        WRITE(*,4000)I
4000    FORMAT(' Input file name for dump file ',I6)
        READ(*,7007)FDUMP(I)
7007    FORMAT(A200)
        WRITE(*,8000)FDUMP(I)(1:SLEN2(FDUMP(I)))
8000    FORMAT(A)
20    CONTINUE
      WRITE(*,4010)
      WRITE(*,*)
     +  ' Output res table (will be appended to existing file)?'
      READ(*,7007)FOUTPAR
      WRITE(*,8000)FOUTPAR(1:SLEN2(FOUTPAR))
      WRITE(*,*)' 3D map file for output?'
      READ(*,7007)F3D
      WRITE(*,8000)F3D(1:SLEN2(F3D))
      WRITE(*,*)' 3D weights file for output?'
      READ(*,7007)FWEIGH
      WRITE(*,8000)FWEIGH(1:SLEN2(FWEIGH))
      WRITE(*,*)' 3D reconstruction halfset 1 for output?'
      READ(*,7007)F3D1
      WRITE(*,8000)F3D1(1:SLEN2(F3D1))
      WRITE(*,*)' 3D reconstruction halfset 2 for output?'
      READ(*,7007)F3D2
      WRITE(*,8000)F3D2(1:SLEN2(F3D2))
      WRITE(*,*)' 3D phase residual file for output?'
      READ(*,7007)FPHA
      WRITE(*,8000)FPHA(1:SLEN2(FPHA))
      WRITE(*,*)' 3D point spread function for output?'
      READ(*,7007)FPOI
      WRITE(*,8000)FPOI(1:SLEN2(FPOI))
      WRITE(*,4010)
C
      CALL CALCSINC(SINCLUT,2000)
      BFACT=0.0
      BMASK=0.0
      INSTAT=11
      IPSTAT=12
      IPOINT=13
      INPIC=14
      I3D1=15
      I3D2=16
      NDOC1=16
C
      DO 10 I=1,N
C
        OPEN(UNIT=10,ACCESS='STREAM',FILE=FDUMP(I),
     .  ACTION='READ',FORM='UNFORMATTED')
C
        WRITE(*,*) 'Reading file',I
        IREC=1
        READ(10,POS=IREC,IOSTAT=ISTAT)NN1,NN2,NN3,NNSTAT,
     .    IPAD,IREDUN,NKSUMM,INTERP,IFSC
        IF (I.EQ.1) THEN
C
          NSAM=NN1
C
          ALLOCATE(KSUM(NN3*NNSTAT+1),ASUM(NN3*NNSTAT+1),
     .    VSUM(NN3*NNSTAT+1),PSUM(NN3*NNSTAT+1),S3DF(NN3),
     .    QCP(NN1/2*NN1*NN1*NNSTAT+NN1*NN1*NNSTAT+1),
     .    A3DF(NN3),
     .    SSNR(NN1/2*NN1*NN1*NNSTAT+NN1*NN1*NNSTAT+1),
     .    B3DF(NN1/2*NN1*NN1+NN1*NN1),
     .    C3DF(NN1/2*NN1*NN1+NN1*NN1),
     .    REF3DV(NN1*NN1*NN1*NNSTAT+2*NN1*NN1*NNSTAT+2),
     .    NPIX(NN1),IC(NN1),
     .    NPIXT(NN1/2+1),RBUF(NN1),NS(NN1),
     .    CBUF(NN1/2+1),PSSNR(NN1/2+1),BF(25*NN1*NN1+50*NN1),
     .    VSN(NN1),VN(NN1),VVSN(NN1),VVN(NN1),RAWBINSM(NN1*2),
     .    RAWBINS(NN1*2),VSNM(NN1),VNM(NN1),VVSNM(NN1),
     .    VVNM(NN1),KSUMM(NN3*NNSTAT+1),ASUMM(NN3*NNSTAT+1),
     .    VSUMM(NN3*NNSTAT+1),PSUMM(NN3*NNSTAT+1),S3DFM(NN3),
     .    A3DFM(NN3),STAT=IERR)
          IF (IERR.NE.0) THEN
            STOP ' ERROR: Memory allocation failed'
          ENDIF
C
          IF (IFSC.EQ.0) THEN
            ALLOCATE(V3DF(NN3),D3DF(NN3),
     +        V3DFM(NN3),D3DFM(NN3),
     +        STAT=IERR)
            IF (IERR.NE.0) THEN
              STOP ' ERROR: Memory allocation failed'
            ENDIF
          ENDIF
C
          IF ((IFSC.EQ.0).OR.(NNSTAT.NE.0)) THEN
            ALLOCATE(PR(NN1/2+1),FSC(NN1/2+1),ACC(NN1/2+1),
     +        QF(NN1/2+1),PF(NN1/2+1),QC(NN1/2+1),STAT=IERR)
            IF (IERR.NE.0) THEN
              STOP ' ERROR: Memory allocation failed'
            ENDIF
          ENDIF
C
          CALL FFTW_PLANS_2D(NN1,BF,BF,FFTW_PLANS(1),FFTW_PLANS(2))
          CALL FFTW_PLANS_3D(NN1,A3DF,A3DF,FFTW_PLANS(3),FFTW_PLANS(4))
C
          IEXCL=0
          PALL=0.0
          OALL=0.0
          ICCOUNT=0
          CCAVE=0.0
          IALL=0
          IWALL=0
          WALL=0.0
          NKSUM=0
          STD=0.0D0
          VTD=0.0D0
          VSN=0.0D0
          VN=0.0D0 
          VVSN=0.0D0
          VVN=0.0D0
          ASUM=0.0
          KSUM=0
          PSUM=0.0
          VSUM=0.0
          A3DF=CMPLX(0.0,0.0)
          S3DF=0.0
          RAWBINS=0.0
          IF (IFSC.EQ.0) THEN
            D3DF=CMPLX(0.0,0.0)
            V3DF=0.0
          ENDIF
C
#ifdef _OPENMP
          IF (NSAM/8.LT.IMP) THEN
            IMP=NSAM/8
            WRITE(*,6801) IMP
6801        FORMAT(/' Number of parallel processes too large.',
     .            ' Resetting IMP = ',I5/)
          ENDIF
          CALL OMP_SET_NUM_THREADS(IMP)
          FLUSH(6)
#endif
C
        ELSE
C
          IF(NN1.NE.NSAM) STOP ' ERROR: Incompatible dump files' 
C
        ENDIF
C
        IF (IFSC.EQ.0) THEN
          READ(10,POS=IREC,IOSTAT=ISTAT)NN1,NN2,NN3,NNSTAT,
     .      IPAD,IREDUN,NKSUMM,INTERP,IFSC,NSET,IWALLM,
     .      WC2,RLIM,RI,RIC,
     .      HALFW,ALPHA,RISE,HSTART,PSIZE,MW,DALT,NSYM,JSYM,
     .      SYMOP,IALLM,NA,WALLM,CFORM,ASYM,VX,FCREF,FBEAUT,
     .      FBFACT,F3DM,FWEIGHM,F3D1M,F3D2M,FPHAM,FPOIM,FOUTPARM,
     .      STDM,VTDM,VSNM,VNM,VVSNM,VVNM,ASUMM,KSUMM,PSUMM,VSUMM,
     .      QCP,A3DFM,S3DFM,D3DFM,V3DFM,
     .      RAWBINSM,IC,IEXCLM,PALLM,OALLM,ICCOUNTM,CCAVEM
        ELSE
          READ(10,POS=IREC,IOSTAT=ISTAT)NN1,NN2,NN3,NNSTAT,
     .      IPAD,IREDUN,NKSUMM,INTERP,IFSC,NSET,IWALLM,
     .      WC2,RLIM,RI,RIC,
     .      HALFW,ALPHA,RISE,HSTART,PSIZE,MW,DALT,NSYM,JSYM,
     .      SYMOP,IALLM,NA,WALLM,CFORM,ASYM,VX,FCREF,FBEAUT,
     .      FBFACT,F3DM,FWEIGHM,F3D1M,F3D2M,FPHAM,FPOIM,FOUTPARM,
     .      STDM,VTDM,VSNM,VNM,VVSNM,VVNM,ASUMM,KSUMM,PSUMM,VSUMM,
     .      QCP,A3DFM,S3DFM,RAWBINSM,IC,IEXCLM,PALLM,
     .      OALLM,ICCOUNTM,CCAVEM
        ENDIF
        IF (ISTAT.NE.0) THEN
          STOP ' ERROR: Dump file could not be read'
        ENDIF
C
        CLOSE(10)
C
        IEXCL=IEXCL+IEXCLM
        OALL=OALL+OALLM*IALLM/100.0
        OALLM=OALLM*IALLM/100.0
        PALL=PALL+PALLM*OALLM/100.0
        ICCOUNT=ICCOUNT+ICCOUNTM
        CCAVE=CCAVE+CCAVEM*ICCOUNTM
        IALL=IALL+IALLM
        WALL=WALL+WALLM*IWALLM
        IWALL=IWALL+IWALLM
        NKSUM=NKSUM+NKSUMM
        STD=STD+STDM
        VTD=VTD+VTDM
        VSN=VSN+VSNM
        VN=VN+VNM
        VVSN=VVSN+VVSNM
        VVN=VVN+VVNM
        ASUM=ASUM+ASUMM
        KSUM=KSUM+KSUMM
        PSUM=PSUM+PSUMM
        VSUM=VSUM+VSUMM
        A3DF=A3DF+A3DFM
        S3DF=S3DF+S3DFM
        RAWBINS=RAWBINS+RAWBINSM
        IF (IFSC.EQ.0) THEN
          D3DF=D3DF+D3DFM
          V3DF=V3DF+V3DFM
        ENDIF
C
10    CONTINUE
C
      I=NDOC1+1
      OPEN(UNIT=I,FILE=FOUTPAR,ACTION='WRITE',POSITION='APPEND')
C
      CALL IOPEN(F3D,INPIC,CFORM,NSAM,NSAM,NSAM,'NEW',ASYM,PSIZE,VX)
      CALL IOPEN(FWEIGH,INSTAT,CFORM,NSAM,NSAM,NSAM,'NEW',ASYM,
     .           PSIZE,VX)
C
      NSAMH=NSAM/2
      JC=NSAM/2+1
C
      WRITE(*,6500)NSET,IEXCL
6500  FORMAT(/' PARTICLES BELOW THRESHOLD (NOT INCLUDED)',
     +        ' IN SET ',I3,': ',I10)
      FLUSH(6)
C
      IF (OALL.NE.0.0) PALL=PALL/OALL*100.0
      IF (IALL.NE.0) OALL=OALL/IALL*100.0
      IF (IWALL.NE.0) WALL=WALL/IWALL
      WRITE(*,16401)IALL,PALL,OALL
16401 FORMAT(I10,' PARTICLES FOR OVERALL SCORE CALCULATION',
     +       /6X,' SCORE (between resolution limits) =',F11.6,
     +       /6X,' AVERAGE OCCUPANCY                 =',F11.6)
      FLUSH(6)
      WRITE(NDOC1+1,16402)IALL,PALL,OALL
16402 FORMAT('C  Total particles included, overall score, ',
     +       'average occupancy',I12,2F11.6)
C
      IF (ICCOUNT.NE.0) THEN
        CCAVE=CCAVE/ICCOUNT
        WRITE (*,6403)ICCOUNT,CCAVE
6403    FORMAT(/I6,' PARTICLES REFINED, AVERAGE CC',F12.8/)
        FLUSH(6)
        WRITE(NDOC1+1,6404)ICCOUNT,CCAVE
6404    FORMAT('C No. of particles refined, ave. CC',I6,F12.8)
      ENDIF
C
C     STD=STD/NSAM/NSAM/NSAM		! invalid for low resolution
      STD=STD/NKSUM
      STD=STD*WC2
C     VTD=VTD/NSAM/NSAM/NSAM		! invalid for low resolution
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
     +                QCP,ASUM,KSUM,PSUM,VSUM,SSNR,A3DF,B3DF,
     +                C3DF,D3DF,S3DF,V3DF,REF3DV,IPAD,NNSTAT,
     +                IFSC)
C     TEST RESOLUTION - ISTP is step size in pixels for resolution statistics
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
     +                  AFMASK,MAVE,MSTD)
        ELSE
          CALL MASK3D(NSAM,B3DF,ABS(RI),RIC,HALFW,
     +                AFMASK,MAVE,MSTD)
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
     +                  DUMMY,DUMMY,DUMMY)
        ELSE
          CALL MASK3D(NSAM,C3DF,ABS(RI),RIC,HALFW,
     +                  DUMMY,DUMMY,DUMMY)
        ENDIF
!$OMP END PARALLEL SECTIONS
        CALL OPMAPS2(F3D1,F3D2,I3D1,I3D2,
     .           JC,NSAM,NSAMH,PSIZE,B3DF,C3DF,
     .           CFORM,ASYM,VX,RBUF)
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
     +		QCP,NSTP,PR,FSC,ACC,QF,QC,PF,NS,NPIX,NPIXT,CBUF,
     +		IREDUN,NNSTAT,PSSNR,S3DF,V3DF,AFMASK)
C
        CALL FIND_FPART(STD,VTD,JC,NSAM,NSAMH,
     +                  A3DF,B3DF,C3DF,D3DF,S3DF,V3DF,
     +                  IREDUN,ASYM,ABS(RI),RIC,B3DF,C3DF,
     +                  PSIZE,AFPART,FSC,MW,DALT,PSSNR,FFTW_PLANS)
C
        WRITE(*,9995) AFMASK,AFPART
9995    FORMAT(/' Fraction_mask, Fraction_particle =',2F10.4/)
        WRITE(*,9994) (PSIZE*NSAM)**3*AFPART*DALT/1000.0
9994    FORMAT(/' This correpsonds to a particle molecular mass of',
     +           F10.4/)
      ELSE
        WRITE(*,9996)IFSC
9996    FORMAT(/' *** IFSC set to',I2,' to save memory:',
     +          ' FSC statistics not done ***'/)
        FLUSH(6)
      ENDIF
C
      IF (.NOT.FCREF) AFPART=0.0
      CALL SHIFTVOL(NSAM,STD,VTD,JC,RLIM,NSAMH,IREDUN,
     +     A3DF,D3DF,S3DF,V3DF,ASUM,QCP,NNSTAT,
     +     IFSC,PSSNR,AFPART,FSC,AFMASK)
C     IF ((FCREF).AND.(IFSC.EQ.0)) THEN
C !$OMP PARALLEL SECTIONS
C !$OMP SECTION
C       CALL APPLYCREF(NSAM,A3DF,FSC)
C !$OMP SECTION
C       CALL APPLYCREF(NSAM,D3DF,FSC) 
C !$OMP END PARALLEL SECTIONS
C     ENDIF
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
C
      IF(FBEAUT) THEN
        IF (ASYM(1:1).EQ.'H') THEN
          CALL HEXTEND(NSAM,ALPHA,RISE*NSAM/PI/2.0,
     +                 A3DF,B3DF,ABS(RI),HSTART,IZ1,IZ2,BF)
        ELSE
          CALL BEAUTIFY(NSAM,NSYM,SYMOP,JSYM,A3DF,B3DF,
     +                         ABS(RI),HALFW,IMP)
        ENDIF
      ENDIF
      IF (RI.GT.0.0) THEN
        IF (ASYM(1:1).EQ.'H') THEN
          CALL MASK3D_C(NSAM,A3DF,RI,RIC,HALFW,AFMASK,
     +                  MAVE,MSTD)
        ELSE
          CALL MASK3D(NSAM,A3DF,RI,RIC,HALFW,AFMASK,
     +                MAVE,MSTD)
        ENDIF
      ENDIF
C
      IF (FBFACT) CALL BFACTORSUB(NSAM,A3DF,PSIZE,FSC,
     +            IFSC,AFMASK,AFPART,BFACT,BMASK,NS,BF,FFTW_PLANS)
C
      IF (IFSC.EQ.0) THEN
	CALL OPRESSTATHALF(NDOC1,VSN,VN,RAWBINS,PF,PR,ISTP,NSTP,NSAM,
     .                   RLIM,PSIZE,FSC,QF,QC,NPIX,NPIXT,CDATE,
     .                   CTIME,CZONE,DTVAL,VX,NN1,IC,VVSN,VVN,IALL,
     .                   NA,PSSNR,AFMASK,AFPART,BFACT,BMASK,
     .                   WALL)
      ENDIF
      CALL OPMAPS(FPOI,FPHA,INSTAT,IPSTAT,
     .            IPOINT,INPIC,PSUM,
     .            NN1,JC,NSAM,NSAMH,PSIZE,A3DF,D3DF,S3DF,
     .            CFORM,ASYM,VX,RBUF,NNSTAT,IFSC)
      CLOSE(NDOC1+1)
C
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
9997    FORMAT(/' *** FSTAT set to F to save memory:',
     +          ' some statistics not done ***'/)
        FLUSH(6)
      ENDIF

9999  CONTINUE
C
      CALL DATE_AND_TIME(CDATE,CTIME,CZONE,DTVAL)
      WRITE(*,9998) CDATE(7:8),CDATE(5:6),CDATE(1:4),
     .	       CTIME(1:2),CTIME(3:4),CTIME(5:8)
9998  FORMAT(' Final date/time ',A2,'-',A2,'-',A4,'/',A2,':',A2,':',A4)
      FLUSH(6)
C
      WRITE(*,*)
      WRITE(*,*) ' Normal termination of merge_3d'
C
      END
