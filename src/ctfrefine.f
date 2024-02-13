C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE CTFREFINE(DFMID1,DFMID2,ANGAST,IFIRSTF,LASTF,FASTIG,
     +			NP,NSAM,MAXR1,MAXR2,MASK,DDFMID1,DDFMID2,
     +			DANGAST,C3DF,IRAD,PBUF,
     +			SHX,SHY,IPAD,CS,WL,WGH1,WGH2,THETATR,CTFF,
     +                  AMAG,RIH,HALFW,RI2,RI3,RI4,PHI,THETA,PSI,RMAX1,
     +		        RMAX2,XSTD,MBUF,ILST,DATC,IBUF,B3DV,
     +			DATD,SINCLUT,IVFLAG,RBFACT,OUTD,OUTC,QBUC,
     +                  CI,RBUF,FSCT,FILM,ILAST,IEWALD,TX,TY,SIG2,
     +                  DFSTD,ASYM,PWEIGHTS,PSSNR,FFTW_PLANS,SM)
C**************************************************************************
C  Refines defocus and astigmatism using Powell minimiser on 3 parameters
C  Calls VA04.
C  Used in LMAIN.
C**************************************************************************
C
      USE ISO_C_BINDING
C
      IMPLICIT NONE

      INTEGER NSAM,MAXR1,MAXR2,IRAD,IBUF,IVFLAG
      INTEGER NCYCLS,IFIRSTF,LASTF,I,NP,MASK(5),ILST(*)
      INTEGER IEWALD,IPAD,ILAST,FILM(*)
      PARAMETER (NCYCLS=50)
      REAL XSTD,SINCLUT(*),RBFACT,OUTD(*),RBUF(*),TX,TY
      REAL SHX,SHY,PBUF(*),RI2,RI3,RI4,DDFMID1,DDFMID2
      REAL XPAR(6),EPAR(6),ESCALE,PRESA,DFMID1(*),DFMID2(*)
      REAL ANGAST(*),XP(6),PHI,THETA,PSI,B3DV(*),FSCT(*)
      REAL CS,WL,WGH1,WGH2,THETATR,RIH,HALFW,MBUF(*)
      REAL AMAG(*),RMAX1,RMAX2,DATD(*),DANGAST,SIG2
      REAL DFSTD,PWEIGHTS(*),PSSNR(*),SM(9)
      COMPLEX C3DF(*),CTFF(*)
      COMPLEX OUTC(*),DATC(*),QBUC(*)
      LOGICAL FASTIG,CI
      CHARACTER ASYM*3
      TYPE(C_PTR) FFTW_PLANS(*)
      DATA EPAR/1.0,1.0,1.0,0.0,0.0,0.0/
      DATA ESCALE/1000.0/
C**************************************************************************
      IF (.NOT.CI) THEN
        WRITE(*,1000)IFIRSTF,LASTF,FASTIG,CI,DFSTD
1000    FORMAT('Entering CTFREFINE with range',2I10,2L4,F11.3)
        XPAR(1)=0.0
        XPAR(2)=0.0
        XPAR(3)=0.0
      ELSE
C        WRITE(*,1000)ILST(1),ILST(1),FASTIG
        XPAR(1)=0.0
        XPAR(2)=0.0
        XPAR(3)=0.0
      ENDIF
      IF(FASTIG) CALL VA04A(XPAR,EPAR,3,PRESA,ESCALE,0,1,NCYCLS,
     .   XP,NP,NSAM,MAXR1,MAXR2,MASK,C3DF,
     .	 IRAD,PBUF,SHX,SHY,IPAD,CS,WL,
     .	 WGH1,WGH2,THETATR,CTFF,AMAG,RIH,HALFW,RI2,RI3,RI4,
     .	 PHI,THETA,PSI,RMAX1,RMAX2,XSTD,MBUF,ILST,DATC,
     .   IBUF,B3DV,DATD,SINCLUT,IVFLAG,RBFACT,OUTD,OUTC,QBUC,
     .   RBUF,FSCT,DFMID1,DFMID2,ANGAST,IEWALD,TX,TY,0.0,0.0,
     .   0.0,0.0,SIG2,DFSTD,ASYM,0.0,1.0,0.0,1.0,0.0,PWEIGHTS,
     .   PSSNR,FFTW_PLANS,SM)
      IF(.NOT.FASTIG) CALL VA04A(XPAR,EPAR,1,PRESA,ESCALE,0,1,
     .   NCYCLS,XP,NP,NSAM,MAXR1,MAXR2,MASK,C3DF,
     .   IRAD,PBUF,SHX,SHY,IPAD,CS,WL,
     .   WGH1,WGH2,THETATR,CTFF,AMAG,RIH,HALFW,RI2,RI3,RI4,
     .   PHI,THETA,PSI,RMAX1,RMAX2,XSTD,MBUF,ILST,DATC,
     .   IBUF,B3DV,DATD,SINCLUT,IVFLAG,RBFACT,OUTD,OUTC,QBUC,
     .   RBUF,FSCT,DFMID1,DFMID2,ANGAST,IEWALD,TX,TY,0.0,0.0,
     .   0.0,0.0,SIG2,DFSTD,ASYM,0.0,1.0,0.0,1.0,0.0,PWEIGHTS,
     .   PSSNR,FFTW_PLANS,SM)
      IF(.NOT.FASTIG) XPAR(2)=XPAR(1)
      IF (.NOT.CI) THEN
        DO 10 I=IFIRSTF,ILAST
          IF (FILM(IFIRSTF).EQ.FILM(I)) THEN
            DFMID1(I)=DFMID1(I)+XPAR(1)
            DFMID2(I)=DFMID2(I)+XPAR(2)
            ANGAST(I)=ANGAST(I)+XPAR(3)
          ELSE
            GOTO 30
          ENDIF
10      CONTINUE
30      CONTINUE
        WRITE(*,1010)XPAR(1),XPAR(2),XPAR(3)*180.0/3.1415926535897
1010    FORMAT('Exit from CTFREFINE, defocus change',2F9.1,F8.2)
      ELSE
          DFMID1(ILST(1))=DFMID1(ILST(1))+XPAR(1)
          DFMID2(ILST(1))=DFMID2(ILST(1))+XPAR(2)
          ANGAST(ILST(1))=ANGAST(ILST(1))+XPAR(3)
      ENDIF
      DDFMID1=XPAR(1)
      DDFMID2=XPAR(2)
      DANGAST=XPAR(3)
C
      RETURN
      END
