C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE MAGREFINE(AMAGP,DFMID1,DFMID2,ANGAST,IFIRSTF,
     +		LASTF,NP,NSAM,MAXR1,MAXR2,MASK,
     +		C3DF,IRAD,PBUF,
     +          SHX,SHY,IPAD,CS,WL,WGH1,WGH2,THETATR,CTFF,
     +          AMAG,RIH,HALFW,RI2,RI3,RI4,PHI,THETA,PSI,RMAX1,
     +		RMAX2,XSTD,MBUF,ILST,DATC,IBUF,B3DV,
     +		DATD,SINCLUT,IVFLAG,RBFACT,OUTD,OUTC,QBUC,
     +          RBUF,FSCT,FILM,ILAST,IEWALD,TX,TY,ASYM,PWEIGHTS,
     +          PSSNR,FFTW_PLANS,SM)
C**************************************************************************
C  Refines MAGNIFICATION using Powell minimiser on just 1 parameter
C  Calls VA04.
C  Used in LMAIN.
C**************************************************************************
C
      USE ISO_C_BINDING
C
      IMPLICIT NONE

      INTEGER NSAM,MAXR1,MAXR2,IRAD,IBUF
      INTEGER NCYCLS,IFIRSTF,LASTF,I,NP,MASK(5),ILST(*)
      INTEGER IEWALD,IPAD,ILAST,IVFLAG,FILM(*)
      PARAMETER (NCYCLS=50)
      REAL XSTD,SINCLUT(*),RBFACT,OUTD(*),TX,TY
      REAL SHX(*),SHY(*),PBUF(*),RI2,RI3,RI4,RBUF(*)
      REAL PWEIGHTS(*),PSSNR(*),SM(9)
      REAL XPAR(6),EPAR(6),ESCALE,PRESA,XP(6),MBUF(*)
      REAL AMAGP(*),DFMID1(*),DFMID2(*),ANGAST(*),B3DV(*)
      REAL CS,WL,WGH1,WGH2,THETATR,RIH,HALFW,RMAX1,DATD(*)
      REAL AMAG,PHI(*),THETA(*),PSI(*),RMAX2,FSCT(*)
      COMPLEX C3DF(*),CTFF(*)
      COMPLEX OUTC(*),DATC(*),QBUC(*)
      CHARACTER ASYM*3
      TYPE(C_PTR) FFTW_PLANS(*)
      DATA EPAR/0.01,0.0,0.0,0.0,0.0,0.0/
      DATA ESCALE/100.0/
C**************************************************************************
C above assumes relative magnification is near 1.0, and steps more than 0.1
C   are not allowed with final accuracy +/- 0.001 (i.e. 1 part in 1000)

      WRITE(*,*)'Entering MAGREFINE with range',IFIRSTF,LASTF
      XPAR(1)=AMAGP(IFIRSTF)
      XPAR(2)=DFMID1(IFIRSTF)
      XPAR(3)=DFMID2(IFIRSTF)
      XPAR(4)=ANGAST(IFIRSTF)
      CALL VA04A(XPAR,EPAR,1,PRESA,ESCALE,0,1,NCYCLS,XP,NP,
     .   NSAM,MAXR1,MAXR2,MASK,C3DF,
     .   IRAD,PBUF,SHX,SHY,IPAD,CS,WL,WGH1,
     .   WGH2,THETATR,CTFF,AMAG,RIH,HALFW,RI2,RI3,RI4,
     .	 PHI,THETA,PSI,RMAX1,RMAX2,XSTD,MBUF,ILST,DATC,
     .   IBUF,B3DV,DATD,SINCLUT,IVFLAG,RBFACT,OUTD,OUTC,QBUC,
     .   RBUF,FSCT,DFMID1,DFMID2,ANGAST,IEWALD,TX,TY,
     .   0.0,0.0,0.0,0.0,0.0,0.0,ASYM,0.0,1.0,0.0,1.0,0.0,
     .   PWEIGHTS,PSSNR,FFTW_PLANS,SM)
      DO 10 I=IFIRSTF,ILAST
        IF (FILM(IFIRSTF).EQ.FILM(I)) THEN
          AMAGP(I)=XPAR(1)
        ELSE
          GOTO 20
        ENDIF
10    CONTINUE
20    CONTINUE
      WRITE(*,*)'Exit from MAGREFINE, relative magnification',XPAR(1)

      RETURN
      END
