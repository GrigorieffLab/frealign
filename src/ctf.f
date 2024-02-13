C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      COMPLEX FUNCTION CTF(CS,WL,WGH1,WGH2IN,DFMID1,DFMID2,ANGAST,
     +                  THETATR,IX,IY,TX,TY)
C**************************************************************************
C  Calculates contrast transfer function, including the contribution from
C     amplitude contrast WGH2 - WGH1 is the resulting phase contrast
C  Used in CTFAPPLY and CTFAPPLY_PHASE_ONLY
C**************************************************************************
      IMPLICIT NONE
      INTEGER MAXSET,IX,IY
      REAL TWOPI
      PARAMETER (TWOPI=6.2831853071796)
      REAL CS,WL,WGH1,WGH2,DFMID1,DFMID2,WGH2IN
      REAL ANGAST,THETATR,RAD,ANGLE,ANGSPT,C1,C2,ANGDIF
      REAL CCOS,DF,CHI,SCHI,CCHI,TX,TY
C**************************************************************************
      WGH2=ABS(WGH2IN)
C
      IF (WGH2IN.NE.-1.0) THEN
C
      RAD=IX**2+IY**2
      IF (RAD.NE.0.0) THEN
        RAD=SQRT(RAD)
        ANGLE=RAD*THETATR
        ANGSPT=ATAN2(REAL(IY),REAL(IX))
        C1=TWOPI*ANGLE*ANGLE/(2.0*WL)
        C2=C1*CS*ANGLE*ANGLE/2.0
        ANGDIF=ANGSPT-ANGAST
        CCOS=COS(2.0*ANGDIF)
C       Positive defocus in FREALIGN is underfocus, convention used
C       in this formula is that underfocus should be negative.
C       => Take negative of defocus in following line
        DF=-0.5*(DFMID1+DFMID2+CCOS*(DFMID1-DFMID2))
        CHI=C1*DF+C2
C       Beam tilt
        CHI=CHI+TWOPI*CS*ANGLE*ANGLE*(IX*TX+IY*TY)*THETATR/WL/1000.0
        SCHI=SIN(CHI)
        CCHI=COS(CHI)
        CTF=CMPLX(WGH1*SCHI-WGH2*CCHI,-WGH1*CCHI-WGH2*SCHI)
      ELSE
C       CTF=CMPLX(-WGH2,-WGH1)
C       The imaginary part, -WGH1, will cancel when the left and right
C       beams are added later on => set to zero
        CTF=CMPLX(-WGH2,0.0)
      ENDIF
C
      IF (WGH2IN.LT.0.0) THEN
        CTF=-CTF
      ENDIF
C
      ELSE
C
      CTF=CMPLX(1.0,0.0)
C
      ENDIF

      RETURN
      END
