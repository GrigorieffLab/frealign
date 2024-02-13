C**************************************************************************
      SUBROUTINE CTFAPPLY_PHASE_ONLY(NSAM,SPEC,SHX,SHY,CS,WL,WGH1,
     +		 WGH2,DFMID1,DFMID2,ANGAST,THETATR,CTFF,OUTD,
     +		 OUTC,AMAGP,RIH,HALFW,RI2,RI3,RI4,
     +		 DATD,DATC,B3DV,PHI,THETA,PSI,W,XSTD,
     +           IBUF,TX,TY,ASYM,FFTW_PLANS)
C**************************************************************************
C  Applies calculated CTF to a particle transform and also masks the 
C  particle in real space with a cosine-bell-edged circle of radius RI
C  Calls FFTW and MASKENV.
C  Uses function CTF.
C  Used in CALCFX.
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE

      INTEGER NSAM,JC,L,LL,M,MM,ID,I,II,NMID,ISUM,IBUF
      INTEGER IC,JJ
      REAL PHASE,OUTD(*),FI1,AMAGP,CRAD2,SRAD2,RIH,W(*)
      REAL RI2,RI3,HALFW,PSHFTR,DATD(*),PI,TX,TY
      REAL RI4,SHX,SHY,CS,WL,WGH1,WGH2,DFMID1,DFMID2,ANGAST
      REAL THETATR,THETATRP
      REAL RIHP,SCAL
      REAL B3DV(*),PHI,THETA,PSI,XSTD
      PARAMETER  (PI=3.1415926535897)
      COMPLEX CTF,CTFV,CTFR,CTFL,CTFF(*)
      COMPLEX PSHFT,SPEC(*),OUTC(*),DATC(*)
      CHARACTER ASYM*3
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
C      write(*,*)'Now inside CTFAPPLY',AMAGP
      THETATRP=THETATR/AMAGP
      RIHP=RIH/AMAGP
      SCAL=1.0/NSAM/NSAM
C	  THETATRP is the diffraction angle for the first pixel in 
C	  the transform of a particle sampled at the input pixel size
      JC=NSAM/2+1
C      NSAMH=NSAM/2
      NMID=NSAM/2+1
      IC=NSAM*JC
C
C     CORRECT FOR CTF AND APPLY SHIFTS
      DO 111 L=1,JC
        LL=L-1
        DO 111 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          PHASE=SHX*LL+SHY*MM
          PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
          CTFR=CTF(CS,WL,WGH1,WGH2,DFMID1,DFMID2,
     +             ANGAST,THETATRP,LL,MM,TX,TY)
          CTFL=CTF(CS,WL,WGH1,WGH2,DFMID1,DFMID2,
     +             ANGAST,THETATRP,-LL,-MM,TX,TY)
C 	in this unweighted CTF-correction subroutine, the value of 
C 	CTFV can only be either +1.0 or -1.0
          CTFV=CTFR+CONJG(CTFL)
      	  IF(ABS(CTFV).NE.0.0) CTFV=CTFV/ABS(CTFV)
      	    ID=L+JC*(M-1)
C            OUTC(ID)=SPEC(ID)*CTFV*PSHFT
            OUTC(ID)=SPEC(ID)*PSHFT
            DATC(ID)=SPEC(ID)*PSHFT
      	    CTFF(ID)=CTFR
      	    CTFF(ID+IC)=CTFL
111   CONTINUE
C
C      write(*,*)'before FFTW'
      CALL FFTW_BWD(OUTD,OUTD,FFTW_PLANS(2))
      DO 21 I=0,NSAM-1
        ID=(NSAM+2)*I
        DO 21 II=1,NSAM
          JJ=II+ID
          OUTD(JJ)=OUTD(JJ)*SCAL
C        DATD(JJ)=OUTD(JJ)
21    CONTINUE

C     MASK PARTICLE, COSINE EDGE
      IF (ASYM(1:1).EQ.'H') THEN
        CALL MASKCOS_C(NSAM,OUTD,RI2,RI3,RIH,
     +                 HALFW,AMAGP,PSI)
      ELSE
        CALL MASKCOS(NSAM,OUTD,RI2,RI3,RIH,
     +               HALFW,AMAGP)
      ENDIF

      IF (XSTD.LT.0.0) THEN
C	Apply tight mask to image used for refinement.
      	CALL MASKENV(NSAM,RIHP,OUTD,B3DV,0.0,0.0,
     +               PHI,THETA,PSI,W,IBUF)
      ENDIF

      CALL FFTW_FWD(OUTD,OUTD,FFTW_PLANS(1))
C      write(*,*)'after FFTW'

C     SHIFT PARTICLE TO ORIGIN
      DO 110 L=1,JC
        LL=L-1
        DO 110 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
      	  ISUM=(LL+MM)
      	  PSHFTR=1.0
      	  IF (MOD(ISUM,2).NE.0) PSHFTR=-1.0
      	    ID=L+JC*(M-1)
            DATC(ID)=DATC(ID)*PSHFTR
            OUTC(ID)=OUTC(ID)*PSHFTR
110   CONTINUE
C
      RETURN
      END
