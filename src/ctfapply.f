C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE CTFAPPLY(NSAM,SPEC,SHX,SHY,CS,WL,WGH1,
     +		 WGH2,DFMID1,DFMID2,ANGAST,THETATR,CTFF,OUTD,
     +		 OUTC,AMAGP,RIH,HALFW,RI2,RI3,RI4,
     +		 DATD,DATC,B3DV,PHI,THETA,PSI,W,XSTD,IBUF,
     +           TX,TY,ASYM,PSSNR,RBUF,PWEIGHTS,ICTF,FFTW_PLANS)
C**************************************************************************
C  Applies calculated CTF to a particle transform and also masks the 
C  particle in real space with a cosine-bell-edged circle of radius RI
C  Calls FFTW and MASKENV.
C  Uses function CTF.
C  Used in LMAIN.
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE

      INTEGER NSAM,JC,L,LL,M,MM,ID,I,II,JJ,NMID
      INTEGER ISUM,IBUF,IC,ICTF
      REAL PHASE,OUTD(*),FI1,AMAGP,CRAD2,SRAD2,RIH,W(*)
      REAL RI2,RI3,HALFW,PSHFTR,DATD(*),PI,TX,TY
      REAL RI4,SHX,SHY,CS,WL,WGH1,WGH2,DFMID1,DFMID2
      REAL THETATR,THETATRP,ANGAST,PX,PSSNR(*)
      REAL RIHP,SCAL,DF1OLD,DF2OLD,ANGOLD,PWEIGHTS(*)
      REAL B3DV(*),PHI,THETA,PSI,XSTD,RBUF(*)
      PARAMETER  (PI=3.1415926535897)
      COMPLEX CTF,CTFR,CTFL,CTFF(*)
      COMPLEX PSHFT,SPEC(*),OUTC(*),DATC(*)
      CHARACTER ASYM*3
      LOGICAL NEWCTF
      TYPE(C_PTR) FFTW_PLANS(*)
      SAVE DF1OLD,DF2OLD,ANGOLD
      DATA DF1OLD,DF2OLD,ANGOLD/0.0,0.0,0.0/
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
C     CORRECT FOR CTF AND APPLY SHIFTS
      NEWCTF=.FALSE.
      IF ((DFMID1.NE.DF1OLD).OR.
     +    (DFMID2.NE.DF2OLD).OR.
     +    (ANGAST.NE.ANGOLD)) NEWCTF=.TRUE.
C
      IF (NEWCTF) THEN
        DF1OLD=DFMID1
        DF2OLD=DFMID2
        ANGOLD=ANGAST
C
      DO 111 L=1,JC
        LL=L-1
        PX=SHX*LL
        DO 111 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          PHASE=PX+SHY*MM
          PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
          CTFR=CTF(CS,WL,WGH1,WGH2,DFMID1,DFMID2,
     +             ANGAST,THETATRP,LL,MM,TX,TY)
          CTFL=CTF(CS,WL,WGH1,WGH2,DFMID1,DFMID2,
     +             ANGAST,THETATRP,-LL,-MM,TX,TY)
      	    ID=L+JC*(M-1)
C            OUTC(ID)=SPEC(ID)*(CTFR+CONJG(CTFL))*PSHFT
            OUTC(ID)=SPEC(ID)*PSHFT
            DATC(ID)=SPEC(ID)*PSHFT
      	    CTFF(ID)=CTFR
      	    CTFF(ID+IC)=CTFL
111   CONTINUE
C
      ELSE
C
      DO 112 L=1,JC
        LL=L-1
        PX=SHX*LL
        DO 112 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          PHASE=PX+SHY*MM
          PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
      	    ID=L+JC*(M-1)
      	    CTFR=CTFF(ID)
      	    CTFL=CTFF(ID+IC)
C            OUTC(ID)=SPEC(ID)*(CTFR+CONJG(CTFL))*PSHFT
            OUTC(ID)=SPEC(ID)*PSHFT
            DATC(ID)=SPEC(ID)*PSHFT
112   CONTINUE
C
      ENDIF
C
c      CALL PWEIGHT(NSAM,OUTC,OUTC,
c     +            PSSNR,CTFF,PWEIGHTS,RBUF)
C
      CALL FFTW_BWD(OUTD,OUTD,FFTW_PLANS(2))
      DO 21 I=0,NSAM-1
        ID=(NSAM+2)*I
        DO 21 II=1,NSAM
          JJ=II+ID
          OUTD(JJ)=OUTD(JJ)*SCAL
21    CONTINUE
C
C     MASK PARTICLE, COSINE EDGE
      IF (ASYM(1:1).EQ.'H') THEN
        CALL MASKCOS_C(NSAM,OUTD,RI2,RI3,RIH,
     +                 HALFW,AMAGP,PSI)
      ELSE
        CALL MASKCOS(NSAM,OUTD,RI2,RI3,RIH,
     +               HALFW,AMAGP)
      ENDIF
C
      IF (XSTD.LT.0.0) THEN
C	Apply tight mask to image used for refinement.
      	CALL MASKENV(NSAM,RIHP,OUTD,B3DV,0.0,0.0,
     +               PHI,THETA,PSI,W,IBUF)
      ENDIF
C
C      CALL IOPEN("testimage.mrc",66,"M",NSAM,NSAM,1,'NEW',
C     +           "C1 ",2.44,"123456789012345")
C      L=0
C      DO 50 M=0,NSAM-1
C        L=L+1
C        ID=1+(NSAM+2)*M
C        CALL IWRITE(66,outd(ID),L)
C50    CONTINUE
C      CALL ICLOSE(66)
C
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
      CALL PWEIGHT(NSAM,OUTC,OUTC,
     +        PSSNR,CTFF,PWEIGHTS,RBUF,ICTF)
C
      RETURN
      END
