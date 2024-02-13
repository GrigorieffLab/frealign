C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PINSERT(NSAM,IRAD,RI,RREC,AMAG,SPEC,CTFF,
     +		 A3DF,S3DF,B3DF,T3DF,
     +		 TTD,PHI,THETA,PSI,SHX,SHY,PRESA,PBC,BOFF,
     +		 ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +           IPAD,THET,IEWALD,NSYM,ISYMAX,JSYM,SYMOP,
     +           NNSTAT,NKSUM,IMP,PSIZE,NKSUM1,TTD1,NS,NS1,OCC,
     +           INTERP,SM,NONSYM)
C**************************************************************************
C  Adds an image transform into the accumulating 3D transform and also 
C  accumulates the CTF**2 for use in later normalisation.
C  Uses PINSERT_S.
C  Used in A3D3.
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,JC,I,J,IRAD,IC,IEWALD,NS1(*),NS(*)
      INTEGER L,NSAMH,ID,IS,IRREC,IRREC2,NSYM,IMP
      INTEGER KSUM(*),IPAD,IRADA,ISYMAX,JSYM(*),NNSTAT
      INTEGER NL,NJ,INTERP,ISYM(ISYMAX)
      INTEGER*8 NKSUM,NKSUM1(*)
      REAL DM(9),SHX,SHY,RREC,TLIM,TM(9),SYMOP(3,3,*)
      REAL CTFV2,CTFV2S2,S3DF(*),CPHI,SPHI,RI,ST
      REAL CTHE,STHE,CPSI,SPSI,PHI,THETA,PSI,PSIZE
      REAL PHASE,PBC,BOFF,WGT,PI,TMP,AMAG,OCC
      REAL RII,PREF,T3DF(*),PRESA,SM(9)
      REAL VSUM(*),PSUM(*),ASUM(*),SINCLUT(*),THET
      DOUBLEPRECISION TTD,TTD1(*)
      PARAMETER (PI=3.1415926535897)
      PARAMETER (TLIM=2.0)
      COMPLEX SAMP,SPEC(*),BSAMP,CTFF(*)
      COMPLEX A3DF(*),B3DF(*),C3DF(*)
      LOGICAL NONSYM
C**************************************************************************
      IF (OCC.NE.0.0) THEN
C
      	PREF=PBC/100.0
        JC=NSAM/2+1
        IC=NSAM*JC
      	NSAMH=NSAM/2
      	RII=2.0*PI*RI/NSAM
      	IRREC=INT(REAL(NSAM)*RREC)
C      	IF (IRREC.GT.NSAMH) IRREC=NSAMH
      	IRREC2=IRREC**2
        CPHI=COS(PHI)
        SPHI=SIN(PHI)
        CTHE=COS(THETA)
        STHE=SIN(THETA)
        CPSI=COS(PSI)
        SPSI=SIN(PSI)
      	TMP=(PRESA/PREF-BOFF/PBC)/PSIZE**2
        IF (TMP.GT.TLIM) THEN
      	  TMP=TLIM
C      	  PRINT *,'TLIM IN ACTION!!!!!'
      	ENDIF
      	TMP=2.0*TMP/REAL(NSAMH)**2
C
        TM(1)=(CPHI*CTHE*CPSI-SPHI*SPSI)/ABS(AMAG)
        TM(2)=(SPHI*CTHE*CPSI+CPHI*SPSI)/ABS(AMAG)
        TM(3)=-STHE*CPSI/ABS(AMAG)
        TM(4)=(-CPHI*CTHE*SPSI-SPHI*CPSI)/ABS(AMAG)
        TM(5)=(-SPHI*CTHE*SPSI+CPHI*CPSI)/ABS(AMAG)
        TM(6)=STHE*SPSI/ABS(AMAG)
        TM(7)=STHE*CPHI/ABS(AMAG)
        TM(8)=STHE*SPHI/ABS(AMAG)
        TM(9)=CTHE/ABS(AMAG)
C
        DO 940 I=1,NSYM
          ISYM(I)=1
940     CONTINUE
C
900     CONTINUE
C
        DO 910 I=1,9
          DM(I)=TM(I)
910     CONTINUE
        IF (NONSYM) THEN
          CALL MATMUL_T(DM,SM,DM)
        ELSE
          DO 920 I=1,NSYM
            IF (ISYM(I).NE.JSYM(I)) THEN
              DO 930 J=1,ISYM(I)
                CALL MATMUL_T(DM,SYMOP(1,1,I),DM)
930           CONTINUE
            ENDIF
920       CONTINUE
        ENDIF
C
        ST=REAL(2*JC-1)/REAL(IMP)
        IF ((IEWALD.EQ.0).AND.(NNSTAT.EQ.0)) THEN
C
        IF (INTERP.EQ.0) THEN
          IF (IMP.GT.1) THEN
!$OMP PARALLEL DO
          DO 20 I=1,IMP
            CALL PINSERT_NEARESTB(NSAM,RII,IRREC2,SPEC,CTFF,
     +             OCC,B3DF,T3DF,TTD1(I),SHX,SHY,
     +             ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +             IPAD,THET,IEWALD,NNSTAT,NKSUM1(I),JC,DM,
     +             TMP,NSAMH,IC,0.0,0.0,ST,I,NS1)
20        CONTINUE
          ELSE
            CALL PINSERT_NEAREST(NSAM,RII,IRREC2,SPEC,CTFF,
     +             OCC,B3DF,T3DF,TTD1(1),SHX,SHY,
     +             ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +             IPAD,THET,IEWALD,NNSTAT,NKSUM1(1),JC,DM,
     +             TMP,NSAMH,IC,0.0,0.0,ST,1,NS1)
          ENDIF
        ELSE
          IF (IMP.GT.1) THEN
!$OMP PARALLEL DO
          DO 21 I=1,IMP
            CALL PINSERT_TRILINB(NSAM,RII,IRREC2,SPEC,CTFF,
     +             OCC,B3DF,T3DF,TTD1(I),SHX,SHY,
     +             ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +             IPAD,THET,IEWALD,NNSTAT,NKSUM1(I),JC,DM,
     +             TMP,NSAMH,IC,0.0,0.0,ST,I,NS1)
21        CONTINUE
          ELSE
            CALL PINSERT_TRILIN(NSAM,RII,IRREC2,SPEC,CTFF,
     +             OCC,B3DF,T3DF,TTD1(1),SHX,SHY,
     +             ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +             IPAD,THET,IEWALD,NNSTAT,NKSUM1(1),JC,DM,
     +             TMP,NSAMH,IC,0.0,0.0,ST,1,NS1)
          ENDIF
        ENDIF
C
        ELSE
C
        IF (INTERP.EQ.0) THEN
          IF (IMP.GT.1) THEN
!$OMP PARALLEL DO
          DO 22 I=1,IMP
            CALL PINSERT_NEARESTB_SE(NSAM,RII,
     +             IRREC2,SPEC,CTFF,
     +             OCC,B3DF,T3DF,TTD1(I),SHX,SHY,
     +             ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +             IPAD,THET,IEWALD,NNSTAT,NKSUM1(I),JC,DM,
     +             TMP,NSAMH,IC,0.0,0.0,ST,I,NS1)
22        CONTINUE
          ELSE
            CALL PINSERT_NEAREST_SE(NSAM,RII,
     +             IRREC2,SPEC,CTFF,
     +             OCC,B3DF,T3DF,TTD1(1),SHX,SHY,
     +             ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +             IPAD,THET,IEWALD,NNSTAT,NKSUM1(1),JC,DM,
     +             TMP,NSAMH,IC,0.0,0.0,ST,1,NS1)
          ENDIF
        ELSE
          IF (IMP.GT.1) THEN
!$OMP PARALLEL DO
          DO 23 I=1,IMP
            CALL PINSERT_TRILINB_SE(NSAM,RII,
     +             IRREC2,SPEC,CTFF,
     +             OCC,B3DF,T3DF,TTD1(I),SHX,SHY,
     +             ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +             IPAD,THET,IEWALD,NNSTAT,NKSUM1(I),JC,DM,
     +             TMP,NSAMH,IC,0.0,0.0,ST,I,NS1)
23        CONTINUE
          ELSE
            CALL PINSERT_TRILIN_SE(NSAM,RII,
     +             IRREC2,SPEC,CTFF,
     +             OCC,B3DF,T3DF,TTD1(1),SHX,SHY,
     +             ASUM,VSUM,PSUM,C3DF,KSUM,SINCLUT,IRADA,
     +             IPAD,THET,IEWALD,NNSTAT,NKSUM1(1),JC,DM,
     +             TMP,NSAMH,IC,0.0,0.0,ST,1,NS1)
          ENDIF
        ENDIF
C
        ENDIF
C
        IF (.NOT.NONSYM) THEN
          I=1
990       CONTINUE
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
      ENDIF
C
      RETURN
      END
