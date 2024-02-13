C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------------
	SUBROUTINE A3D3(NSET,ISWITCH,NSAM,IRAD,RI,RREC,
     +           K,ILAST,AMAGP,DATC,CTFF,
     +           ASUM,VSUM,PSUM,C3DF,KSUM,
     +           PHI,THETA,PSI,DSHX,DSHY,PRESA,PBC,BOFF,
     +           A3DF,S3DF,STD,D3DF,V3DF,
     +           VTD,NN1,MAXSET,
     +           SINCLUT,IRADA,IPAD,THETATR,IEWALD,
     +           NSYM,ISYMAX,JSYM,SYMOP,NNSTAT,NKSUM,
     +           IFSC,IMP,NVOL,NVOL1,NVOL2,OCC,SIG2N,PSIZE,
     +           ASYM,ALPHA,RISE,NU,HSTART,NKSUM1,TTD1,TTD2,
     +           NS,NS1,INTERP,RECBF,RECB,ISW,
     +           IBUF,CTFBF,NN2,NN3,SM,NONSYM,ISM)
C -----------------------------------------------------------------
C         Put alternating particles into A3 or D3
C	  Used in LMAIN
C	  Calls A3D3_S
C -----------------------------------------------------------------
	IMPLICIT NONE
	INTEGER NN1,PBC,IRAD,NSAM,ISWITCH,MAXSET,NSET
        INTEGER IPAD,IRADA,NSYM,ISYMAX,NU,NS1(*),IS,K
        INTEGER IEWALD,IMP,JSYM(*),IFSC,NNSTAT,HSTART
        INTEGER NS(*),INTERP,ISW(*),IBUF,JC,I
        INTEGER NVOL,NVOL1,NVOL2,NN2,NN3,ISM
	INTEGER KSUM(NN3*NNSTAT+1,NVOL),ILAST
        INTEGER*8 NKSUM,NKSUM1(*)
	REAL RREC(*),AMAGP,DSHX,DSHY,BOFF,RI
        REAL S3DF(NN3,NVOL1),V3DF(NN3,NVOL2)
	REAL SYMOP(3,3,*),PSIZE,OCC,SIG2N,SM(9)
        REAL ASUM(NN3*NNSTAT+1,NVOL)
        REAL VSUM(NN3*NNSTAT+1,NVOL)
        REAL PSUM(NN3*NNSTAT+1,NVOL)
	REAL PSI,THETA,PHI,PRESA,SINCLUT(*)
        REAL THETATR,ALPHA,RISE,RECB(10,NVOL)
	DOUBLE PRECISION STD,VTD,TTD1(*),TTD2(*)
	COMPLEX CTFF(*),DATC(*),C3DF(*)
        COMPLEX A3DF(NN3,NVOL1),D3DF(NN3,NVOL2)
        COMPLEX RECBF(NN1/2*NN1,NVOL)
        COMPLEX CTFBF(NN1*NN1,NVOL)
        CHARACTER ASYM*3
        LOGICAL NONSYM
C -----------------------------------------------------------------
        JC=NSAM/2+1
        IF (IBUF.LT.NVOL) THEN
C
        IF (ISWITCH.EQ.0) THEN
C
          IF (OCC.NE.0.0) THEN
            IBUF=IBUF+1
            ISW(IBUF)=ISWITCH
            DO 10 I=1,NSAM*JC
              RECBF(I,IBUF)=DATC(I)
10          CONTINUE
            DO 13 I=1,NSAM*(NSAM+2)
              CTFBF(I,IBUF)=CTFF(I)
13          CONTINUE
            RECB(1,IBUF)=PHI
            RECB(2,IBUF)=THETA
            RECB(3,IBUF)=PSI
            RECB(4,IBUF)=DSHX
            RECB(5,IBUF)=DSHY
            RECB(6,IBUF)=AMAGP
            RECB(7,IBUF)=OCC
            RECB(8,IBUF)=SIG2N
            RECB(9,IBUF)=PRESA
            RECB(10,IBUF)=RREC(NSET)
          ENDIF
C
          IF ((IFSC.LE.2).AND.(ISM.EQ.1)) ISWITCH=1
        ELSE
          IF (IFSC.EQ.0) THEN
C
            IF (OCC.NE.0.0) THEN
              IBUF=IBUF+1
              ISW(IBUF)=ISWITCH
              DO 20 I=1,NSAM*JC
                RECBF(I,IBUF)=DATC(I)
20            CONTINUE
              DO 23 I=1,NSAM*(NSAM+2)
                CTFBF(I,IBUF)=CTFF(I)
23            CONTINUE
              RECB(1,IBUF)=PHI
              RECB(2,IBUF)=THETA
              RECB(3,IBUF)=PSI
              RECB(4,IBUF)=DSHX
              RECB(5,IBUF)=DSHY
              RECB(6,IBUF)=AMAGP
              RECB(7,IBUF)=OCC
              RECB(8,IBUF)=SIG2N
              RECB(9,IBUF)=PRESA
              RECB(10,IBUF)=RREC(NSET)
            ENDIF
C
          ENDIF
          IF (ISM.EQ.1) ISWITCH=0
C
        ENDIF
C
        ENDIF
C
        IF ((IBUF.EQ.NVOL).OR.(K.EQ.ILAST)) THEN
C
        IF (IBUF.GT.1) THEN
C          DO 800 I=1,NVOL
C            NKSUM1(I)=0   
C            TTD1(I)=0.0D0
C            IF (IFSC.EQ.0) TTD2(I)=0.0D0
C800       CONTINUE
!$OMP PARALLEL DO
        DO 30 I=1,IBUF
          CALL A3D3_S(ISW,NSAM,IRAD,RI,RECB,
     +           RECBF,CTFBF,
     +           ASUM,VSUM,PSUM,C3DF,KSUM,PBC,BOFF,
     +           A3DF,S3DF,STD,D3DF,V3DF,
     +           VTD,NN1,MAXSET,SINCLUT,IRADA,
     +           IPAD,THETATR,IEWALD,NSYM,ISYMAX,JSYM,
     +           SYMOP,NNSTAT,NKSUM,IFSC,1,PSIZE,
     +           ASYM,ALPHA,RISE,NU,HSTART,NKSUM1,TTD1,
     +           TTD2,NS,NS1,INTERP,I,NVOL,NVOL1,NVOL2,
     +           NN2,NN3,SM,NONSYM)
30      CONTINUE
        ELSEIF (IBUF.EQ.1) THEN
          CALL A3D3_S(ISW,NSAM,IRAD,RI,RECB,
     +           RECBF,CTFBF,
     +           ASUM,VSUM,PSUM,C3DF,KSUM,PBC,BOFF,
     +           A3DF,S3DF,STD,D3DF,V3DF,
     +           VTD,NN1,MAXSET,SINCLUT,IRADA,
     +           IPAD,THETATR,IEWALD,NSYM,ISYMAX,JSYM,
     +           SYMOP,NNSTAT,NKSUM,IFSC,IMP,PSIZE,
     +           ASYM,ALPHA,RISE,NU,HSTART,NKSUM1,TTD1,
     +           TTD2,NS,NS1,INTERP,1,NVOL,NVOL1,NVOL2,
     +           NN2,NN3,SM,NONSYM)
        ENDIF
C
        IBUF=0
C
        ENDIF
C
	RETURN
	END
