C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------------
	SUBROUTINE A3D3_S(ISW,NSAM,IRAD,RI,RECB,
     +           DATC,CTFF,
     +           ASUM,VSUM,PSUM,C3DF,KSUM,PBC,BOFF,
     +           A3DF,S3DF,STD,D3DF,V3DF,
     +           VTD,NN1,MAXSET,SINCLUT,IRADA,
     +           IPAD,THETATR,IEWALD,NSYM,ISYMAX,JSYM,
     +           SYMOP,NNSTAT,NKSUM,IFSC,IMP,PSIZE,
     +           ASYM,ALPHA,RISE,NU,HSTART,NKSUM1,TTD1,
     +           TTD2,NS,NS1,INTERP,I,NVOL,NVOL1,NVOL2,
     +           NN2,NN3,SM,NONSYM)
C -----------------------------------------------------------------
C         Put alternating particles into A3 or D3
C	  Used in A3D3
C	  Calls PINSERT
C -----------------------------------------------------------------
	IMPLICIT NONE
	INTEGER NN1,PBC,IRAD,NSAM,ISW(*),MAXSET
        INTEGER IPAD,IRADA,NSYM,ISYMAX,NU,NS1(*),NN2,NN3
        INTEGER IEWALD,IMP,JSYM(*),IFSC,NNSTAT,HSTART
        INTEGER NS(*),INTERP,I,IS,NVOL,NVOL1,NVOL2
	INTEGER KSUM(NN3*NNSTAT+1,NVOL)
        INTEGER*8 NKSUM,NKSUM1(*)
        REAL S3DF(NN3,NVOL1),V3DF(NN3,NVOL2)
	REAL SYMOP(3,3,*),PSIZE,ON,BOFF,RI,SM(9)
        REAL ASUM(NN3*NNSTAT+1,NVOL)
        REAL VSUM(NN3*NNSTAT+1,NVOL)
        REAL PSUM(NN3*NNSTAT+1,NVOL)
	REAL SINCLUT(*)
        REAL THETATR,ALPHA,RISE,RECB(10,NVOL)
	DOUBLE PRECISION STD,VTD,TTD1(*),TTD2(*)
	COMPLEX CTFF(NN1*NN1,NVOL)
	COMPLEX DATC(NN1/2*NN1,NVOL),C3DF(*)
        COMPLEX A3DF(NN3,NVOL1),D3DF(NN3,NVOL2)
        CHARACTER ASYM*3
        LOGICAL NONSYM
C -----------------------------------------------------------------
        ON=RECB(7,I)/RECB(8,I)
        IS=I
        IF (IFSC.EQ.0) IS=INT((I-1)/2)+1
C
        IF (ISW(I).EQ.0) THEN
C
          IF (ASYM(1:1).EQ.'H') THEN
C
          CALL PINSERT_C(NSAM,IRAD,RI,RECB(10,I),
     +    RECB(6,I),DATC(1,I),CTFF(1,I),
     +    A3DF(1,IS),S3DF(1,IS),
     +    A3DF(1,IS),S3DF(1,IS),
     +    STD,RECB(1,I),RECB(2,I),RECB(3,I),RECB(4,I),
     +    RECB(5,I),RECB(9,I),PBC,BOFF,ASUM(1,I),VSUM(1,I),
     +    PSUM(1,I),C3DF,KSUM(1,I),SINCLUT,IRADA,IPAD,
     +    THETATR/RECB(6,I),IEWALD,NSYM,ISYMAX,JSYM,
     +    SYMOP,NNSTAT,NKSUM,IMP,PSIZE,ALPHA,RISE,NU,HSTART,
     +    NKSUM1(I),TTD1(IS),NS,NS1,ON,INTERP,SM,NONSYM)
C
          ELSE
C
          CALL PINSERT(NSAM,IRAD,RI,RECB(10,I),
     +    RECB(6,I),DATC(1,I),CTFF(1,I),
     +    A3DF(1,IS),S3DF(1,IS),
     +    A3DF(1,IS),S3DF(1,IS),
     +    STD,RECB(1,I),RECB(2,I),RECB(3,I),RECB(4,I),
     +    RECB(5,I),RECB(9,I),PBC,BOFF,ASUM(1,I),VSUM(1,I),
     +    PSUM(1,I),C3DF,KSUM(1,I),SINCLUT,IRADA,IPAD,
     +    THETATR/RECB(6,I),IEWALD,NSYM,ISYMAX,JSYM,
     +    SYMOP,NNSTAT,NKSUM,IMP,PSIZE,
     +    NKSUM1(I),TTD1(IS),NS,NS1,ON,INTERP,SM,NONSYM)
C
          ENDIF
C
        ELSE
C
          IF (ASYM(1:1).EQ.'H') THEN
C
          CALL PINSERT_C(NSAM,IRAD,RI,RECB(10,I),
     +    RECB(6,I),DATC(1,I),CTFF(1,I),
     +    D3DF(1,IS),V3DF(1,IS),
     +    D3DF(1,IS),V3DF(1,IS),
     +    VTD,RECB(1,I),RECB(2,I),RECB(3,I),RECB(4,I),
     +    RECB(5,I),RECB(9,I),PBC,BOFF,ASUM(1,I),VSUM(1,I),
     +    PSUM(1,I),C3DF,KSUM(1,I),SINCLUT,IRADA,IPAD,
     +    THETATR/RECB(6,I),IEWALD,NSYM,ISYMAX,JSYM,
     +    SYMOP,NNSTAT,NKSUM,IMP,PSIZE,ALPHA,RISE,NU,HSTART,
     +    NKSUM1(I),TTD2(IS),NS,NS1,ON,INTERP,SM,NONSYM)
C
          ELSE
C
          CALL PINSERT(NSAM,IRAD,RI,RECB(10,I),
     +    RECB(6,I),DATC(1,I),CTFF(1,I),
     +    D3DF(1,IS),V3DF(1,IS),
     +    D3DF(1,IS),V3DF(1,IS),
     +    VTD,RECB(1,I),RECB(2,I),RECB(3,I),RECB(4,I),
     +    RECB(5,I),RECB(9,I),PBC,BOFF,ASUM(1,I),VSUM(1,I),
     +    PSUM(1,I),C3DF,KSUM(1,I),SINCLUT,IRADA,IPAD,
     +    THETATR/RECB(6,I),IEWALD,NSYM,ISYMAX,JSYM,
     +    SYMOP,NNSTAT,NKSUM,IMP,PSIZE,
     +    NKSUM1(I),TTD2(IS),NS,NS1,ON,INTERP,SM,NONSYM)
C
          ENDIF
        ENDIF
	RETURN
	END
