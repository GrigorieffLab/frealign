C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
c ---------------------------------------------------------------------
        SUBROUTINE OPRESSTATMAPS(NSTP,ISTP,NSAM,RLIM,PSIZE,PR,
     +                  FSC,ACC,NPIX,NPIXT,NN1)
c ---------------------------------------------------------------------
C	O/P RESOLUTION STATISTICS BETWEEN TWO MAPS
C	Used in Frealign.
c ---------------------------------------------------------------------
        IMPLICIT NONE
        REAL PI
        PARAMETER  (PI=3.1415926535897)
        INTEGER NSTP,ISTP,NSAM,NN1,I
	INTEGER NPIX(*),NPIXT(*)
        REAL RLIM,PSIZE,RINGRAD
	REAL PR(*),FSC(*),ACC(*)
c ---------------------------------------------------------------------
        WRITE(*,7017)
7017    FORMAT('C Statistics of new merged data versus reference data'/
     +         'C  NO. RESOL  RING RAD   FSPR   FSC ',
     +         '  RFACT NONZERO  TOTVOX')
        DO 63 I=2,NSTP
          RINGRAD=REAL(I-1)*ISTP/NSAM
          IF(RINGRAD.GT.RLIM) GO TO 63
          WRITE(*,7016)I,PSIZE/RINGRAD,RINGRAD,PR(I)*180.0/PI,
     +          FSC(I),ACC(I),NPIX(I),NPIXT(I)
7016      FORMAT('C',I4,F7.1,F10.4,F7.2,F7.3,F7.3,2I8)
63      CONTINUE
	RETURN
	END

