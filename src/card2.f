C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
	SUBROUTINE CARD2(RI,RIC,PSIZE,MW,WGH,XSTD,PBC,BOFF,
     +			 DANGIN,IPMAX,ITMAX,ITMAXA,IPMAXX) 
C**************************************************************************
C Card 2: global parameters
C         RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX
C Called by Frealign.
C**************************************************************************
	IMPLICIT NONE
        INTEGER ITMAX,ITMAXA,IPMAX,IPMAXX
        REAL RI,RIC,PSIZE,WGH,XSTD,PBC,BOFF,DANGIN,DUMMY,MW
        CHARACTER*80 INLINE
C**************************************************************************
      	WRITE(*,*)' RO,RI,PSIZE,MW,WGH,XSTD,PBC,BOFF,DANG,',
     +            'ITMAX,IPMAX?'
        READ(*,7007)INLINE
7007    FORMAT(A80)
        READ(INLINE,*,ERR=99)RI,RIC,PSIZE,MW,WGH,XSTD,PBC,BOFF,DANGIN,
     +                       ITMAX,IPMAX
        IF(PBC.NE.0.0)GOTO 98
99      CONTINUE
        WRITE(*,*) 'Card 2 error. Trying old Card 2 input...'
        READ(INLINE,*,ERR=97)RI,RIC,PSIZE,WGH,XSTD,PBC,BOFF,DANGIN,
     +                       ITMAX,IPMAX
        MW=0.0
        IF(PBC.NE.0.0)GOTO 98
97      CONTINUE
        READ(INLINE,*)RI,PSIZE,WGH,DUMMY,XSTD,PBC,BOFF,DANGIN,ITMAX
        RIC=0.0
        IPMAX=10
        MW=0.0
98      CONTINUE
      	WRITE(*,17002)RI,RIC,PSIZE,MW,WGH,XSTD,PBC,BOFF,DANGIN,
     +                ITMAX,IPMAX
        IF (RI.EQ.0.0) STOP ' ERROR: RO = 0 (Card 2)'
      		RI=RI/PSIZE		  ! convert to RI in pixels
      		RIC=RIC/PSIZE		  ! convert to RIC in pixels
      		IF(ITMAX.EQ.0) ITMAX=ITMAXA
      		IF(IPMAX.EQ.0) IPMAX=1
      		IF(IPMAX.GT.IPMAXX)
     +   STOP ' ERROR: IPMAX too large (Card 2)'
17002	FORMAT(2F8.1,F6.2,F12.3,2F7.3,3F6.1,2I5)
	RETURN
	END
