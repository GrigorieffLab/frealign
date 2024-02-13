C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE CARD6(FALL,NSET,MAXSET,RELMAG,DSTEP,
     .			TARGET,THRESH,CS,AKV,WL,TX,TY)
C**************************************************************************
C Read input card 6 
C  RELMAG,DSTEP,TARGET,THRESH,CS,AKV
C              Relative magnification to apply         (1.0)
C              Densitometer step size in microns       (7.0)
C              Target score for search/refine          (15.0)
C              Worst score for inclusion               (90.0)
C              CS, KV                                  (2.0, 120)
C Called by Frealign.
C*************************************************************************
        IMPLICIT NONE
        INTEGER NSET,MAXSET
        REAL RELMAG(*),DSTEP(*),TARGET(*)
        REAL THRESH(*),CS(*),WL(*),TX(*),TY(*)
        REAL AKV
        LOGICAL FALL
        CHARACTER*80 INLINE
C*************************************************************************
	WRITE(*,*) ' MAGNIFICATION, STEPSIZE, TARGET & THRESH SCORE',
     .                  ' LIMIT, CS, KV, BEAM TILT X,Y ?'
        IF(NSET.GT.MAXSET) STOP '  too many datasets for MAXSET limit!'
        READ(*,7007)INLINE
7007    FORMAT(A80)
        READ(INLINE,*,ERR=99) RELMAG(NSET),DSTEP(NSET),TARGET(NSET),
     .            THRESH(NSET),CS(NSET),AKV,TX(NSET),TY(NSET)
        GOTO 98
99      CONTINUE
        WRITE(*,*) 'Card 6 error. Trying old CARD 6 input...'
        READ(INLINE,*) RELMAG(NSET),DSTEP(NSET),TARGET(NSET),
     .            THRESH(NSET),CS(NSET),AKV
        TX(NSET)=0.0
        TY(NSET)=0.0
98      CONTINUE
        WRITE(*,17004) RELMAG(NSET),DSTEP(NSET),TARGET(NSET),
     .          THRESH(NSET),CS(NSET),AKV,TX(NSET),TY(NSET)
17004   FORMAT(8F12.1)
        IF (RELMAG(NSET).LT.0.0) FALL=.FALSE.
c       IF (RELMAG(NSET).EQ.0.0.OR.(.NOT.FALL)) GOTO 1997
        IF (RELMAG(NSET).NE.0.0.AND.FALL) THEN
        WRITE(*,7030) NSET,RELMAG(NSET),DSTEP(NSET)
7030    FORMAT('  NSET=',I2,', RELMAG=',F7.3,', DSTEP=',F8.2)
        AKV=1000.0*AKV
        WL(NSET)=12.26/SQRT(AKV+0.9785*AKV**2/(10.0**6.0))
        ENDIF
1997	RETURN
	END
