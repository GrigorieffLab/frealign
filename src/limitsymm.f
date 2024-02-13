C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE LIMITSYMM(ASYM,NASYM,THETAMAX,PHIMAX,JMAX,IQUADMAX)
C**************************************************************************
C  Restricts the number of search angles required for any requested point 
C    group symmetry.
C       ----------------- point group ----------------------
C              Cn        Dn         T         O         I
C
C phimax      360/N     360/N     180        90       180
C thetamax     90        90        45->54.7  45->54.7  20.8->31.7
C psimax       90        90        90        90        90
C jmax          2         1         1         1         1
C iquadmax      4         4         4         4         4

C  Used in SEARCHANG.
C**************************************************************************
      IMPLICIT NONE

      INTEGER J,JMAX,IQUADMAX,NASYM,N_CDTOI
      REAL THETAMAX,PHIMAX
      REAL THETASTORE(5),PHISTORE(5),JSTORE(5),IQUADSTORE(5)
      CHARACTER*1 ASYM
      CHARACTER*16 ASYMTEST
      DATA  ASYMTEST/' CDTOI0123456789'/
      DATA  THETASTORE/90.0,  90.0,  54.7,  54.7,  31.7/
      DATA  PHISTORE/360.0,  360.0, 180.0,  90.0, 180.0/
      DATA  JSTORE/2,1,1,1,1/
      DATA  IQUADSTORE/4,4,4,4,4/
C**************************************************************************
      	WRITE(*,*)' Entering LIMITSYMM with ASYM,NASYM  ',ASYM,NASYM
      		N_CDTOI = 0
      	DO 40 J = 2,6
        	IF(ASYM.EQ.ASYMTEST(J:J)) N_CDTOI = J-1
40    	CONTINUE
      IF(N_CDTOI.EQ.0) STOP ' Invalid call to LIMITSYMM'
      		THETAMAX=THETASTORE(N_CDTOI)
      		JMAX=JSTORE(N_CDTOI)
      		IQUADMAX=IQUADSTORE(N_CDTOI)
      	IF(N_CDTOI.LE.2) THEN
      		PHIMAX=PHISTORE(N_CDTOI)/NASYM
      	ELSE
      		PHIMAX=PHISTORE(N_CDTOI)
      	ENDIF
      RETURN
      END
