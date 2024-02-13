C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
c -----------------------------------------------------------------------------------
	SUBROUTINE SEARCHANG(RMM,RI,DANGIN,NSANG,IQUADMAX,
     .	          ASYM,NASYM,SANG,NDOC1,NSET,ASYMSTORE,MASK,STORE)
c -----------------------------------------------------------------------------------
c Calculate search angles from DANG
c Calls LIMITSYMM to limit the search angles
c Used in Frealign.
c -----------------------------------------------------------------------------------
#ifdef _NAG
        USE F90_UNIX
#endif
        IMPLICIT NONE
        REAL PI
	PARAMETER  (PI=3.1415926535897)
	CHARACTER*3 ASYM
	INTEGER NSANG,NDOC1,NASYM
	INTEGER MASK(5),IQUADMAX,NSET,J,JMAX
        REAL SANG(3,*),RI,RMM,THETA1,THETAMAX,THETAMAX1,PHIMAX,PSIMAX
        REAL DPHI,DPHI1,PSI1,PHI1,DANG,DANGIN,PTTT,A25,A23,A253,D,C
        REAL PCIRC1,TMP,THETAMIN
	LOGICAL ASYMSTORE,STORE
c -----------------------------------------------------------------------------------
      	  DANG=1.0/RMM/RI/PI*180.0
      	  IF(DANGIN.NE.0.0) DANG=DANGIN
      	  TMP=INT(90.0/DANG+0.5)
      	  TMP=AMAX1(0.25,TMP)
      	  DANG=90.0/TMP

C   This part of the loop has been changed to restrict the search angles to those 
C   required for the symmetry described by the ASYM input parameter.
C   This reduces the search time required if there is symmetry in the particle.
C   The original general values are set before the symmetry subroutine

      	NSANG=0
      	THETAMAX=90.0*REAL(MASK(2))
      	THETAMIN=0.0
        IF (ASYM(1:1).EQ.'H') THETAMAX=90.0+6.0*REAL(MASK(2))
        IF (ASYM(1:1).EQ.'H') THETAMIN=90.0-6.0*REAL(MASK(2))
C   Set THETAMAX1 to 90.0 to cover full angular space (previously was set to 0.0)
      	THETAMAX1=90.0
C      	PSIMAX=89.999		! set below
      	PHIMAX=359.999*REAL(MASK(3))
      	JMAX=2
        IF (ASYM(1:1).EQ.'H') JMAX=1
      	IQUADMAX=4
      	IF(ASYMSTORE) THEN
      	  CALL LIMITSYMM(ASYM(1:1),NASYM,THETAMAX,PHIMAX,JMAX,IQUADMAX)
      	ENDIF

      	DO 187 THETA1=THETAMIN,THETAMAX,DANG
      	    IF ((THETA1.EQ.0.0).OR.(THETA1.EQ.180.0)) THEN
      	      DPHI=360.0
      	    ELSE
              DPHI=DANG/SIN(THETA1*PI/180.)
              DPHI=360.0/MAX0(INT(359.99/DPHI)-1,1)
C		angular sampling was adapted from subroutine VOEA (Paul Penczek)
C		SIND is the intrinsic sine function with the argument in degrees
      	    ENDIF
      	  DO 186 PHI1=0.0,PHIMAX,DPHI
            IF ((ABS(THETA1-90.0).GT.1.0).OR.(PHI1.LT.179.0)) THEN
      		IF(ASYM(1:1).EQ.'T'.OR.ASYM(1:1).EQ.'O') THEN
                        PTTT=ABS(AMOD(PHI1,90.0))
                        PTTT=AMIN1(PTTT,90.-PTTT)
                        THETAMAX1=ATAN2(1.0,COS(PTTT*PI/180.))*180./PI    ! TETRAHEDRAL,OCTAHEDRAL
C			 thetamax varies from 45 to 54.7 depending on phi
      		ELSEIF(ASYM(1:1).EQ.'I') THEN
      			A25=ATAN2(2.0,1.0+SQRT(5.0))
      			A23=ASIN(TAN(30.0*PI/180.)*TAN(A25))
      			A253=ATAN2(TAN(30.0*PI/180.)*COS(A23),1.0)
      			IF(NASYM.EQ.2) THEN
      				PCIRC1 = ABS(AMOD(PHI1,180.0))
      				PCIRC1 = AMIN1(PCIRC1,180.0-PCIRC1)*PI/180.0
      				D = SIN(A25)/(COS(PCIRC1)+SIN(PCIRC1)/TAN(A253))
      				C = COS(A25)-D*SIN(PCIRC1)*TAN(A23)
      				THETAMAX1=ATAN2(D,C)*180.0/PI
C				  thetamax varies from 20.8 to 31.7 depending on phi
      			ELSE
      				PCIRC1 = ABS(AMOD(PHI1+90.0,180.0))
      				PCIRC1 = AMIN1(PCIRC1,180.0-PCIRC1)*PI/180.0
      				D = SIN(A25)/(COS(PCIRC1)+SIN(PCIRC1)/TAN(A253))
      				C = COS(A25)-D*SIN(PCIRC1)*TAN(A23)
      				THETAMAX1=ATAN2(D,C)*180.0/PI
C				  thetamax varies from 20.8 to 31.7 depending on phi
      			ENDIF
      		ENDIF
      		IF(THETA1.GT.THETAMAX1) GO TO 186
      		PSIMAX=90.0*REAL(MASK(1))
            PSI1=0.0
C      	    DO 185 PSI1=0.0,PSIMAX-DANG,DANG
C		the psi angle is here searched only from 0 to 90 because 
C		the other 3 quadrants are generated inside PSEARCH/CCP by
C		index permutation of the x,y coordinates of the projection
C		using the PSEARCH variable IQUAD=1,4 or IQUADMAX
C		here the index J=2 gives the mirror image
C
      	      DO 185 J=1,JMAX
      	        NSANG=NSANG+1
C      		   IF(NSANG.GT.IANG) THEN
C      			WRITE(*,*) 'Program dimension IANG too small',
C     .					NSANG,'>',IANG
C      			STOP
C      		   ENDIF
C
                IF (STORE) THEN
C
      		IF (J.EQ.1) THEN
      	          SANG(1,NSANG)=PHI1
      	          SANG(2,NSANG)=THETA1
      	          SANG(3,NSANG)=PSI1
      		ELSE
      	          SANG(1,NSANG)=PHI1+180.0
      	          SANG(2,NSANG)=180.0-THETA1
      	          SANG(3,NSANG)=PSI1+180.0
      		ENDIF
      		IF (SANG(1,NSANG).GT.180.0) SANG(1,NSANG)=SANG(1,NSANG)-360.0 
      		IF (SANG(2,NSANG).GT.180.0) SANG(2,NSANG)=SANG(2,NSANG)-360.0 
      		IF (SANG(3,NSANG).GT.180.0) SANG(3,NSANG)=SANG(3,NSANG)-360.0 
C		write(*,*)'search phi,th,psi',NSANG,(SANG(M,NSANG),M=1,3)
C		  convert search angles into radians
     	        SANG(1,NSANG)=SANG(1,NSANG)/180.0*PI
      	        SANG(2,NSANG)=SANG(2,NSANG)/180.0*PI
      	        SANG(3,NSANG)=SANG(3,NSANG)/180.0*PI
C
                ENDIF
C
185	    CONTINUE
            ENDIF
186	  CONTINUE
187     CONTINUE

C      	  CLOSE(77)
C
          IF (STORE) THEN
C
      	  WRITE(*,6010) DANG,NSANG
      	  WRITE(NDOC1+NSET,6010) DANG,NSANG
	  CALL FLUSH(NDOC1+NSET)
6010      FORMAT('C Search angle incr.',F8.2,
     +           ' No. of search angles:',I10)
C
          ENDIF
C
	RETURN
	END
