C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE CARDS8AND9(NSET,MAXSET,FINPAT1,FINPAT2,
     .			INPROJ,CFORM,ASYM,
     .		        PSIZE,VX,NSAM,ICMP,RLIM,RMM,
     .			RMAX1,RMAX2,RREC,RI,MSIZE,RCLAS)
C**************************************************************************
C Read input cards 8 and 9 
C  Card 8      FINPAT1 - PARTICLE IMAGE STACK FILENAME
C  Card 9      FINPAT2 - MATCHING PROJECTIONS STACK (for O/P if FMATCH=T)
C Calls IOPEN and ICLOSE.
C Called by Frealign.
C*************************************************************************
      IMPLICIT NONE
      INTEGER N2,NSTACK,MSIZE
      INTEGER NSET,INPROJ
      INTEGER NSAM,MAXSET,ICMP
      REAL RLIM,RMM,RTEMP,RI
      REAL PSIZE,PSIZE1
      REAL RREC(*),RMAX1(*),RMAX2(*),RCLAS(*)
      CHARACTER*3 ASYM
      CHARACTER*15 VX
      CHARACTER CFORM
      CHARACTER*200 FINPAT,FINPAT1(*),FINPAT2(*)
C**************************************************************************
        WRITE(*,*)' PARTICLE IMAGE STACK ?'
        READ(*,7006) FINPAT1(NSET)
        WRITE(*,17006) FINPAT1(NSET)
        IF (NSET.EQ.1) THEN
          CALL IOPEN(FINPAT1(NSET),INPROJ,CFORM,NSAM,N2,NSTACK,
     +                'OLD',ASYM,PSIZE1,VX)
          CALL ICLOSE(INPROJ)
          IF (2*(NSAM/2).NE.NSAM) THEN
            STOP 'ERROR: Particle dimensions must be even!!!'
          ENDIF
          IF (NSAM.NE.N2) THEN
            STOP 'ERROR: Particle X, Y dimensions must be the same!'
          ENDIF
        ENDIF
C       check and convert resolution from Angstroms to reciprocal fraction
        IF(RMAX1(NSET).LT.RMAX2(NSET)) THEN     ! swap
                        RTEMP=RMAX1(NSET)
                        RMAX1(NSET)=RMAX2(NSET)
                        RMAX2(NSET)=RTEMP
        ENDIF
        IF(RREC(NSET).LT.0.5.OR.RREC(NSET).GT.100.0)
     +    STOP 'ERROR: Unrealistic resolution for RREC (Card 7)'
        IF(RMAX1(NSET).LT.10.0.OR.RMAX1(NSET).GT.1000.0)
     +    STOP 'ERROR: Unrealistic resolution for RMIN (Card 7)'
        IF(RMAX2(NSET).LT.0.5.OR.RMAX2(NSET).GT.200.0)
     +    STOP 'ERROR: Unrealistic resolution for RMAX (Card 7)'
        IF(RCLAS(NSET).LT.0.5.OR.RCLAS(NSET).GT.200.0)
     +    STOP 'ERROR: Unrealistic resolution for RCLAS (Card 7)'
        RREC(NSET)=PSIZE/RREC(NSET)
        IF (RREC(NSET).GT.0.5) RREC(NSET)=0.5
C        RLIM=RREC(NSET)-2.0/NSAM
        RLIM=RREC(NSET)
        RMAX1(NSET)=PSIZE/RMAX1(NSET)
        RMAX2(NSET)=PSIZE/RMAX2(NSET)
        RCLAS(NSET)=PSIZE/RCLAS(NSET)
C        IF (RREC(NSET).GT.RLIM) RLIM=RREC(NSET)   ! highest resolution of datasets
        IF (RMAX2(NSET).GT.RLIM) RMAX2(NSET)=RLIM
        IF (RMAX2(NSET).GT.RMM) RMM=RMAX2(NSET)   ! only used if DANGIN not set
        IF (RCLAS(NSET).GT.RLIM) RCLAS(NSET)=RLIM
        ICMP=3.0*RI/MSIZE
        IF (3.0*RI.GT.NSAM) ICMP=NSAM/MSIZE
        IF (ICMP.EQ.0) ICMP=1
C
C        IF (FMATCH) THEN
C now always read in
          WRITE(*,*)' MATCHING PROJECTIONS IMAGE STACK ?'
          READ(*,7006) FINPAT2(NSET)
          WRITE(*,17006) FINPAT2(NSET)
C        ENDIF
7006            FORMAT(A200)
17006           FORMAT(3X,A200)
	RETURN
	END
