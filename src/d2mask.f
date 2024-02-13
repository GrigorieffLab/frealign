C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE D2MASK(NSAM,B3DV,B3DF,XSTD,FFTW_PLANS)
C**************************************************************************
C When XSTD is negative (see CALL), produces a sharp (0/1) mask boundary 
C defined by density being above the mean level of a low pass filtered
C map by an amount equal to XSTD times the standard deviation of the map.

C Calls FFTW.
C Used by MASK.
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE

      INTEGER NSAM,L,LL,M,MM,N,NN,ID,IS,I,JC,NSAM3
      INTEGER NSAM32,J,NM
      REAL RAD2,FRAD,B3DV(*),TLEVEL,XSTD,SCAL
      PARAMETER (FRAD=0.1)
      DOUBLE PRECISION DMEAN,DSTD
      COMPLEX B3DF(*)
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
      JC=NSAM/2+1
C      NSAMH=NSAM/2
      NSAM3=NSAM*NSAM*NSAM
      NSAM32=NSAM*NSAM*(NSAM+2)
      SCAL=1.0/NSAM/NSAM/NSAM

      CALL FFTW_FWD(B3DV,B3DV,FFTW_PLANS(3))
C
C     Low-pass filter 3D map....
      DO 88 L=1,JC
        LL=L-1
        DO 88 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          DO 88 N=1,NSAM
            NN=N-1
            IF (NN.GE.JC) NN=NN-NSAM
            RAD2=REAL(LL**2+MM**2+NN**2)/NSAM**2/FRAD**2
              ID=L+JC*((M-1)+NSAM*(N-1))
              B3DF(ID)=B3DF(ID)*EXP(-RAD2)
88    CONTINUE

      CALL FFTW_BWD(B3DV,B3DV,FFTW_PLANS(4))
C
C     Calculate mean, STD
      DSTD=0.0D0
      DMEAN=0.0D0
      DO 87 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 87 I=1,NSAM
          ID=IS+I
          B3DV(ID)=B3DV(ID)*SCAL
          DMEAN=DMEAN+B3DV(ID)
          DSTD=DSTD+B3DV(ID)**2
87    CONTINUE
      DMEAN=DMEAN/NSAM3
      DSTD=DSTD/NSAM3
      DSTD=DSTD-DMEAN**2
      IF (DSTD.GT.0.0D0) THEN
        DSTD=SQRT(DSTD)
      ELSE
        DSTD=0.0D0
      ENDIF
C
      WRITE(*,1000) REAL(DMEAN),REAL(DSTD)
1000  FORMAT(' Mean, STD of low-pass filtered map: ',2F14.8)
C
C     Threshold filtered map using XSTD x STD above mean
      NM=0
      TLEVEL=DMEAN+ABS(XSTD)*DSTD
      DO 86 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 86 I=1,NSAM
          ID=IS+I
          IF (B3DV(ID).GE.TLEVEL) THEN
            B3DV(ID)=1.0
            NM=NM+1
          ELSE
            B3DV(ID)=0.0
          ENDIF
86    CONTINUE
C
      WRITE(*,1100) NM
1100  FORMAT(' Number of voxels inside mask:       ',I14,/)
C
      RETURN
      END
