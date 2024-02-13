C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE D3MASK(NSAM,C3DV,B3DV,B3DF,XSTD,DMEAN,
     +                  FFTW_PLANS)
C**************************************************************************
C When XSTD is positive (see CALL), produces 5-pixel-cosine-bell smoothed 
C mask boundary defined by density being above the mean level of a low pass 
C filtered map by an amount equal to XSTD times the standard deviation of 
C the map.

C Calls FFTW.
C Used by MASK.
C**************************************************************************
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE

      INTEGER NSAM,L,LL,M,MM,N,NN,ID,IS,I,J,K,IRAD,JC,NSAM3
      INTEGER NM,NSAM32
      REAL RAD2,FRAD,B3DV(*),TLEVEL,XSTD,EDGE,PI,C3DV(*),SCAL
      PARAMETER (FRAD=0.1)
      PARAMETER (IRAD=5)
      PARAMETER (PI=3.1415926535897)
      DOUBLE PRECISION DMEAN,DSTD
      COMPLEX B3DF(*)
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
      JC=NSAM/2+1
C      NSAMH=NSAM/2
      NSAM3=NSAM*NSAM*NSAM
      NSAM32=NSAM*NSAM*(NSAM+2)
      SCAL=1.0/NSAM3

      DO 10 I=1,NSAM32
      	B3DV(I)=C3DV(I)
10    CONTINUE
      CALL FFTW_FWD(B3DV,B3DV,FFTW_PLANS(3))

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
      TLEVEL=DMEAN+ABS(XSTD)*DSTD
      NM=0
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
C     Smoothing the edge.....
      DO 11 I=1,NSAM
      	DO 11 J=1,NSAM
      	  DO 11 K=1,NSAM
      	    ID=K+NSAM*((J-1)+NSAM*(I-1))
      	    IF (B3DV(ID).EQ.1.0) THEN
      	      DO 12 L=-IRAD,IRAD
      	        DO 12 M=-IRAD,IRAD
      	          DO 12 N=-IRAD,IRAD
      	            RAD2=L**2+M**2+N**2
      	            RAD2=SQRT(RAD2)
                    EDGE=(1.0+COS(PI*RAD2/IRAD))/2.0
      		    ID=K+N+(NSAM+2)*((J+M-1)+NSAM*(I+L-1))
      	            IF ((ID.GE.1).AND.(ID.LE.NSAM32)
     +			.AND.(RAD2.LE.IRAD)) THEN
      	              IF (B3DV(ID).LT.EDGE) B3DV(ID)=EDGE
      	            ENDIF
12	      CONTINUE
      	    ENDIF
11    CONTINUE

C     Calculate mean within edge
      DMEAN=0.0D0
      NM=0
      DO 21 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 21 I=1,NSAM
          ID=IS+I
          IF ((B3DV(ID).LT.1.0).AND.(B3DV(ID).GT.0.0)) THEN
            DMEAN=DMEAN+C3DV(ID)
            NM=NM+1
          ENDIF
21    CONTINUE
      IF (NM.NE.0) THEN
        DMEAN=DMEAN/NM
      ENDIF

      DO 20 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 20 I=1,NSAM
          ID=IS+I
C          C3DV(ID)=C3DV(ID)*B3DV(ID)/NSAM3
      	  C3DV(ID)=C3DV(ID)*B3DV(ID)+(1.0-B3DV(ID))*DMEAN
20    CONTINUE

      RETURN
      END
