C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------------------------
      SUBROUTINE CARDS13AND14(F3D,FWEIGH,FMATCH,FDEF,FMAG,
     .			INPIC,INSTAT,IFLAG,ISTAT,
     .			NSAM,NSAMH,NN1,JC,PSIZE,
     .			UTD,STDD,
     .			B3DV,C3DV,
     .			CFORM,ASYM,VX,RI,HALFW,XSTD,
     .			IPAD,IEWALD,FALL,FDUMP,FFTW_PLANS)
C -----------------------------------------------------------------------------
C Cards 13 and 14 are filenames for the 3D structure I/O (overwrites input)
C Card 13         F3D    - 3D MAP FILE FOR INPUT and OUTPUT, Input 3D reference
C                        reconstruction, overwritten by output reconstruction!
C Card 14         FWEIGH - 3D WEIGHTS FILE FOR INPUT and OUTPUT, Input 3D file
C                        containing the sum of weights (as defined in the JMB
C                        reference), will be overwritten with new file!
C Calls IOPEN, IREAD, ICLOSE, MASKING, SHIFT.
C Called by Frealign.
C
C Output 3D arrays: B3D = binary mask for masking input particles (if XSTD < 0)
C                   C3D = reference volume for refinement
C
C -----------------------------------------------------------------------------
C
          USE ISO_C_BINDING
          USE FFTW33
C
          IMPLICIT NONE

          INTEGER NN1,INPIC,INSTAT,IFLAG,NSAM,NSAMH,JC,ISUM
          INTEGER N1,N2,N3,I,J,K,ID,NUTD,IPAD,I1,I2,ID4
          INTEGER ISTAT,IEWALD,SLEN
          REAL B3DV(*),C3DV(*),CIRC,STDC,PSIZE1
          REAL PSIZE,RI,XSTD,UTD,STDD,HALFW
          DOUBLE PRECISION STD,MEAN
          CHARACTER*200 F3D,F3DT,FWEIGH
          CHARACTER*1 CFORM
          CHARACTER*3 ASYM
          CHARACTER*15 VX
          LOGICAL FMATCH,FDEF,FMAG,EX,FALL,FDUMP
          TYPE(C_PTR) FFTW_PLANS(*)
c----------------------------------------------------------------------------------
        ISTAT=1

        IF (FALL) THEN
      	  WRITE(*,*)' 3D MAP FILE FOR INPUT (OVERWRITTEN BY OUTPUT)?'
        ELSE
      	  WRITE(*,*)' 3D MAP FILE FOR INPUT?'
        ENDIF
      	READ(*,7006)F3D
7006	FORMAT(A200)
      	WRITE(*,17006)F3D
17006	FORMAT(3X,A200)
C
        INQUIRE(FILE=F3D,EXIST=EX)
        IF ((CFORM.EQ.'I').AND.(.NOT.EX)) THEN
          F3DT=F3D(1:SLEN(F3D))//'.hed'
          INQUIRE(FILE=F3DT,EXIST=EX)
        ENDIF
C
      	N1=NSAM
      	N2=NSAM
      	N3=NSAM
        IF (.NOT.EX) THEN

      	  IF ((IFLAG.EQ.0).AND.(XSTD.EQ.0.0).AND.(.NOT.FMATCH).AND.
     +	      (.NOT.FDEF).AND.(.NOT.FMAG)) THEN
            WRITE(*,*)' No 3D map on input -- cannot do some',
     +                ' statistics at end of run'
            ISTAT=0
            IF (ABS(IEWALD).GT.1) THEN
              WRITE(*,*)' Cannot do model-based Ewald sphere',
     +                  ' correction, IEWALD set to 1'
              IEWALD=1
            ENDIF
      	    CALL IOPEN(F3D,INPIC,CFORM,N1,N2,N3,'NEW',ASYM,
     +		       PSIZE,VX)

          ELSE

            STOP ' 3D input map required!!!'

          ENDIF

        ELSE

      	  CALL IOPEN(F3D,INPIC,CFORM,N1,N2,N3,'OLD',ASYM,
     +		   PSIZE1,VX)
C          IF (ABS(PSIZE-PSIZE1)/PSIZE.GT.0.01)
C     +      WRITE(*,*)
C     +      ' ***WARNING*** Input and file pixel sizes differ!'

          IF (NSAM.NE.N1) THEN
          STOP ' Particle dimensions must match those of 3D reference!'
          ENDIF

      	  IF ((RI+HALFW/2.0).GT.REAL(NSAM/2)) THEN
      	    WRITE(*,*)
     +      ' ***WARNING*** Particle mask exceeds window size!'
      	  ENDIF

          K=0
          ISUM=0
          CIRC=0.0
          MEAN=0.0D0
C
          DO 701 I=1,NSAM
            DO 701 J=1,NSAM
              K=K+1
      	      ID=1+(NSAM+2)*((J-1)+NSAM*(I-1))
      	      CALL IREAD(INPIC,C3DV(ID),K)
              DO 711 I2=1,NSAM
      	        ID=I2+(NSAM+2)*((J-1)+NSAM*(I-1))
                MEAN=MEAN+C3DV(ID)
                IF ((I.EQ.1).OR.(J.EQ.1).OR.(I.EQ.NSAM).OR.(J.EQ.NSAM))
     +          THEN
                  CIRC=CIRC+C3DV(ID)
                  ISUM=ISUM+1
                ELSEIF ((I2.EQ.1).OR.(I2.EQ.NSAM)) THEN
                  CIRC=CIRC+C3DV(ID)
                  ISUM=ISUM+1
                ENDIF
711           CONTINUE
701       CONTINUE
C
          CALL ICLOSE(INPIC)
          CIRC=CIRC/ISUM
          MEAN=MEAN/NSAM/NSAM/NSAM
C
          STD=0.0D0
          STDC=0.0
C
          DO 741 I=1,NSAM
            DO 741 J=1,NSAM
              DO 741 I2=1,NSAM
      	        ID=I2+(NSAM+2)*((J-1)+NSAM*(I-1))
                STD=STD+(MEAN-C3DV(ID))**2
                IF ((I.EQ.1).OR.(J.EQ.1).OR.(I.EQ.NSAM).OR.(J.EQ.NSAM))
     +          THEN
                  STDC=STDC+(CIRC-C3DV(ID))**2
                ELSEIF ((I2.EQ.1).OR.(I2.EQ.NSAM)) THEN
                  STDC=STDC+(CIRC-C3DV(ID))**2
                ENDIF
741       CONTINUE
C
          STD=SQRT(STD/NSAM/NSAM/NSAM)
          STDC=SQRT(STDC/ISUM)
          IF (10.0*STDC.GT.STD) WRITE(*,*)
     +        ' ***WARNING*** Circumference STD',
     +        ' significant compared with volume STD, could indicate',
     +        ' noisy reference or density gradient !!!'
          WRITE(*,1000) REAL(MEAN),REAL(STD),CIRC,STDC
1000      FORMAT(/,' Mean density of 3D volume:     ',F14.8,
     +           /,' STD of 3D volume:              ',F14.8,
     +           /,' Mean density of circumference: ',F14.8,
     +           /,' STD of circumference:          ',F14.8,/)
C
C         Mask 3D file.....
C
          IF (XSTD.GT.0.0) THEN
            WRITE(*,6301)
6301        FORMAT(' Masking 3D model...')
            CALL D3MASK(NSAM,C3DV,B3DV,B3DV,XSTD,MEAN,FFTW_PLANS)
            CIRC=MEAN
          ENDIF
C
          IF ((XSTD.LT.0.0).OR.(IPAD.NE.1)) THEN
            DO 702 I=1,NSAM*NSAM*(NSAM+2)
              B3DV(I)=C3DV(I)
702         CONTINUE
          ENDIF
C
C         Pad volume.....
C
          IF (IPAD.NE.1) THEN
            DO 72 I=1,IPAD*IPAD*(IPAD*NSAM+2)*NSAM*NSAM
              C3DV(I)=CIRC
72          CONTINUE
            I1=INT(REAL(IPAD-1)/2.0*NSAM)
            DO 73 I=1,NSAM
              DO 73 J=1,NSAM
                DO 73 I2=1,NSAM
      	          ID=I2+(NSAM+2)*((J-1)+NSAM*(I-1))
                  ID4=I1+I2+(IPAD*NSAM+2)*((I1+J-1)+IPAD*NSAM*(I1+I-1))
                  C3DV(ID4)=B3DV(ID)
73          CONTINUE
          ENDIF

      	  N1=NSAM
      	  N2=NSAM
      	  N3=NSAM
          IF (FALL) THEN
      	    CALL IOPEN(F3D,INPIC,CFORM,N1,N2,N3,'NEW',ASYM,PSIZE,
     +			VX)
          ENDIF
C
C         Fourier-transform reference volume....
C
          CALL FFTW_FWD(C3DV,C3DV,FFTW_PLANS(5))
C
C         Generate binary 3D mask file for 2D masks.....
C
          IF (XSTD.LT.0.0) THEN
            WRITE(*,6302)
6302        FORMAT(' Using 3D model to generate 2D mask...')
            CALL D2MASK(NSAM,B3DV,B3DV,XSTD,FFTW_PLANS)
          ENDIF
C
          CALL SHIFT(IPAD*NSAM,C3DV)
C
      	ENDIF

      	WRITE(*,*)' 3D WEIGHTS FILE FOR OUTPUT?'
     	READ(*,7006)FWEIGH

      	N1=JC
      	N2=NSAM
      	N3=NSAM
        IF ((FALL).AND.(.NOT.FDUMP)) THEN
      	  CALL IOPEN(FWEIGH,INSTAT,CFORM,N1,N2,N3,'NEW',ASYM,
     .			PSIZE,VX)
        ENDIF

	UTD=0.0
        NUTD=0
        STDD=0.0

	END
