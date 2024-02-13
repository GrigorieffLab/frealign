C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------------------------
	SUBROUTINE SHIFT(NSAM,C3DF)
C -----------------------------------------------------------------------------
C	  SHIFT VOLUME TO ORIGIN
C	  C3D     - 3D map I/O
C	  U3D     - 3D weights I/O
C         Used in CARD13AND14
C -----------------------------------------------------------------------------
          IMPLICIT NONE
C
          INTEGER NSAM,JC,L,LL,M,MM,N,NN,ID,ISUM
	  REAL PSHFTR
          COMPLEX C3DF(*)
C -----------------------------------------------------------------------------
C          NSAMH=NSAM/2
          JC=NSAM/2+1
          DO 80 L=1,JC
      	    LL=L-1
      	    DO 80 M=1,NSAM
      	      MM=M-1
      	      IF (MM.GE.JC) MM=MM-NSAM
      	      DO 80 N=1,NSAM
      	        NN=N-1
      	        IF (NN.GE.JC) NN=NN-NSAM
      		ISUM=(LL+MM+NN)
      		PSHFTR=1.0
      		IF (MOD(ISUM,2).NE.0) PSHFTR=-1.0
      	          ID=L+JC*((M-1)+NSAM*(N-1))
      		  C3DF(ID)=C3DF(ID)*PSHFTR
80	  CONTINUE
	RETURN
	END
