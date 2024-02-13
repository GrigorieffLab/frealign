C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      COMPLEX FUNCTION AINTERPO3DBIG(NSAMA,IPAD,A3DF,X,Y,Z)
C**************************************************************************
C     Used in funct PRES, and subr PRESB, CCP, PEXTRACT.	
C**************************************************************************
      IMPLICIT NONE
      INTEGER NSAM,JC,IPAD,NSAMA
      INTEGER L,M,N,LL,MM,NN,ID
      REAL X,Y,Z
      COMPLEX A3DF(*),SAMP
C**************************************************************************
      NSAM=IPAD*NSAMA
      SAMP=CMPLX(0.0,0.0)
      JC=NSAM/2+1
      L=NINT(IPAD*X)
      M=NINT(IPAD*Y)
      N=NINT(IPAD*Z)
      IF (L.GE.0) THEN
        LL=L+1
        MM=M+1
        IF (MM.LT.1) MM=MM+NSAM
        NN=N+1
        IF (NN.LT.1) NN=NN+NSAM
      	ID=LL+JC*((MM-1)+NSAM*(NN-1))
        SAMP=SAMP+A3DF(ID)
      ELSE
        LL=-L+1
        MM=-M+1
        IF (MM.LT.1) MM=MM+NSAM
        NN=-N+1
        IF (NN.LT.1) NN=NN+NSAM
      	ID=LL+JC*((MM-1)+NSAM*(NN-1))
        SAMP=SAMP+CONJG(A3DF(ID))
      ENDIF
      AINTERPO3DBIG=SAMP
      RETURN
      END
