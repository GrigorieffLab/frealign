C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      REAL FUNCTION TRILINMAP(NSAM,A3DV,X,Y,Z)
C**************************************************************************
C   Trilinear interpolation of 3D realspace map used in BEAUTIFY
C   Used in BEAUTIFY.
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,LS,LT,NS,NT,MS,MT
      INTEGER L,M,N,LL,MM,NN,IS,ID
      REAL X,Y,Z,WGT,VMIN,VMAX
      REAL A3DV(*),SAMP
C**************************************************************************
      VMIN=1.0
      VMAX=FLOAT(NSAM) + 0.0001
      IF((X.LT.VMIN).OR.(Y.LT.VMIN).OR.(Z.LT.VMIN).OR.
     .          (X.GT.VMAX).OR.(Y.GT.VMAX).OR.(Z.GT.VMAX)) THEN
        write(*,10) X,Y,Z,VMIN,VMAX
10      FORMAT(' TRILINMAP: X,Y or Z too small or too big'/
     .          '   X,Y,Z =',3F12.5,'    MIN,MAX =',2F12.5)
      	STOP ' ERROR: This should never occur. Please report bug'
      ENDIF

      	SAMP=0.0
      LS=INT(X)
      LS=MAX(LS,1)
      LT=INT(X)+1
      LT=MIN(LT,NSAM)
      	MS=INT(Y)
      	MS=MAX(MS,1)
      	MT=INT(Y)+1
      	MT=MIN(MT,NSAM)
      NS=INT(Z)
      NS=MAX(NS,1)
      NT=INT(Z)+1
      NT=MIN(NT,NSAM)

      DO 40 L=LS,LT
        DO 40 M=MS,MT
          DO 40 N=NS,NT
C	    TRILINEAR INTERPOLATION.....
      	    WGT=(1.0-ABS(X-L))*(1.0-ABS(Y-M))*(1.0-ABS(Z-N))
      		ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
                SAMP=SAMP+A3DV(ID)*WGT
40    CONTINUE
      TRILINMAP=SAMP
      RETURN
      END
