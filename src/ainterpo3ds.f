C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      COMPLEX FUNCTION AINTERPO3DS(NSAM,IRAD,A3DF,
     .                             XX,YY,ZZ,SINCLUT,IPAD)
C**************************************************************************
C     Uses function  BOX FT_LUT.
C     Used in funct PRES, and subr PRESB, CCP, PEXTRACT and APPLYSYMC.	
C**************************************************************************
      IMPLICIT NONE
      INTEGER NSAM,IRAD,LS,LT,NS,NT,MS,MT,JC
      INTEGER L,M,N,MM1,NN,IPAD,SIGN,I,J
      REAL X,Y,Z,PI,BOXFT,BOXFTV,BOXFT_LUT,SINCLUT(*)
      REAL A3(8), A2(8), A1(8)
      INTEGER P2(8), P3(8), LL(8), MM(8), ID(8)
      REAL XX,YY,ZZ,P1,RBUF(8)
      PARAMETER  (PI=3.1415926535897)
      COMPLEX A3DF(*),SAMP,CBUF(8), R(8)
C**************************************************************************
      SAMP=CMPLX(0.0,0.0)
      JC=NSAM/2+1
      IF (XX.LT.0.0) THEN
        X=-IPAD*XX
        Y=-IPAD*YY
        Z=-IPAD*ZZ
        SIGN=-1
      ELSE
        X=IPAD*XX
        Y=IPAD*YY
        Z=IPAD*ZZ
        SIGN=1
      ENDIF
      LS=INT(X)
      LT=LS+1
      IF (LT.GE.JC) LT=JC-1
      MS=INT(Y)
      MT=MS+1
      IF (Y.LT.0.0) THEN
        MS=MS-1
        MT=MT-1
      ENDIF
      IF (MS.LE.-JC) MS=-JC+1
      IF (MT.GE.JC) MT=JC-1
      NS=INT(Z)
      NT=NS+1
      IF (Z.LT.0.0) THEN
        NS=NS-1
        NT=NT-1
      ENDIF
      IF (NS.LE.-JC) NS=-JC+1
      IF (NT.GE.JC) NT=JC-1

      A3(1:4) = 1.0 - ABS(Z-NS)
      A3(5:8) = 1.0 - ABS(Z-NT)
      A2(1:2) = 1.0 - ABS(Y-MS)
      A2(3:4) = 1.0 - ABS(Y-MT)
      A2(5:6) = 1.0 - ABS(Y-MS)
      A2(7:8) = 1.0 - ABS(Y-MT)
      A1(1:7:2) = 1.0 - ABS(X-LS)
      A1(2:8:2) = 1.0 - ABS(X-LT)

      NN = NS + 1
      IF (NN.LT.1) NN=NN+NSAM
      P2(1:4) = (NN - 1)
      NN = NT + 1
      IF (NN.LT.1) NN=NN+NSAM
      P2(5:8) = (NN - 1)
      P2 = P2 * NSAM
      MM1 = MS + 1
      IF (MM1.LT.1) MM1=MM1+NSAM
      MM(1:2) = (MM1 - 1)
      MM(5:6) = (MM1 - 1)
      MM1 = MT + 1
      IF (MM1.LT.1) MM1=MM1+NSAM
      MM(3:4) = (MM1 - 1)
      MM(7:8) = (MM1 - 1)
      P3 = JC * (MM + P2)

      RBUF = A1 * A2 * A3

      LL(1:7:2) = LS + 1
      LL(2:8:2) = LT + 1
      ID = LL + P3
      CBUF = A3DF(ID)

      R = CBUF * RBUF
      SAMP = SUM(R)

      IF (SIGN.LT.0) THEN
        AINTERPO3DS=CONJG(SAMP)
      ELSE
        AINTERPO3DS=SAMP
      ENDIF
      RETURN
      END
