C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE MASKCOS_C(NSAM,OUTD,RI2,RI3,RIH,
     +                     HALFW,AMAGP,PSI)
C**************************************************************************
C  Applies cosine-edged mask to 2D image
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,ID,I,II,NMID,ISUM
      REAL OUTD(*),AMAGP,CRAD2,SRAD2,RIH
      REAL RI2,RI3,HALFW,PI,SUM,EDGE,X,Y
      REAL RIHP,RI2P,RI3P,PSI,CANG,SANG
      REAL RXI,RXH,RX2,RX3,XR2,SRX2
      PARAMETER (PI=3.1415926535897)
C**************************************************************************
      RIHP=RIH/AMAGP
      RI2P=RI2/AMAGP/AMAGP
      RI3P=RI3/AMAGP/AMAGP
      NMID=NSAM/2+1
      CANG=COS(PSI)
      SANG=SIN(PSI)
      RXI=NSAM/2-1-HALFW/2.0
      RXH=NSAM/2-1-3.0*HALFW/2.0
      RX2=RXH**2
      RX3=RXI**2
C
      SUM=0.0
      ISUM=0
      DO 20 I=1,NSAM
        DO 20 II=1,NSAM
C          X=-SANG*(I-NMID)+CANG*(II-NMID)
          Y=CANG*(I-NMID)+SANG*(II-NMID)
          CRAD2=Y**2
      	  ID=II+(NSAM+2)*(I-1)
          IF ((CRAD2.GT.RI2P).AND.(CRAD2.LE.RI3P)) THEN
            SUM=SUM+OUTD(ID)
            ISUM=ISUM+1
          ENDIF
20    CONTINUE
      IF (ISUM.NE.0) SUM=SUM/ISUM
C
C     MASK PARTICLE, COSINE EDGE
      DO 21 I=1,NSAM
        DO 21 II=1,NSAM
          X=-SANG*(I-NMID)+CANG*(II-NMID)
          Y=CANG*(I-NMID)+SANG*(II-NMID)
          XR2=X**2
          CRAD2=Y**2
      	  ID=II+(NSAM+2)*(I-1)
          IF ((CRAD2.GT.RI2P).AND.(CRAD2.LE.RI3P)) THEN
            SRAD2=SQRT(CRAD2)
            EDGE=(1.0+COS(PI*(SRAD2-RIHP)/HALFW))/2.0
            OUTD(ID)=OUTD(ID)*EDGE+(1.0-EDGE)*SUM
      	  ELSEIF (CRAD2.GT.RI3P) THEN
            OUTD(ID)=SUM
          ENDIF
cc          IF ((XR2.GT.RX2).AND.(XR2.LE.RX3)) THEN
cc            SRX2=SQRT(XR2)
cc            EDGE=(1.0+COS(PI*(SRX2-RXH)/HALFW))/2.0
cc            OUTD(ID)=OUTD(ID)*EDGE+(1.0-EDGE)*SUM
cc          ELSEIF (XR2.GT.RX3) THEN
cc            OUTD(ID)=SUM
cc          ENDIF
21    CONTINUE
C
      RETURN
      END
