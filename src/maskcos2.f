C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE MASKCOS2(NSAM,OUTD,RI2,RI3,RIH,HALFW,AMAGP)
C**************************************************************************
C  Applies cosine-edged mask to 2D image
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,ID,I,J,II,JJ,NMID,ISUM,NSAMH
      REAL OUTD(*),FI1,AMAGP,CRAD2,SRAD2,RIH
      REAL RI2,RI3,HALFW,PI,SUM,EDGE
      REAL RIHP,RI2P,RI3P
      PARAMETER (PI=3.1415926535897)
C**************************************************************************
      RIHP=RIH/AMAGP
      RI2P=RI2/AMAGP/AMAGP
      RI3P=RI3/AMAGP/AMAGP
      NMID=1
      NSAMH=NSAM/2
C
      SUM=0.0
      ISUM=0
      DO 20 I=1,NSAM
        II=I
        IF (II.GT.NSAMH) II=II-NSAM
        FI1=REAL(II-NMID)**2
        DO 20 J=1,NSAM
          JJ=J
          IF (JJ.GT.NSAMH) JJ=JJ-NSAM
      	  ID=J+(NSAM+2)*(I-1)
          CRAD2=REAL(JJ-NMID)**2+FI1
          IF ((CRAD2.GT.RI2P).AND.(CRAD2.LE.RI3P)) THEN
            SUM=SUM+OUTD(ID)
            ISUM=ISUM+1
          ENDIF
20    CONTINUE
      IF (ISUM.NE.0) SUM=SUM/ISUM
C
C     MASK PARTICLE, COSINE EDGE
      DO 21 I=1,NSAM
        II=I
        IF (II.GT.NSAMH) II=II-NSAM
        FI1=REAL(II-NMID)**2
        DO 21 J=1,NSAM
          JJ=J
          IF (JJ.GT.NSAMH) JJ=JJ-NSAM
      	  ID=J+(NSAM+2)*(I-1)
          CRAD2=REAL(JJ-NMID)**2+FI1
          IF ((CRAD2.GT.RI2P).AND.(CRAD2.LE.RI3P)) THEN
            SRAD2=SQRT(CRAD2)
            EDGE=(1.0+COS(PI*(SRAD2-RIHP)/HALFW))/2.0
            OUTD(ID)=OUTD(ID)*EDGE+(1.0-EDGE)*SUM
      	  ELSEIF (CRAD2.GT.RI3P) THEN
            OUTD(ID)=SUM
          ENDIF
21    CONTINUE
C
      RETURN
      END
