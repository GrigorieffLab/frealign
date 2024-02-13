C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------------------------
      SUBROUTINE SMOOTH_PSSNR(NSAM,S3DF,V3DF,
     +                        RLIM,PSSNR,SUM,ISUM)
C -----------------------------------------------------------------------------
C       Smooth PSSNR curve to reduce noise
C       Used in Frealign.
C -----------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER L,LL,M,MM,N,NN,JC,NSAM,IS,ID
      INTEGER JRAD2,IRREC2,I,ISUM(*)
      REAL S3DF(*),V3DF(*),PSSNR(*),RLIM,SUM(*)
C -----------------------------------------------------------------------------
      JC=NSAM/2+1
      IRREC2=(INT(REAL(NSAM)*RLIM))**2
C
      DO 10 I=1,NSAM
        SUM(I)=0.0
        ISUM(I)=0
10    CONTINUE
C
      DO 61 L=1,JC
        LL=L-1
        DO 61 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          DO 61 N=1,NSAM
            NN=N-1
            IF (NN.GE.JC) NN=NN-NSAM
            JRAD2=LL**2+MM**2+NN**2
            IF(JRAD2.LE.IRREC2) THEN
              I=INT(SQRT(REAL(LL**2+MM**2+NN**2))+0.5)+1
                ID=L+JC*((M-1)+NSAM*(N-1))
                SUM(I)=SUM(I)+S3DF(ID)+V3DF(ID)
61    CONTINUE
C
      RETURN
      END
