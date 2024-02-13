C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------------
      SUBROUTINE COMBINE_ARRAYS(NSAM,ASUM,VSUM,PSUM,KSUM,
     +           A3DF,S3DF,STD,D3DF,V3DF,
     +           VTD,NN1,NNSTAT,NKSUM,IFSC,IMP,NVOL,
     +           NVOL1,NVOL2,NKSUM1,TTD1,TTD2,NN2,NN3)
C -----------------------------------------------------------------
C  Combine parallel arrays populated by A3D3
C  Used in FREALIGN
C -----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NN1,NSAM,IS,NNSTAT,IFSC,IMP,I,JC
      INTEGER NVOL,NVOL1,NVOL2,NN2,NN3
      INTEGER KSUM(NN3*NNSTAT+1,NVOL)
      INTEGER*8 NKSUM,NKSUM1(*)
      REAL S3DF(NN3,NVOL1),V3DF(NN3,NVOL2)
      REAL ASUM(NN3*NNSTAT+1,NVOL)
      REAL VSUM(NN3*NNSTAT+1,NVOL)
      REAL PSUM(NN3*NNSTAT+1,NVOL)
      DOUBLE PRECISION STD,VTD,TTD1(*),TTD2(*)
      COMPLEX A3DF(NN3,NVOL1),D3DF(NN3,NVOL2)
C -----------------------------------------------------------------
      JC=NSAM/2+1
!$OMP PARALLEL DO
      DO 101 I=1,NSAM*NSAM*JC
        DO 101 IS=2,NVOL1
          A3DF(I,1)=A3DF(I,1)+A3DF(I,IS)
          S3DF(I,1)=S3DF(I,1)+S3DF(I,IS)
101   CONTINUE
      IF (IFSC.EQ.0) THEN
!$OMP PARALLEL DO
        DO 103 I=1,NSAM*NSAM*JC
          DO 103 IS=2,NVOL2
            D3DF(I,1)=D3DF(I,1)+D3DF(I,IS)
            V3DF(I,1)=V3DF(I,1)+V3DF(I,IS)
103     CONTINUE
      ENDIF
      IF (NNSTAT.NE.0) THEN
!$OMP PARALLEL DO
        DO 104 I=1,NSAM*NSAM*JC
          DO 104 IS=2,NVOL
            KSUM(I,1)=KSUM(I,1)+KSUM(I,IS)
            PSUM(I,1)=PSUM(I,1)+PSUM(I,IS)
            VSUM(I,1)=VSUM(I,1)+VSUM(I,IS)
            ASUM(I,1)=ASUM(I,1)+ASUM(I,IS)
104     CONTINUE
      ENDIF
      DO 810 I=1,IMP
        NKSUM=NKSUM+NKSUM1(I)
810   CONTINUE
      IF (NVOL1.NE.1) THEN
        DO 820 I=1,NVOL1
          STD=STD+TTD1(I)
          IF (IFSC.EQ.0) VTD=VTD+TTD2(I)
820     CONTINUE
      ELSE
        DO 821 I=1,IMP
          STD=STD+TTD1(I)
          IF (IFSC.EQ.0) VTD=VTD+TTD2(I)
821     CONTINUE
      ENDIF
C
      RETURN
      END
