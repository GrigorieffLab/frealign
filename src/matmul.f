C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE MATMUL(A,B,C)
C**************************************************************************
C     3x3 matrix multiplication
C     Used in APPLYSYMC,APPLYSYMR, BEAUTIFY and CHECKSYM.
C**************************************************************************
      IMPLICIT NONE
      INTEGER I,J
      REAL A(3,3),B(3,3),C(3,3),T(3,3)
C**************************************************************************
      DO 10 I=1,3
        DO 10 J=1,3
      	  T(I,J)=A(I,1)*B(1,J)+A(I,2)*B(2,J)+A(I,3)*B(3,J)
10    CONTINUE
      DO 20 I=1,3
        DO 20 J=1,3
      	  C(I,J)=T(I,J)
20    CONTINUE
      RETURN
      END
