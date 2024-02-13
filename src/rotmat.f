C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE ROTMAT(PHI,THETA,PSI,AMAG,RM)
C**************************************************************************
C   Sets up rotation matrix
C**************************************************************************
      IMPLICIT NONE
C
      REAL PHI,THETA,PSI,AMAG
      REAL RM(9),CPHI,SPHI,CTHE,STHE,CPSI,SPSI
C**************************************************************************
C
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      CTHE=COS(THETA)
      STHE=SIN(THETA)
      CPSI=COS(PSI)
      SPSI=SIN(PSI)
C
      RM(1)=(CPHI*CTHE*CPSI-SPHI*SPSI)/ABS(AMAG)
      RM(2)=(SPHI*CTHE*CPSI+CPHI*SPSI)/ABS(AMAG)
      RM(3)=-STHE*CPSI/ABS(AMAG)
      RM(4)=(-CPHI*CTHE*SPSI-SPHI*CPSI)/ABS(AMAG)
      RM(5)=(-SPHI*CTHE*SPSI+CPHI*CPSI)/ABS(AMAG)
      RM(6)=STHE*SPSI/ABS(AMAG)
      RM(7)=STHE*CPHI/ABS(AMAG)
      RM(8)=STHE*SPHI/ABS(AMAG)
      RM(9)=CTHE/ABS(AMAG)
C
      RETURN
      END
