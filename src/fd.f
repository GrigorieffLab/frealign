C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      REAL FUNCTION FD(SIG2,D1,D2,DFSTD)
C**************************************************************************
C  Calculates probability distribution for shifts
C  Used in CALCFX and PREFINE
C**************************************************************************
      IMPLICIT NONE
      REAL SD,D1,D2,SIG2,DFSTD
C**************************************************************************
C
      FD=SIG2*(-D1**2/2.0/DFSTD**2-D2**2/2.0/DFSTD**2)
C
      RETURN
      END
