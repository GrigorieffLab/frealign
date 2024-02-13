C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      REAL FUNCTION FSH(SIG2,X,Y,XM,YM,SX,SY)
C**************************************************************************
C  Calculates probability distribution for shifts
C  Used in CALCFX and PREFINE
C**************************************************************************
      IMPLICIT NONE
      REAL XM,YM,SX,SY,X,Y,SIG2
C**************************************************************************
C
      FSH=SIG2*(-(X-XM)**2/2.0/SX**2-(Y-YM)**2/2.0/SY**2)
C
      RETURN
      END
