C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PHASEFLIP(NSAM,SPEC,CTFF)
C**************************************************************************
C   Phase-flips the input image FFT (SPEC) according to the CTF (CTFF).
C   Used in SIGMA.
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,JC,I,J,II,JJ,ID,IC
      REAL CTFAMP
      COMPLEX SPEC(*),CTFF(*),CTF
C**************************************************************************
      JC=NSAM/2+1
      IC=NSAM*JC
C
      DO 30 I=0,JC-1
      	II=I+1
      	DO 30 J=-JC+1,JC-1
      	  JJ=J+1
      	  IF (JJ.LT.1) JJ=JJ+NSAM
            ID=II+JC*(JJ-1)
            CTF=CTFF(ID)+CONJG(CTFF(ID+IC))
            CTFAMP=CABS(CTF)
            IF (CTFAMP.NE.0.0) SPEC(ID)=SPEC(ID)*CTF/CTFAMP
30    CONTINUE
C
      RETURN
      END
