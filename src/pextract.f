C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PEXTRACT(NSAM,IRAD,OUTC,A3DF,
     +			  PHI,THETA,PSI,SINCLUT,IPAD,THET,
     +                    IEWALD,CTFF,SM)
C**************************************************************************
C   Extracts projection from a 3D transform (A3DF) and places
C   the resulting transform points into OUTC.  At present (6.2.99)
C   Uses function AINTERPO3DS or AINTERPO3DBIG.
C   Used in MATCH and LMAIN.
C   MW uses function EWALDEX if IEWALD is true
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,JC,I,J,II,JJ,IRAD,ID,NSAMH,IR2LIM,IR2,IPAD,IC
      REAL DM(9),PHI,THETA,PSI,CPHI,SPHI,CTHE,STHE,CPSI,SPSI
      REAL X3,Y3,Z3,SINCLUT(*),SM(9)
      COMPLEX OUTC(*),AINPO
      COMPLEX A3DF(*),AINTERPO3DBIG,AINTERPO3DS
C MW Ewald sphere correction
      REAL THET
      COMPLEX EWALDEX,CTFF(*),CTFR,CTFL
      INTEGER IEWALD
C**************************************************************************
      JC=NSAM/2+1
      IC=NSAM*JC
      NSAMH=NSAM/2
      IR2LIM=NSAMH**2
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      CTHE=COS(THETA)
      STHE=SIN(THETA)
      CPSI=COS(PSI)
      SPSI=SIN(PSI)

      DM(1)=CPHI*CTHE*CPSI-SPHI*SPSI
      DM(2)=SPHI*CTHE*CPSI+CPHI*SPSI
      DM(3)=-STHE*CPSI
      DM(4)=-CPHI*CTHE*SPSI-SPHI*CPSI
      DM(5)=-SPHI*CTHE*SPSI+CPHI*CPSI
      DM(6)=STHE*SPSI
C MW next lines uncommented for Ewald sphere correction
      DM(7)=STHE*CPHI
      DM(8)=STHE*SPHI
      DM(9)=CTHE
      CALL MATMUL_T(DM,SM,DM)

      DO 30 I=0,JC-1
      	II=I+1
      	DO 30 J=-JC+1,JC-1
          IR2=I**2+J**2
          IF (IR2.LE.IR2LIM) THEN
      	    JJ=J+1
      	    IF (JJ.LT.1) JJ=JJ+NSAM
      	      ID=II+JC*(JJ-1)
      	      CTFR=CTFF(ID)
      	      CTFL=CTFF(ID+IC)
            IF (IEWALD.NE.0) THEN
              IF (IEWALD.LT.0) THEN
                CTFR=CONJG(CTFR)
                CTFL=CONJG(CTFL)
              ENDIF
C MW call ewald sphere corrected extraction
              AINPO=EWALDEX(NSAM,IRAD,A3DF,SINCLUT,IPAD,
     +                         I,J,DM,THET,CTFR,CTFL)
            ELSE
      	      X3=DM(1)*I+DM(4)*J
      	      Y3=DM(2)*I+DM(5)*J
      	      Z3=DM(3)*I+DM(6)*J
              IF (IRAD.EQ.0) THEN
                AINPO=AINTERPO3DBIG(NSAM,IPAD,A3DF,X3,Y3,Z3)
              ELSE
                AINPO=AINTERPO3DS(IPAD*NSAM,IRAD,A3DF,
     +                            X3,Y3,Z3,SINCLUT,IPAD)
              ENDIF
                AINPO=AINPO*(CTFR+CONJG(CTFL))
            ENDIF
C      	    IF (II.NE.JC) THEN
      	      OUTC(ID)=AINPO
C      	    ELSE
C      	      OUTQ(JJ)=AINPO
C      	    ENDIF
          ELSE
      	    JJ=J+1
      	    IF (JJ.LT.1) JJ=JJ+NSAM
C      	    IF (II.NE.JC) THEN
      	      ID=II+JC*(JJ-1)
      	      OUTC(ID)=CMPLX(0.0,0.0)
C      	    ELSE
C      	      OUTQ(JJ)=CMPLX(0.0,0.0)
C      	    ENDIF
          ENDIF
30    CONTINUE

      RETURN
      END
