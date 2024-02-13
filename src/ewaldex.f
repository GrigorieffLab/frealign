C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      COMPLEX FUNCTION EWALDEX(NSAM,IRAD,A3DF,SINCLUT,IPAD,
     +                         I,J,DM,THET,CTFR,CTFL)
C**************************************************************************
C   Extraction of voxel on Ewald sphere from reference 3d volume
C   Uses function AINTERPO3DS or AINTERPO3DBIG.
C   Used in PEXTRACT, PRESB, CC3 and CC3M, if IEWALD is set
C   M.WOLF, Aug.2004
C ** MW still have to check, if all calling routines need THEATATR/AMAGP(K)
C**************************************************************************

      IMPLICIT NONE

      INTEGER NSAM,IRAD,IPAD,I,J
      REAL X3,Y3,Z3,SINCLUT(*),DM(*),THET
      COMPLEX A3DF(*),AINTERPO3DBIG,AINTERPO3DS
C MW Curvature of the Ewald sphere
      REAL THETAH,GPIX,X,Y,Z,ROBS,IOBS
      COMPLEX AINPO1,AINPO2,CTFR,CTFL

C**************************************************************************

      GPIX=SQRT(REAL(I**2+J**2))   ! length of resol.vector g in pixel
      THETAH=GPIX*THET/2           ! THET=(WL/(PSIZE*NSAM))/AMAGP
      X=I*COS(THETAH)              ! THETAH=scattering angle/2.
      Y=J*COS(THETAH)
      Z=GPIX*SIN(THETAH)
      X3=DM(1)*X+DM(4)*Y+DM(7)*Z
      Y3=DM(2)*X+DM(5)*Y+DM(8)*Z
      Z3=DM(3)*X+DM(6)*Y+DM(9)*Z
      IF (IRAD.EQ.0) THEN          ! AINPO1 is first beam
        AINPO1=AINTERPO3DBIG(NSAM,IPAD,A3DF,X3,Y3,Z3)
      ELSE
        AINPO1=AINTERPO3DS(IPAD*NSAM,IRAD,A3DF,
     +                      X3,Y3,Z3,SINCLUT,IPAD)
      ENDIF                 
      X3=DM(1)*X+DM(4)*Y-DM(7)*Z   ! AINPO2 is second beam
      Y3=DM(2)*X+DM(5)*Y-DM(8)*Z   ! = conjg(x,y,-z)
      Z3=DM(3)*X+DM(6)*Y-DM(9)*Z   
      IF (IRAD.EQ.0) THEN
        AINPO2=AINTERPO3DBIG(NSAM,IPAD,A3DF,X3,Y3,Z3)
      ELSE
        AINPO2=AINTERPO3DS(IPAD*NSAM,IRAD,A3DF,
     +                           X3,Y3,Z3,SINCLUT,IPAD)
      ENDIF

      EWALDEX=(AINPO1*CTFR+AINPO2*CONJG(CTFL))     ! combine beams

      RETURN
      END
