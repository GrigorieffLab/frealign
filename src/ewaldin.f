C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE EWALDIN(NSAM,IRAD,A3DF,SINCLUT,IPAD,
     +                   I,J,DM,THET,XYZ,FOR,FOL,FOBS,PREL,
     +                   IEWALD,CTFR,CTFL)
C**************************************************************************
C   Reference based derivation of left and right beam for Ewald insertion
C   Uses function AINTERPO3DS or AINTERPO3DBIG.
C   Used in PINSERT2. A3DF (reference volume) is C3DF elsewhere.
C   M.WOLF, Okt.2004
C**************************************************************************

      IMPLICIT NONE

      INTEGER NSAM,IRAD,IPAD,I,J,IEWALD
      REAL X3,Y3,Z3,SINCLUT(*),DM(*),THET
      REAL CS,WL,WGH1,WGH2,DF1,DF2,ANGAST
      COMPLEX A3DF(*),AINTERPO3DBIG,AINTERPO3DS
C Ewald sphere
      REAL THETAH,GPIX,X,Y,Z,XYZ(*),AOBS,AREF,PREL
      COMPLEX FOBS,FREF
      COMPLEX FRR,FRL,FOR,FOL,FREL,CTFR,CTFL
C**************************************************************************

C First get reference amplitude and phase for left and right beam
C from previous 3D model

      GPIX=SQRT(REAL(I**2+J**2))             ! length of resol.vector g in pixel
      THETAH=GPIX*THET/2                     ! THET=(WL/(PSIZE*NSAM))/AMAGP
      X=I*COS(THETAH)                        ! THETAH=scattering angle/2.
      Y=J*COS(THETAH)
      Z=GPIX*SIN(THETAH)
      XYZ(1)=DM(1)*X+DM(4)*Y+DM(7)*Z
      XYZ(2)=DM(2)*X+DM(5)*Y+DM(8)*Z
      XYZ(3)=DM(3)*X+DM(6)*Y+DM(9)*Z
      XYZ(4)=DM(1)*X+DM(4)*Y-DM(7)*Z
      XYZ(5)=DM(2)*X+DM(5)*Y-DM(8)*Z
      XYZ(6)=DM(3)*X+DM(6)*Y-DM(9)*Z
      IF (ABS(IEWALD).GT.1) THEN
        X3=XYZ(1)
        Y3=XYZ(2)
        Z3=XYZ(3)
        IF (IRAD.EQ.0) THEN                      ! FRR is first beam
          FRR=AINTERPO3DBIG(NSAM,IPAD,A3DF,X3,Y3,Z3)
        ELSE
          FRR=AINTERPO3DS(IPAD*NSAM,IRAD,A3DF,
     +                      X3,Y3,Z3,SINCLUT,IPAD)
        ENDIF
        X3=XYZ(4)
        Y3=XYZ(5)
        Z3=XYZ(6)
        IF (IRAD.EQ.0) THEN
          FRL=AINTERPO3DBIG(NSAM,IPAD,A3DF,X3,Y3,Z3)
        ELSE
          FRL=AINTERPO3DS(IPAD*NSAM,IRAD,A3DF,
     +                           X3,Y3,Z3,SINCLUT,IPAD)
        ENDIF

C we only use the relative phase between the two extracted beams AINPO1/2 to
C reconstitute the unknowns FOR and FOL while maintaining amplitude and phase
C of the observed sample FOBS.

        AOBS=CABS(FOBS)                          ! amplitude, observed sf
        FREF=FRR*CTFR+FRL*CONJG(CTFL)            ! beam combination acc. to  DDR
        AREF=CABS(FREF)                          ! amplitude, reference sf
        PREL=0.0

        IF (AOBS.EQ.0.0) THEN
          FOR=CMPLX(0.0,0.0)
          FOL=CMPLX(0.0,0.0)
        ELSE
         IF (AREF.NE.0.0) THEN
          FREL=FOBS/FREF*CABS(CTFR+CONJG(CTFL))**2            ! difference vector
          PREL=ABS(ATAN2(AIMAG(FREL),REAL(FREL))) ! phase of diff. vector
          FOR=FRR*FREL                            ! right derived sf
          FOL=FRL*FREL                            ! left derived sf
         ELSE
          FOR=FOBS*CONJG(CTFR)                    ! without reference just use
          FOL=FOBS*CTFL                           ! the observed values 
         ENDIF
        ENDIF
C
      ELSE
C
        FOR=FOBS*CONJG(CTFR)                      ! without reference just use
        FOL=FOBS*CTFL                             ! the observed values 
C
      ENDIF
  
      RETURN

      END
