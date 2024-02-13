C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE ROTATE2D(NSAM,IPADA,SPEC,
     +                    OUTC,PSI,BUFC,FFTW_PLANS)
C**************************************************************************
C     Used in function PSEARCH
C**************************************************************************
C
      USE ISO_C_BINDING
C
      IMPLICIT NONE
C
      INTEGER NSAM,JC4,IPAD,IPADA,NSAM4,I,J,II,JJ,ID4
      INTEGER L,M,LL,MM,ID,NSAMH,IRAD2,NSAMH4,ITEMP,JC
      REAL X,Y,PSI,DM(9),CPSI,SPSI
      COMPLEX SPEC(*),BUFC(*),OUTC(*)
      TYPE(C_PTR) FFTW_PLANS(*)
C**************************************************************************
C
      IPAD=ABS(IPADA)
      NSAM4=IPAD*NSAM
      NSAMH=NSAM/2
      JC=NSAM/2+1
      NSAMH4=NSAM4/2
      IRAD2=(NSAMH-1)**2
      JC4=NSAM4/2+1
      DO 10 I=1,NSAM*JC
        OUTC(I)=CMPLX(0.0,0.0)
10    CONTINUE
C
      CPSI=COS(PSI)
      SPSI=SIN(PSI)
C
      DM(1)=CPSI
      DM(2)=SPSI
      DM(4)=-SPSI
      DM(5)=CPSI
C
      IF (IPADA.LT.0) 
     +  CALL PAD(NSAM,IPAD,SPEC,BUFC,FFTW_PLANS)
C
      DO 20 I=0,NSAMH
        II=I+1
        DO 20 J=-NSAMH,NSAMH
          ITEMP=I**2+J**2
          IF (ITEMP.LT.IRAD2) THEN
            JJ=J+1
            IF (JJ.LT.1) JJ=JJ+NSAM
            ID=II+JC*(JJ-1)
            X=DM(1)*I+DM(4)*J
            Y=DM(2)*I+DM(5)*J
C
            L=NINT(IPAD*X)
            M=NINT(IPAD*Y)
            IF (L.GE.0) THEN
              LL=L+1
              MM=M+1
              IF (MM.LT.1) MM=MM+NSAM4
      	      ID4=LL+JC4*(MM-1)
              OUTC(ID)=BUFC(ID4)
            ELSE
              LL=-L+1
              MM=-M+1
              IF (MM.LT.1) MM=MM+NSAM4
      	      ID4=LL+JC4*(MM-1)
              OUTC(ID)=CONJG(BUFC(ID4))
            ENDIF
          ENDIF
20    CONTINUE
C
      RETURN
      END
