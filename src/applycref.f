C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------------------------
      SUBROUTINE APPLYCREF(NSAM,A3DF,FSC)
C -----------------------------------------------------------------------------
C       APPLY WINER FILTER TO VOLUME
C       Used in Frealign
C -----------------------------------------------------------------------------
        IMPLICIT NONE
C
        INTEGER L,LL,M,MM,N,NN,JC,NSAM,ID
        INTEGER JRAD2,ILIM
        REAL FSC(*),F
        COMPLEX A3DF(*)
C -----------------------------------------------------------------------------
C        NSAMH=NSAM/2
        JC=NSAM/2+1
        ILIM=(NSAM/2)**2
        DO 61 L=1,JC
          LL=L-1
          DO 61 M=1,NSAM
            MM=M-1
            IF (MM.GE.JC) MM=MM-NSAM
            DO 61 N=1,NSAM
              NN=N-1
              IF (NN.GE.JC) NN=NN-NSAM
              JRAD2=LL**2+MM**2+NN**2
              F=0.0
              IF (JRAD2.LE.ILIM) F=ABS(FSC(INT(SQRT(REAL(JRAD2)))+1))
                ID=L+JC*((M-1)+NSAM*(N-1))
                IF(JRAD2.LE.ILIM) THEN
                  A3DF(ID)=A3DF(ID)
     .              *SQRT(2.0*F/(1.0+F))
                ELSE
                  A3DF(ID)=CMPLX(0.0,0.0)
                ENDIF
61      CONTINUE
        RETURN
        END

