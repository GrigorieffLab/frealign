C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C ---------------------------------------------------------------------------------        
      REAL FUNCTION FPART_CC(JC,NSAM,NSAMH,
     +			B3DF,C3DF,D3DF,S3DF,V3DF,
     +                  ASYM,RI,RIC,B3DV,C3DV,
     +                  FPART,SSNR,MAVE,MSTD,FSC,FFTW_PLANS)
C ---------------------------------------------------------------------------------        
C       Find fparticle to scale SSNR for Wiener filter
C       Used in Frealign.
C ---------------------------------------------------------------------------------        
C
      USE ISO_C_BINDING
      USE FFTW33
C
      IMPLICIT NONE
C
      INTEGER L,LL,M,MM,N,NN,JC,NSAM,NSAMH,IS,ID
      INTEGER NSAM3,ISUM,I,J
      REAL PSHFTR,MAVE,MSTD,RI,RIC,SSNR(*)
      DOUBLE PRECISION SUM1,SUM2,SUM3
      REAL S3DF(*),B3DV(*),C3DV(*),FPART
      REAL V3DF(*),MAVE2,MSTD2,FSC(*)
      COMPLEX B3DF(*),C3DF(*),D3DF(*)
      CHARACTER ASYM*3
      TYPE(C_PTR) FFTW_PLANS(*)
C ---------------------------------------------------------------------------------        
C
      NSAM3=NSAM*NSAM*NSAM
C
      DO 70 N=1,NSAM
        NN=N-1
        IF (NN.GE.JC) NN=NN-NSAM
        DO 70 M=1,NSAM
          MM=M-1
          IF (MM.GE.JC) MM=MM-NSAM
          DO 70 L=1,JC
            LL=L-1
            ISUM=(LL+MM+NN)
            PSHFTR=1.0
            IF (MOD(ISUM,2).NE.0) PSHFTR=-1.0
            I=INT(SQRT(REAL(LL**2+MM**2+NN**2))+0.5)+1
              ID=L+JC*((M-1)+NSAM*(N-1))
              IF ((SSNR(I).NE.0.0).AND.(FSC(I).GT.0.5)) THEN
                C3DF(ID)=D3DF(ID)*PSHFTR
     +            /(V3DF(ID)+FPART/SSNR(I))/NSAM3
              ELSE
                C3DF(ID)=CMPLX(0.0,0.0)
              ENDIF
70    CONTINUE
C
      CALL FFTW_BWD(C3DF,C3DF,FFTW_PLANS(4))
C
      ISUM=0
      SUM1=0.0
      SUM2=0.0
      DO 20 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 20 I=1,NSAM
          ID=IS+I
          IF (B3DV(ID).NE.0.0) THEN
            SUM1=SUM1+C3DV(ID)
            SUM2=SUM2+C3DV(ID)**2
            ISUM=ISUM+1
          ENDIF
20    CONTINUE
      IF (ISUM.NE.0) THEN
        SUM1=SUM1/ISUM
        SUM2=SUM2/ISUM
      ENDIF
      MAVE2=SUM1
      MSTD2=SQRT(ABS(SUM2-SUM1**2))
C
      SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      DO 10 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 10 I=1,NSAM
          ID=IS+I
          IF (B3DV(ID).NE.0.0) THEN
            SUM1=SUM1+(B3DV(ID)-MAVE)*(C3DV(ID)-MAVE2)
            SUM2=SUM2+(B3DV(ID)-MAVE)**2
            SUM3=SUM3+(C3DV(ID)-MAVE2)**2
          ENDIF
10    CONTINUE
      FPART_CC=SUM1/SQRT(SUM2*SUM3)
C
      RETURN
      END
