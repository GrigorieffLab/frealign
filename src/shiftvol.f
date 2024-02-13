C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------------------------
	SUBROUTINE SHIFTVOL(NSAM,STD,VTD,JC,RLIM,NSAMH,IREDUN,
     +			    A3DF,D3DF,S3DF,V3DF,ASUM,
     +                      QCP,NNSTAT,
     +                      IFSC,PSSNR,FPART,FSC,FMASK)
C -----------------------------------------------------------------------------
C       SHIFT VOLUME INTO CENTRE
C       Used in Frealign.
C -----------------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER L,LL,M,MM,N,NN,JC,NSAM,WC1,NSAMH,ID
        INTEGER NUTD,ISUM,IREDUN,JRAD2,IRREC2,NNSTAT,IFSC,I
        REAL    PSHFTR,UTD,TMP,RLIM,FPART,A1,FMASK
	DOUBLEPRECISION STD,VTD
	COMPLEX A3DF(*),D3DF(*)
	REAL    S3DF(*),V3DF(*),QCP(*),PSSNR(*)
	REAL    ASUM(*),FSC(*)
C -----------------------------------------------------------------------------
        STD=STD+VTD
        IRREC2=(INT(REAL(NSAM)*RLIM))**2
        DO 61 L=1,JC
          LL=L-1
          DO 61 M=1,NSAM
            MM=M-1
            IF (MM.GE.JC) MM=MM-NSAM
            DO 61 N=1,NSAM
              NN=N-1
              IF (NN.GE.JC) NN=NN-NSAM
              JRAD2=LL**2+MM**2+NN**2
              ISUM=(LL+MM+NN)
              PSHFTR=1.0
              IF (MOD(ISUM,2).NE.0) PSHFTR=-1.0
              I=INT(SQRT(REAL(LL**2+MM**2+NN**2))+0.5)+1
                ID=L+JC*((M-1)+NSAM*(N-1))
                IF (IFSC.EQ.0) THEN
                  A3DF(ID)=A3DF(ID)+D3DF(ID)
                  IF (NNSTAT.NE.0) THEN
                    IF (QCP(ID).NE.0.0) THEN
                      TMP=ASUM(ID)/QCP(ID)-1.0
                    ELSE
                      TMP=0.0
                    ENDIF
                    IF (TMP.GT.0.0) THEN
                      D3DF(ID)=CMPLX(TMP,0.0)*PSHFTR
                    ELSE
                      D3DF(ID)=CMPLX(0.0,0.0)
                    ENDIF
                  ENDIF
                  S3DF(ID)=S3DF(ID)+V3DF(ID)
                ENDIF
                IF(JRAD2.LE.IRREC2) THEN
                  IF (FPART.EQ.0.0) THEN
                    A3DF(ID)=A3DF(ID)*PSHFTR/(S3DF(ID)+STD*IREDUN)
                  ELSEIF ((PSSNR(I).NE.0.0).AND.(S3DF(ID).NE.0.0)) THEN
                    IF (ABS(FSC(I)).LT.1.0) THEN
                      A1=2.0*ABS(FSC(I))/(1.0-ABS(FSC(I)))
                    ELSE
                      A1=1000.0
                    ENDIF
C                    A3DF(ID)=A3DF(ID)*PSHFTR/(S3DF(ID)+1.0/PSSNR(I))
                     A3DF(ID)=A3DF(ID)*PSHFTR/(S3DF(ID)+STD*IREDUN)
     +                        /(1.0+FPART/A1/FMASK)
                  ELSE
                    A3DF(ID)=CMPLX(0.0,0.0)
                  ENDIF
                ELSE
                  A3DF(ID)=CMPLX(0.0,0.0)
                ENDIF
61      CONTINUE
	RETURN
	END
