C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PWEIGHT(NSAM,SPEC,OUTC,
     +          SSNR,CTFF,PWEIGHTS,RBUF,ICTF)
C**************************************************************************
C  Applies SNR weighting to input 2D image
C  Used in CTFAPPLY
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER NSAM,JC,I,J,II,JJ,IRAD,NSAMH,ICTF
      INTEGER MAXR12,MAXR22,ITEMP,IBIN,IC,ID
      REAL SSNR(*),SNR,WGT,PWEIGHTS(*),RBUF(*)
      COMPLEX SAMP,SPEC(*),OUTC(*)
      COMPLEX CTFF(*),CTFR,CTFL,Z
C**************************************************************************
      JC=NSAM/2+1
      IC=NSAM*JC
C      NSAMH=NSAM/2
      MAXR12=0
      MAXR22=JC**2
C
      DO 50 I=1,JC
        PWEIGHTS(I)=0.0
50    CONTINUE
C
      DO 10 I=1,2*JC
        RBUF(I)=0.0
10    CONTINUE
C
      DO 20 I=0,JC-1
        II=I+1
        DO 20 J=-JC+1,JC-1
          ITEMP=I**2+J**2
          IF ((ITEMP.GE.MAXR12).AND.(ITEMP.LT.MAXR22)) THEN
            IBIN=INT(SQRT(FLOAT(ITEMP))+0.5)+1
            JJ=J+1
            IF (JJ.LT.1) JJ=JJ+NSAM
              ID=II+JC*(JJ-1)
              Z=SPEC(ID)
            RBUF(IBIN)=RBUF(IBIN)+CABS(Z)**2
            RBUF(IBIN+JC)=RBUF(IBIN+JC)+1.0
          ENDIF
20    CONTINUE
C
      DO 60 I=1,JC
        IF (RBUF(I+JC).NE.0.0) RBUF(I)=SQRT(RBUF(I)/RBUF(I+JC))
        IF (RBUF(I).EQ.0.0) RBUF(I)=1.0
60    CONTINUE
C RBUF(I) = average amplitude of particle image at resolution I
C
      DO 30 I=0,JC-1
      	II=I+1
      	DO 30 J=-JC+1,JC-1
      	  ITEMP=I**2+J**2
          IF ((ITEMP.GE.MAXR12).AND.(ITEMP.LT.MAXR22)) THEN
            IBIN=INT(SQRT(FLOAT(ITEMP))+0.5)+1
      	    JJ=J+1
      	    IF (JJ.LT.1) JJ=JJ+NSAM
      	      ID=II+JC*(JJ-1)
              IF (ICTF.EQ.1) THEN
       	        CTFR=CTFF(ID)
      	        CTFL=CTFF(ID+IC)
                SNR=SSNR(IBIN)*CABS(CTFR+CONJG(CTFL))**2
              ELSE
                SNR=SSNR(IBIN)
              ENDIF
              WGT=SNR/(SNR+1.0)/RBUF(IBIN)
              Z=SPEC(ID)*WGT
              OUTC(ID)=Z
C              IF (ICTF.EQ.1) THEN
C                SNR=SSNR(IBIN)*CABS(CTFR+CONJG(CTFL))**2
C              ELSE
C                SNR=SSNR(IBIN)
C              ENDIF
C              WGT=SNR/(SNR+1.0)
C              PWEIGHTS(IBIN)=PWEIGHTS(IBIN)+WGT
      	  ELSEIF (ITEMP.GE.MAXR22) THEN
            IBIN=INT(SQRT(FLOAT(ITEMP))+0.5)+1
      	    JJ=J+1
      	    IF (JJ.LT.1) JJ=JJ+NSAM
      	      ID=II+JC*(JJ-1)
              OUTC(ID)=CMPLX(0.0,0.0)
          ENDIF
30    CONTINUE
C OUT(I) = weighted particle image amplitudes (w_i*X_i/X_n,k)
C
C      DO 40 I=1,JC
C        IF (RBUF(I+JC).NE.0.0) PWEIGHTS(I)=PWEIGHTS(I)/RBUF(I+JC)
C40    CONTINUE
C
      RETURN
      END
