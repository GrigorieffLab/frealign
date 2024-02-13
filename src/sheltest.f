C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE SHELTEST(NSAM,ISTP,A3DF,D3DF,ASUM,PSUM,
     +		      	  KSUM,QCP,NSTP,PR,FSC,ACC,QF,QC,PF,
     +			  NS,NPIX,NPIXT,DFSC,IREDUN,NNSTAT,PSSNR,
     +                    S3DF,V3DF,FMASK)
C**************************************************************************
C  Calculates some useful statistics in resolution shells
C  Uses function PDIFF.
C  Used in Frealign.
C**************************************************************************
      IMPLICIT NONE
      INTEGER NSAM,ISTP,NSTP,NPIX(*),I,IS,IL,JC,NQP,KSUM(*),IREDUN
      INTEGER L,M,N,LL,MM,NN,IS2,IL2,NPIXT(*),ID,IT,NS(*)
      INTEGER NNSTAT
      REAL PR(*),FSC(*),FRAD2,A1,A2,QCP(*),QC(*),PSSNR(*)
      REAL PDIFF,PI,PROD,ACC(*),SUMD,SUMS,ASUM(*),PSUM(*),QF(*),PF(*)
      REAL S3DF(*),V3DF(*),FMASK
      PARAMETER (PI=3.1415926535897)
      DOUBLEPRECISION MEAN12,MEAN22,MEANP,DFSC(*)
      COMPLEX A3DF(*),D3DF(*)
C**************************************************************************
      WRITE(*,1000)IREDUN
1000  FORMAT(/' SYMMETRY REDUNDANCY: ',I4/)
      JC=NSAM/2+1
C      NSAMH=NSAM/2
      NSTP=JC/ISTP
      IF (NSTP*ISTP.LT.JC) NSTP=NSTP+1
      DO 10 I=1,NSTP
      	NPIX(I)=0
      	NPIXT(I)=0
      	PR(I)=0.0
      	DFSC(I)=0.0
      	MEAN12=0.0D0
      	MEAN22=0.0D0
      	MEANP=0.0D0
      	FSC(I)=0.0
      	PSSNR(I)=0.0
      	ACC(I)=0.0
        IF (NNSTAT.NE.0) THEN 
     	  QF(I)=0.0
      	  QC(I)=0.0
      	  PF(I)=0.0
      	  NS(I)=0
        ENDIF
      	NQP=0
      	SUMD=0.0
      	SUMS=0.0
      	IS=(I-1)*ISTP
      	IL=I*ISTP
      	IF (I.EQ.NSTP) IL=JC
      	IL2=IL**2
      	IS2=IS**2
        DO 20 L=1,JC
      	  LL=L-1
      	  DO 20 M=1,NSAM
      	    MM=M-1
      	    IF (MM.GE.JC) MM=MM-NSAM
      	    DO 20 N=1,NSAM
      	      NN=N-1
      	      IF (NN.GE.JC) NN=NN-NSAM
      	      FRAD2=LL**2+MM**2+NN**2
      	      IF ((FRAD2.GE.IS2).AND.(FRAD2.LT.IL2)) THEN
      		  ID=L+JC*((M-1)+NSAM*(N-1))
      		  A1=CABS(A3DF(ID))
      		  A2=CABS(D3DF(ID))
      		  PROD=A1*A2
      		  IF (PROD.NE.0.0) THEN
      		    DFSC(I)=DFSC(I)+REAL(A3DF(ID)*CONJG(D3DF(ID)))
      		    MEANP=MEANP+DBLE(A1)*DBLE(A2)
      		    MEAN12=MEAN12+DBLE(A1)**2
      		    MEAN22=MEAN22+DBLE(A2)**2
      		    PR(I)=PR(I)+PDIFF(A3DF(ID),D3DF(ID))*A1*A2
      		    SUMD=SUMD+ABS(A1-A2)
      		    SUMS=SUMS+A1+A2
                    IF (NNSTAT.NE.0) THEN
      		      QF(I)=QF(I)+ASUM(ID)
      		      QC(I)=QC(I)+QCP(ID)
      		      PF(I)=PF(I)+PSUM(ID)
      		      NS(I)=NS(I)+KSUM(ID)
      		      NQP=NQP+1
                    ENDIF
      		    NPIX(I)=NPIX(I)+1
C Caution: PSSNR will be affected by the weighting function for
C individual particles. Needs to add some normalization.
                    PSSNR(I)=PSSNR(I)+S3DF(ID)+V3DF(ID)
      		  ENDIF
      		NPIXT(I)=NPIXT(I)+1
      	      ENDIF
20	CONTINUE
      	IF (NPIX(I).NE.0) THEN
C
C THIS FOURIER SHELL CORRELATION COEFFICIENT WAS DESCRIBED IN 
C VAN HEEL (1987), ULTRAMICROSCOPY 21, 95-100
          DFSC(I)=DFSC(I)/SQRT(MEAN12*MEAN22)
      	  FSC(I)=DFSC(I)
          IF (ABS(FSC(I)).LT.1.0) THEN
            A1=2.0*ABS(DFSC(I))/(1.0-ABS(DFSC(I)))
          ELSE
            A1=1000.0
          ENDIF
          IF (PSSNR(I).NE.0.0)
     +      PSSNR(I)=A1*NPIX(I)/PSSNR(I)*FMASK
c     +      PSSNR(I)=A1*NPIX(I)/PSSNR(I)
c     +              *SQRT(REAL(IREDUN))
      	  PR(I)=PR(I)/SNGL(MEANP)
      	  IF (SUMS.NE.0.0) ACC(I)=SUMD/SUMS*2.0
      	ENDIF
      	IF (NQP.NE.0) THEN
      	  QF(I)=QF(I)/NQP
      	  QC(I)=QC(I)/NQP
      	  PF(I)=PF(I)/NQP
      	  NS(I)=NS(I)/NQP
      	ENDIF
10    CONTINUE
      RETURN
      END
