C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C ---------------------------------------------------------------------------------        
	SUBROUTINE STORESHIFT(STD,VTD,JC,NSAM,NSAMH,IREDUN,NN1,
     +			QCP,ASUM,KSUM,PSUM,VSUM,SSNR,A3DF,B3DF,
     +			C3DF,D3DF,S3DF,V3DF,REF3DF,IPAD,NNSTAT,
     +                  IFSC)
C ---------------------------------------------------------------------------------        
C       SHIFT VOLUME INTO CENTRE AFTER STORING THE REFERENCE TRANSFORM IN REF3DF/S
C       Used in Frealign.
C ---------------------------------------------------------------------------------        
        IMPLICIT NONE
        INTEGER NN1,L,LL,M,MM,N,NN,JC,NSAM,WC1,NSAMH,IS,ID
        INTEGER NUTD,ISUM,IREDUN,JRAD2,IRREC2,NSAM3
	INTEGER KSUM(*),ID4,IS4,IPAD,NNSTAT,IFSC,I
        REAL    PSHFTR,UTD,TMP,RLIM,PI,TMP2
	PARAMETER  (PI=3.1415926535897)
	DOUBLEPRECISION	STD,VTD
	REAL 	QCP(*),SSNR(*)
	REAL 	ASUM(*),VSUM(*),PSUM(*)
	REAL 	S3DF(*),V3DF(*)
        COMPLEX PSHFT,REF3DF(*)
	COMPLEX A3DF(*),B3DF(*)
        COMPLEX C3DF(*),D3DF(*)
C ---------------------------------------------------------------------------------        
       IF (IFSC.EQ.0) THEN
       IF(STD.EQ.0.0.OR.VTD.EQ.0.0) THEN
       WRITE(*,6102) STD,VTD
6102   FORMAT(' Abnormal termination, probably no acceptable particles'/
     . '   either STD or VTD zero, STD,VTD=',2F15.7/
     . '   at least two particles needed')
       STOP 'STOP 6102'
       ENDIF
       ENDIF

       NSAM3=NSAM*NSAM*NSAM
       TMP=(1.0*SQRT(4.0/PI-1.0)+1.0)*SQRT(PI)/2.0

        DO 60 N=1,NSAM
          NN=N-1
          IF (NN.GE.JC) NN=NN-NSAM
          DO 60 M=1,NSAM
            MM=M-1
            IF (MM.GE.JC) MM=MM-NSAM
            DO 60 L=1,JC
              LL=L-1
              ISUM=(LL+MM+NN)
              PSHFTR=1.0
              IF (MOD(ISUM,2).NE.0) PSHFTR=-1.0
              I=INT(SQRT(REAL(LL**2+MM**2+NN**2))+0.5)+1
                ID=L+JC*((M-1)+NSAM*(N-1))
            ID4=IPAD*(L-1)+1+IPAD*(IPAD*NSAMH+1)*((M-1)+IPAD*NSAM*(N-1))
                IF (IFSC.EQ.0) THEN
                  B3DF(ID)=A3DF(ID)*PSHFTR/(S3DF(ID)+STD*IREDUN)/NSAM3
                ENDIF
                IF (NNSTAT.NE.0) REF3DF(ID)=C3DF(ID4)*PSHFTR      ! STORE
                IF (IFSC.EQ.0) THEN
                  C3DF(ID)=D3DF(ID)*PSHFTR/(V3DF(ID)+VTD*IREDUN)/NSAM3
                ENDIF
                IF (NNSTAT.NE.0) THEN
                  IF (KSUM(ID).NE.0) THEN
                    IF(ASUM(ID).NE.0.)
     +                QCP(ID)=TMP*SQRT(VSUM(ID))/ASUM(ID)
                    IF (IFSC.EQ.0) THEN
                      TMP2=VSUM(ID)-CABS(A3DF(ID)+D3DF(ID))**2/KSUM(ID)
                      IF (TMP2.GT.0.0) THEN
                        SSNR(ID)=CABS(A3DF(ID)+D3DF(ID))**2/TMP2
     +                          *(KSUM(ID)-1)/KSUM(ID)
                      ELSE
                        SSNR(ID)=0.0
                      ENDIF
                    ELSE
                      TMP2=VSUM(ID)-CABS(A3DF(ID))**2/KSUM(ID)
                      IF (TMP2.GT.0.0) THEN
                        SSNR(ID)=CABS(A3DF(ID))**2/TMP2
     +                          *(KSUM(ID)-1)/KSUM(ID)
                      ELSE
                        SSNR(ID)=0.0
                      ENDIF
                    ENDIF
C                   IF (SSNR(ID).LT.1.0) SSNR(ID)=0.0  ! removed 12.1.99 (RH)
                    VSUM(ID)=VSUM(ID)-ASUM(ID)**2/KSUM(ID)
                    IF (KSUM(ID).NE.1) THEN
                      IF (VSUM(ID).GE.0.0) VSUM(ID)
     +                                  =SQRT(VSUM(ID)/(KSUM(ID)-1))
                    ELSE
                      IF (VSUM(ID).GE.0.0) VSUM(ID)=SQRT(VSUM(ID))
                    ENDIF
                    IF (IFSC.EQ.0) THEN
                      ASUM(ID)=CABS(A3DF(ID)+D3DF(ID))/ASUM(ID)
                    ELSE
                      ASUM(ID)=CABS(A3DF(ID))/ASUM(ID)
                    ENDIF
C                   PSUM(ID)=SQRT(PSUM(ID)/KSUM(ID))
                    PSUM(ID)=PSUM(ID)/KSUM(ID)
                  ENDIF
                ENDIF
60      CONTINUE
	RETURN
	END
