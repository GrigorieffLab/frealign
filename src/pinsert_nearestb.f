C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PINSERT_NEARESTB(NSAM,RII,IRREC2,SPEC,
     +		   CTFF,OCC,B3DF,T3DF,TTD,
     +		   SHX,SHY,ASUM,VSUM,PSUM,C3DF,KSUM,
     +             SINCLUT,IRADA,IPAD,THET,IEWALD,
     +             NNSTAT,NKSUM,JC,DM,TMP,NSAMH,IC,DX,DY,
     +             ST,IP,NS1)
C**************************************************************************
C Called by PINSERT, performs main loop
C Uses function PDIFF
C shared: NKSUM,TTD
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,JC,I,J,II,JJ,LS,LT,MS,MT,NS,NT
      INTEGER L,M,N,LL,MM,NN,NSAMH,ID,IS,IRREC2
      INTEGER JRAD2,L2,M2,N2,KSUM(*),ID4,IPAD,IRADA,IC
      INTEGER NNSTAT,IRAD,I1,I2,IP,I3,NS1(*),IBIN,I4
      INTEGER IQ,JNM,N1,M1
      INTEGER*8 NKSUM
      PARAMETER (IRAD=1)
      REAL X3,Y3,Z3,DM(9),SHX,SHY
      REAL CTFV2,DX,DY,DM1,DM2,DM3
      REAL PHASE,ST,PX,RI2,RI12,OCC
      REAL RII,DELTA,DX2,DY2,DZ2
      REAL T3DF(*),WGT,PI,TMP,G2,FRAD2
      REAL VSUM(*),PSUM(*),ASUM(*),PDIFF,A
      REAL SINCLUT(*)
      DOUBLEPRECISION TTD
      PARAMETER (PI=3.1415926535897)
      COMPLEX SAMP,SPEC(*),CTFF(*)
      COMPLEX B3DF(*),C3DF(*),PSHFT
      INTEGER IEWALD
      REAL THET,XYZ(6),PREL
      INTEGER BEAM
      COMPLEX FOR,FOL,CTFR,CTFL
C**************************************************************************
        I1=-JC+1+NINT((IP-1)*ST)
        I3=-JC+NINT(IP*ST)
        I2=(I1+I3)/2
C        I4=JC*(IP-1)+1
        RI2=REAL(IRAD)/2.0
        RI12=REAL(IRAD-1)/2.0
C
      	DO 30 I=I1,I2
          DM1=DM(1)*I
          DM2=DM(2)*I
          DM3=DM(3)*I
          IQ=I**2
          PX=(SHX+DX)*I
      	  DO 30 J=-JC+1,JC-1
           G2=IQ+J**2
           IF (G2.LE.IRREC2) THEN
      	      X3=DM1+DM(4)*J
      	      Y3=DM2+DM(5)*J
      	      Z3=DM3+DM(6)*J
      	    WGT=EXP(TMP*G2)*OCC
            PHASE=PX+(SHY+DY)*J
            PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
      	    IF (I.GE.0) THEN
      	      II=I+1
      	      JJ=J+1
      	      IF (JJ.LT.1) JJ=JJ+NSAM
      		ID=II+JC*(JJ-1)
      		CTFR=CTFF(ID)
      		CTFL=CTFF(ID+IC)
                SAMP=SPEC(ID)*PSHFT
                  SAMP=SAMP*(CTFR+CONJG(CTFL))
            ELSE
              II=-I+1
              JJ=-J+1
              IF (JJ.LT.1) JJ=JJ+NSAM
                ID=II+JC*(JJ-1)
       		CTFR=CTFF(ID)
       		CTFL=CTFF(ID+IC)
                SAMP=CONJG(SPEC(ID))*PSHFT
                  SAMP=SAMP*(CTFR+CONJG(CTFL))
      	    ENDIF
      	    CTFV2=0.25*CABS(CTFR+CONJG(CTFL))**2*WGT
            SAMP=SAMP*WGT

            A=X3+0.5-RI12
      	    LS=INT(A)
      	    IF (A.LT.0.0) LS=LS-1
      	    IF (LS.LT.0) LS=0
            A=X3+RI2
      	    LT=INT(A)
      	    IF (A.LT.0.0) LT=LT-1
      	    IF (LT.GT.NSAMH) LT=NSAMH
            A=Y3+0.5-RI12
      	    MS=INT(A)
      	    IF (A.LT.0.0) MS=MS-1
      	    IF (MS.LT.-NSAMH) MS=-NSAMH
            A=Y3+RI2
      	    MT=INT(A)
      	    IF (A.LT.0.0) MT=MT-1
      	    IF (MT.GT.NSAMH) MT=NSAMH
            A=Z3+0.5-RI12
      	    NS=INT(A)
      	    IF (A.LT.0.0) NS=NS-1
      	    IF (NS.LT.-NSAMH) NS=-NSAMH
            A=Z3+RI2
      	    NT=INT(A)
      	    IF (A.LT.0.0) NT=NT-1
      	    IF (NT.GT.NSAMH) NT=NSAMH
      	    DO 40 N=NS,NT
C      	      N2=N**2
      	      NN=N
      	      IF (NN.LT.0) NN=NN+NSAM
              N1=NSAM*NN
      	      DO 40 M=MS,MT
C      	        M2=M**2
      		MM=M
      		IF (MM.LT.0) MM=MM+NSAM
                M1=JC*(MM+N1)
C                JNM=N2+M2
      	        DO 40 L=LS,LT
C      	          L2=L**2
C      		  	JRAD2=JNM+L2
C      		  	IF (JRAD2.LE.IRREC2) THEN
C      	          		LL=L
C      		  		ctfv2s2=ctfv2*boxftv
C      		  		BSAMP=samp*boxftv
      		  		    ID=L+M1+1
      		  		    B3DF(ID)=B3DF(ID)+SAMP
      		  		    T3DF(ID)=T3DF(ID)+CTFV2
                                    NKSUM=NKSUM+1
C                                    IBIN=INT(SQRT(REAL(JRAD2))+0.5)+I4
C                                    NS1(IBIN)=NS1(IBIN)+1
      		  		  TTD=TTD+CTFV2
40	    CONTINUE
           ENDIF
30	CONTINUE
C
C$OMP BARRIER
C
      	DO 31 I=I2+1,I3
          DM1=DM(1)*I
          DM2=DM(2)*I
          DM3=DM(3)*I
          IQ=I**2
          PX=(SHX+DX)*I
      	  DO 31 J=-JC+1,JC-1
           G2=IQ+J**2
           IF (G2.LE.IRREC2) THEN
      	      X3=DM1+DM(4)*J
      	      Y3=DM2+DM(5)*J
      	      Z3=DM3+DM(6)*J
      	    WGT=EXP(TMP*G2)*OCC
            PHASE=PX+(SHY+DY)*J
            PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
      	    IF (I.GE.0) THEN
      	      II=I+1
      	      JJ=J+1
      	      IF (JJ.LT.1) JJ=JJ+NSAM
      		ID=II+JC*(JJ-1)
      		CTFR=CTFF(ID)
      		CTFL=CTFF(ID+IC)
                SAMP=SPEC(ID)*PSHFT
                  SAMP=SAMP*(CTFR+CONJG(CTFL))
            ELSE
              II=-I+1
              JJ=-J+1
              IF (JJ.LT.1) JJ=JJ+NSAM
                ID=II+JC*(JJ-1)
       		CTFR=CTFF(ID)
       		CTFL=CTFF(ID+IC)
                SAMP=CONJG(SPEC(ID))*PSHFT
                  SAMP=SAMP*(CTFR+CONJG(CTFL))
      	    ENDIF
      	    CTFV2=0.25*CABS(CTFR+CONJG(CTFL))**2*WGT
            SAMP=SAMP*WGT

            A=X3+0.5-RI12
      	    LS=INT(A)
      	    IF (A.LT.0.0) LS=LS-1
      	    IF (LS.LT.0) LS=0
            A=X3+RI2
      	    LT=INT(A)
      	    IF (A.LT.0.0) LT=LT-1
      	    IF (LT.GT.NSAMH) LT=NSAMH
            A=Y3+0.5-RI12
      	    MS=INT(A)
      	    IF (A.LT.0.0) MS=MS-1
      	    IF (MS.LT.-NSAMH) MS=-NSAMH
            A=Y3+RI2
      	    MT=INT(A)
      	    IF (A.LT.0.0) MT=MT-1
      	    IF (MT.GT.NSAMH) MT=NSAMH
            A=Z3+0.5-RI12
      	    NS=INT(A)
      	    IF (A.LT.0.0) NS=NS-1
      	    IF (NS.LT.-NSAMH) NS=-NSAMH
            A=Z3+RI2
      	    NT=INT(A)
      	    IF (A.LT.0.0) NT=NT-1
      	    IF (NT.GT.NSAMH) NT=NSAMH
      	    DO 41 N=NS,NT
C      	      N2=N**2
      	      NN=N
      	      IF (NN.LT.0) NN=NN+NSAM
              N1=NSAM*NN
      	      DO 41 M=MS,MT
C      	        M2=M**2
      		MM=M
      		IF (MM.LT.0) MM=MM+NSAM
                M1=JC*(MM+N1)
C                JNM=N2+M2
      	        DO 41 L=LS,LT
C      	          L2=L**2
C      		  	JRAD2=JNM+L2
C      		  	IF (JRAD2.LE.IRREC2) THEN
C      	          		LL=L
C      		  		ctfv2s2=ctfv2*boxftv
C      		  		BSAMP=samp*boxftv
      		  		    ID=L+M1+1
      		  		    B3DF(ID)=B3DF(ID)+SAMP
      		  		    T3DF(ID)=T3DF(ID)+CTFV2
                                    NKSUM=NKSUM+1
C                                    IBIN=INT(SQRT(REAL(JRAD2))+0.5)+I4
C                                    NS1(IBIN)=NS1(IBIN)+1
      		  		  TTD=TTD+CTFV2
41	    CONTINUE
           ENDIF
31	CONTINUE
C
      	RETURN
      	END

