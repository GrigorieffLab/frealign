C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE PINSERT_TRILIN_SE(NSAM,RII,IRREC2,SPEC,
     +		   CTFF,OCC,B3DF,T3DF,
     +		   TTD,SHX,SHY,ASUM,VSUM,PSUM,C3DF,KSUM,
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
      PARAMETER (IRAD=2)
      REAL X3,Y3,Z3,DM(9),SHX,SHY
      REAL CTFV2,DX,DY,CTFV2S2,DM1,DM2,DM3
      REAL PHASE,ST,ARG(3),BOXFTV,PX,RI2,RI12,OCC
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
      COMPLEX FOR,FOL,CTFR,CTFL,BSAMP
C**************************************************************************
        I1=-JC+1+NINT((IP-1)*ST)
        I3=-JC+NINT(IP*ST)
        I2=(I1+I3)/2
        I4=JC*(IP-1)+1
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
C
           DO 50 BEAM=1,2                               ! insert both left and right beam if IEWALD
            IF ((IEWALD.EQ.0).AND.(BEAM.EQ.2)) GOTO 30  ! skip 2nd beam if planar insertion
            IF (IEWALD.EQ.0) THEN
      	      X3=DM1+DM(4)*J                        ! planar insertion
      	      Y3=DM2+DM(5)*J
      	      Z3=DM3+DM(6)*J
            ENDIF
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
                IF (IEWALD.LT.0) THEN
                  CTFR=CONJG(CTFR)
                  CTFL=CONJG(CTFL)
                ENDIF
                SAMP=SPEC(ID)*PSHFT
                IF (IEWALD.NE.0) THEN 
                  IF (BEAM.EQ.1) THEN
                    CALL EWALDIN (NSAM,IRADA,C3DF,SINCLUT,IPAD, ! reference-derived 
     +                    I,J,DM,THET,XYZ,FOR,FOL,SAMP,PREL,IEWALD,  ! structure factors 
     +                    CTFR,CTFL)
                    X3=XYZ(1)
                    Y3=XYZ(2)
                    Z3=XYZ(3)
                    SAMP=FOR                                         ! right beam
                  ELSE
                    X3=XYZ(4)
                    Y3=XYZ(5)
                    Z3=XYZ(6)
                    SAMP=FOL                                         ! left beam
                  ENDIF
                ELSE
                  SAMP=SAMP*(CTFR+CONJG(CTFL))
                ENDIF
            ELSE
              II=-I+1
              JJ=-J+1
              IF (JJ.LT.1) JJ=JJ+NSAM
                ID=II+JC*(JJ-1)
       		CTFR=CTFF(ID)
       		CTFL=CTFF(ID+IC)
                IF (IEWALD.LT.0) THEN
                  CTFR=CONJG(CTFR)
                  CTFL=CONJG(CTFL)
                ENDIF
                SAMP=CONJG(SPEC(ID))*PSHFT
                IF (IEWALD.NE.0) THEN 
                  IF (BEAM.EQ.1) THEN
                    CALL EWALDIN (NSAM,IRADA,C3DF,SINCLUT,IPAD, ! reference-derived 
     +                    I,J,DM,THET,XYZ,FOR,FOL,SAMP,PREL,IEWALD,  ! structure factors 
     +                    CTFR,CTFL)
                    X3=XYZ(1)
                    Y3=XYZ(2)
                    Z3=XYZ(3)
                    SAMP=FOR                                         ! right beam
                  ELSE
                    X3=XYZ(4)
                    Y3=XYZ(5)
                    Z3=XYZ(6)
                    SAMP=FOL                                         ! left beam
                  ENDIF
                ELSE
                  SAMP=SAMP*(CTFR+CONJG(CTFL))
                ENDIF
      	    ENDIF
      	    CTFV2=0.25*CABS(CTFR+CONJG(CTFL))**2*WGT
            SAMP=SAMP*WGT
C
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
              ARG(3)=ABS(Z3-N)
      	      NN=N
      	      IF (NN.LT.0) NN=NN+NSAM
              N1=NSAM*NN
      	      DO 40 M=MS,MT
C      	        M2=M**2
                ARG(2)=ABS(Y3-M)
      		MM=M
      		IF (MM.LT.0) MM=MM+NSAM
                M1=JC*(MM+N1)
C                JNM=N2+M2
      	        DO 40 L=LS,LT
C      	          L2=L**2
                  ARG(1)=ABS(X3-L)
C      		  boxftv=BOXFT_LUT(arg,SINCLUT)
                  boxftv=(1.0-ARG(1))*(1.0-ARG(2))*(1.0-ARG(3))
C      		  	JRAD2=JNM+L2
C      		  	IF (JRAD2.LE.IRREC2) THEN
C      	          		LL=L
      		  		ctfv2s2=ctfv2*boxftv
      		  		BSAMP=samp*boxftv
      		  		    ID=L+M1+1
      		  		    B3DF(ID)=B3DF(ID)+BSAMP
      		  		    T3DF(ID)=T3DF(ID)+CTFV2S2
                                    NKSUM=NKSUM+1
                                    IF (NNSTAT.NE.0) THEN
      		  		      ID4=IPAD*L+1+IPAD*JC
     .                                  *(MM+IPAD*NSAM*NN)
                                      KSUM(ID)=KSUM(ID)+1
                                      VSUM(ID)=VSUM(ID)+CABS(BSAMP)**2
				      PSUM(ID)=PSUM(ID)
     .                                       +PDIFF(BSAMP,C3DF(ID4))
                                      ASUM(ID)=ASUM(ID)+CABS(BSAMP)
                                    ENDIF
C                                    IBIN=INT(SQRT(REAL(JRAD2))+0.5)+I4
C                                    NS1(IBIN)=NS1(IBIN)+1
      		  		  TTD=TTD+CTFV2S2
40	    CONTINUE
C
50         CONTINUE
           ENDIF
30	CONTINUE
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
C
           DO 51 BEAM=1,2                               ! insert both left and right beam if IEWALD
            IF ((IEWALD.EQ.0).AND.(BEAM.EQ.2)) GOTO 31  ! skip 2nd beam if planar insertion
            IF (IEWALD.EQ.0) THEN
      	      X3=DM1+DM(4)*J                        ! planar insertion
      	      Y3=DM2+DM(5)*J
      	      Z3=DM3+DM(6)*J
            ENDIF
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
                IF (IEWALD.LT.0) THEN
                  CTFR=CONJG(CTFR)
                  CTFL=CONJG(CTFL)
                ENDIF
                SAMP=SPEC(ID)*PSHFT
                IF (IEWALD.NE.0) THEN 
                  IF (BEAM.EQ.1) THEN
                    CALL EWALDIN (NSAM,IRADA,C3DF,SINCLUT,IPAD, ! reference-derived 
     +                    I,J,DM,THET,XYZ,FOR,FOL,SAMP,PREL,IEWALD,  ! structure factors 
     +                    CTFR,CTFL)
                    X3=XYZ(1)
                    Y3=XYZ(2)
                    Z3=XYZ(3)
                    SAMP=FOR                                         ! right beam
                  ELSE
                    X3=XYZ(4)
                    Y3=XYZ(5)
                    Z3=XYZ(6)
                    SAMP=FOL                                         ! left beam
                  ENDIF
                ELSE
                  SAMP=SAMP*(CTFR+CONJG(CTFL))
                ENDIF
            ELSE
              II=-I+1
              JJ=-J+1
              IF (JJ.LT.1) JJ=JJ+NSAM
                ID=II+JC*(JJ-1)
       		CTFR=CTFF(ID)
       		CTFL=CTFF(ID+IC)
                IF (IEWALD.LT.0) THEN
                  CTFR=CONJG(CTFR)
                  CTFL=CONJG(CTFL)
                ENDIF
                SAMP=CONJG(SPEC(ID))*PSHFT
                IF (IEWALD.NE.0) THEN 
                  IF (BEAM.EQ.1) THEN
                    CALL EWALDIN (NSAM,IRADA,C3DF,SINCLUT,IPAD, ! reference-derived 
     +                    I,J,DM,THET,XYZ,FOR,FOL,SAMP,PREL,IEWALD,  ! structure factors 
     +                    CTFR,CTFL)
                    X3=XYZ(1)
                    Y3=XYZ(2)
                    Z3=XYZ(3)
                    SAMP=FOR                                         ! right beam
                  ELSE
                    X3=XYZ(4)
                    Y3=XYZ(5)
                    Z3=XYZ(6)
                    SAMP=FOL                                         ! left beam
                  ENDIF
                ELSE
                  SAMP=SAMP*(CTFR+CONJG(CTFL))
                ENDIF
      	    ENDIF
      	    CTFV2=0.25*CABS(CTFR+CONJG(CTFL))**2*WGT
            SAMP=SAMP*WGT
C
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
              ARG(3)=ABS(Z3-N)
      	      NN=N
      	      IF (NN.LT.0) NN=NN+NSAM
              N1=NSAM*NN
      	      DO 41 M=MS,MT
C      	        M2=M**2
                ARG(2)=ABS(Y3-M)
      		MM=M
      		IF (MM.LT.0) MM=MM+NSAM
                M1=JC*(MM+N1)
C                JNM=N2+M2
      	        DO 41 L=LS,LT
C      	          L2=L**2
                  ARG(1)=ABS(X3-L)
C      		  boxftv=BOXFT_LUT(arg,SINCLUT)
                  boxftv=(1.0-ARG(1))*(1.0-ARG(2))*(1.0-ARG(3))
C      		  	JRAD2=JNM+L2
C      		  	IF (JRAD2.LE.IRREC2) THEN
C      	          		LL=L
      		  		ctfv2s2=ctfv2*boxftv
      		  		BSAMP=samp*boxftv
      		  		    ID=L+M1+1
      		  		    B3DF(ID)=B3DF(ID)+BSAMP
      		  		    T3DF(ID)=T3DF(ID)+CTFV2S2
                                    NKSUM=NKSUM+1
                                    IF (NNSTAT.NE.0) THEN
      		  		      ID4=IPAD*L+1+IPAD*JC
     .                                  *(MM+IPAD*NSAM*NN)
                                      KSUM(ID)=KSUM(ID)+1
                                      VSUM(ID)=VSUM(ID)+CABS(BSAMP)**2
				      PSUM(ID)=PSUM(ID)
     .                                       +PDIFF(BSAMP,C3DF(ID4))
                                      ASUM(ID)=ASUM(ID)+CABS(BSAMP)
                                    ENDIF
C                                    IBIN=INT(SQRT(REAL(JRAD2))+0.5)+I4
C                                    NS1(IBIN)=NS1(IBIN)+1
      		  		  TTD=TTD+CTFV2S2
41	    CONTINUE
51         CONTINUE
C
           ENDIF
31	CONTINUE
C
      	RETURN
      	END

