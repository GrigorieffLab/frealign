C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE MASKENV(NSAM,RI,OUTD,B3DV,SHX,SHY,PHI,THETA,PSI,W,
     +			 IBUF)
C**************************************************************************
C If XSTD is negative (see CALL), masks the extracted matching projection 
C  of the 3D structure using the D2MASK binary mask, and a cosine bell 
C  shaped ring at radius RI.

C Used in MATCH, CTFAPLLY and CTFAPPLY_PHASE_ONLY. 
C**************************************************************************
      IMPLICIT NONE

      INTEGER NSAM,I,J,K,ID,IXP,IYP,IXPP1,IYPP1,IRAD,M,N,NSAM2
      INTEGER IBUF,JJ,KK
      REAL SHX,SHY,PHI,THETA,PSI,OUTD(*),B3DV(*),RI,PSHX,PSHY,W(*)
      REAL CPHI,SPHI,CTHE,STHE,CPSI,SPSI,XH,YH,ZH,RI2,RAD2,PI
      REAL DDX,DDY,DDZ,CX1,CX2,CY1,CY2,CZ1,CZ2,DMI(9),X2,Y2,EDGE
      PARAMETER  (PI=3.1415926535897)
      PARAMETER  (IRAD=15)
C**************************************************************************
      IF (IBUF.LT.0) GOTO 100
      CPHI=COS(-PSI)
      SPHI=SIN(-PSI)
      CTHE=COS(-THETA)
      STHE=SIN(-THETA)
      CPSI=COS(-PHI)
      SPSI=SIN(-PHI)
      XH=REAL(NSAM/2)+1.0
      YH=REAL(NSAM/2)+1.0
      ZH=REAL(NSAM/2)+1.0
      RI2=RI**2
      PSHX=-SHX*NSAM/PI/2.0
      PSHY=-SHY*NSAM/PI/2.0
      NSAM2=NSAM*(NSAM+2)

      DMI(1)=CPHI*CTHE*CPSI-SPHI*SPSI
      DMI(2)=SPHI*CTHE*CPSI+CPHI*SPSI
C      DMI(3)=-STHE*CPSI
      DMI(4)=-CPHI*CTHE*SPSI-SPHI*CPSI
      DMI(5)=-SPHI*CTHE*SPSI+CPHI*CPSI
C      DMI(6)=STHE*SPSI
      DMI(7)=STHE*CPHI
      DMI(8)=STHE*SPHI
C      DMI(9)=CTHE

      DO 9 I=1,NSAM*NSAM
      	W(I)=0.0
9     CONTINUE

      DO 10 K=1,NSAM
        DDZ=(K-ZH)**2
        KK=NSAM*(K-1)
        DO 10 J=1,NSAM
          DDY=(J-YH)**2
          JJ=(NSAM+2)*(J-1+KK)
          DO 10 I=1,NSAM
            DDX=(I-XH)**2
            RAD2=DDX+DDY+DDZ
            IF (RAD2.LE.RI2) THEN
              ID=I+JJ
      	      IF (B3DV(ID).NE.0.0) THEN
                CX1=DMI(1)*(I-XH)
                CX2=DMI(2)*(I-XH)
                CY1=DMI(4)*(J-YH)
                CY2=DMI(5)*(J-YH)
                CZ1=DMI(7)*(K-ZH)
                CZ2=DMI(8)*(K-ZH)
                X2=CX1+CY1+CZ1
                Y2=CX2+CY2+CZ2
                X2=X2+XH+PSHX
                Y2=Y2+YH+PSHY
                IXP=INT(X2)
                IYP=INT(Y2)
                IXPP1=IXP+1
                IYPP1=IYP+1
                IF ((IXP .GE. 1).AND.(IYP .GE. 1).AND.
     +             (IXPP1 .LE. NSAM).AND.(IYPP1 .LE. NSAM)) THEN

                ID=IXP+(NSAM+2)*(IYP-1)
      		W(ID)=1.0
                ID=IXP+(NSAM+2)*(IYPP1-1)
      		W(ID)=1.0
                ID=IXPP1+(NSAM+2)*(IYP-1)
      		W(ID)=1.0
                ID=IXPP1+(NSAM+2)*(IYPP1-1)
      		W(ID)=1.0

                ENDIF
              ENDIF
      	    ENDIF
10    CONTINUE

      DO 11 J=1,NSAM
        JJ=(NSAM+2)*(J-1)
        DO 11 I=1,NSAM
          ID=I+JJ
      	  IF (W(ID).EQ.1.0) THEN
      	    DO 12 M=-IRAD,IRAD
      	      DO 12 N=-IRAD,IRAD
      	        RAD2=M**2+N**2
      	        RAD2=SQRT(RAD2)
                EDGE=(1.0+COS(PI*RAD2/IRAD))/2.0
                ID=I+M+(NSAM+2)*(J+N-1)
      	        IF ((ID.GE.1).AND.(ID.LE.NSAM2).AND.(RAD2.LE.IRAD))THEN
      	          IF (W(ID).LT.EDGE) W(ID)=EDGE
      	        ENDIF
12	    CONTINUE
      	  ENDIF
11    CONTINUE

100   CONTINUE
      DO 20 I=1,NSAM
        ID=(NSAM+2)*(I-1)
        DO 20 J=1,NSAM
          JJ=J+ID
          OUTD(JJ)=OUTD(JJ)*W(JJ)
20    CONTINUE
      RETURN
      END
