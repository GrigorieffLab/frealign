C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE MASK3D_C(NSAM,A3DV,RIO,RIC,HALFW,FMASK,
     +                    AVE,STD)
C**************************************************************************
C Applies cylindrical mask with outer radius RIO, inner radius RIC, and
C cosine edge.
C Used in Frealix.
C**************************************************************************
      IMPLICIT NONE

      INTEGER L,M,N,NSAM,ID
      INTEGER I,J,K,NBG,NA
      REAL A3DV(*),XC,YC,ZC,XL,YM,ZN,RM,RIO
      REAL FMASK,AVE,STD,RI4,RI5
      REAL RAD,RAD2,HALFW,RI,RII,RIH,RI2,EDGE,PI,RI3,RIC
      DOUBLE PRECISION BG,SUM1,SUM2,SUM3
      PARAMETER  (PI=3.1415926535897)
C**************************************************************************
      RI=RIO
      RI4=RIO**2
      RI5=RIC**2
C
      RM=MIN(RI,FLOAT(NSAM/2)-1.0-HALFW/2.0)
      RII=RM+HALFW/2.0
      RIH=RM-HALFW/2.0
      IF (RIH.LT.0.0) RIH=0.0
      RI2=RIH**2
      RI3=RII**2
      BG=0.0
      NBG=0
      XC=1.0 + 0.5*FLOAT(NSAM)
      YC=XC
      ZC=XC
      SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      NA=0
      DO 50 L=1,NSAM
      	XL=L-XC
      	DO 50 M=1,NSAM
      	  YM=M-YC
      	  DO 50 N=1,NSAM
      	    ZN=N-ZC
      	    ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
      	    RAD2=XL**2+YM**2
      	    IF ((RAD2.GT.RI2).AND.(RAD2.LT.RI3)) THEN
              BG=BG+A3DV(ID)
              NBG=NBG+1
            ENDIF
            IF ((RAD2.LE.RI4).AND.(RAD2.GE.RI5)) THEN
              SUM1=SUM1+A3DV(ID)
              SUM2=SUM2+A3DV(ID)**2
              NA=NA+1
            ENDIF
50    CONTINUE
      IF (NBG.NE.0) BG=BG/NBG
      IF (NA.NE.0) THEN
        SUM1=SUM1/NA
        SUM2=SUM2/NA
      ENDIF
      AVE=SUM1
      STD=SQRT(ABS(SUM2-SUM1**2))
C
      DO 51 L=1,NSAM
      	XL=L-XC
      	DO 51 M=1,NSAM
      	  YM=M-YC
      	  DO 51 N=1,NSAM
      	    ZN=N-ZC
      	    ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
      	    RAD2=XL**2+YM**2
            EDGE=1.0
      	    IF (RAD2.GT.RI2) THEN
      	      RAD=SQRT(RAD2)
      	      IF(RAD.GT.RIH.AND.RAD.LT.RII) THEN
      	        EDGE=(1.0+COS(PI*(RAD-RIH)/HALFW))/2.0
      	      ELSE
      	        EDGE=0.0
      	      ENDIF
      	      A3DV(ID)=A3DV(ID)*EDGE+BG*(1.0-EDGE)
            ENDIF
            SUM3=SUM3+EDGE**2
51    CONTINUE
C
      IF (RIC.GT.0.0) THEN
C
      RI=RIC
C
      RM=MAX(RI,HALFW/2.0)
      RII=RM+HALFW/2.0
      RIH=RM-HALFW/2.0
      RI2=RIH**2
      RI3=RII**2
      BG=0.0
      NBG=0
      DO 60 L=1,NSAM
      	XL=L-XC
      	DO 60 M=1,NSAM
      	  YM=M-YC
      	  DO 60 N=1,NSAM
      	    ZN=N-ZC
      	    ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
      	    RAD2=XL**2+YM**2
      	    IF ((RAD2.GT.RI2).AND.(RAD2.LT.RI3)) THEN
              BG=BG+A3DV(ID)
              NBG=NBG+1
            ENDIF
60    CONTINUE
      IF (NBG.NE.0) BG=BG/NBG
C
      DO 61 L=1,NSAM
      	XL=L-XC
      	DO 61 M=1,NSAM
      	  YM=M-YC
      	  DO 61 N=1,NSAM
      	    ZN=N-ZC
      	    ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
      	    RAD2=XL**2+YM**2
            EDGE=1.0
      	    IF (RAD2.LT.RI3) THEN
      	      RAD=SQRT(RAD2)
      	      IF(RAD.GT.RIH.AND.RAD.LT.RII) THEN
      	        EDGE=(1.0+COS(PI*(RII-RAD)/HALFW))/2.0
      	      ELSE
      	        EDGE=0.0
      	      ENDIF
      	      A3DV(ID)=A3DV(ID)*EDGE+BG*(1.0-EDGE)
            ENDIF
            IF (EDGE.NE.1.0) SUM3=SUM3-1.0+EDGE**2
61    CONTINUE
C
      ENDIF
C
      FMASK=SUM3/NSAM/NSAM/NSAM
      RETURN
      END
