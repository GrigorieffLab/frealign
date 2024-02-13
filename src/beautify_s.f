C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE BEAUTIFY_S(NSAM,A3DV,B3DV,RI3,TM,XC,YC,ZC,
     +                      ST,IP)
C**************************************************************************
C Called by BEAUTIFY, performs main loop.
C Uses function TRILINMAP.
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER L,M,N,NSAM,ID,I1,I2,IP,I3
      REAL A3DV(*),B3DV(*),ST
      REAL X3,Y3,Z3,XC,YC,ZC,XL,YM,ZN
      REAL AINPO,TRILINMAP,TM(3,3)
      REAL RAD2,RI3
C**************************************************************************
      I1=1+(IP-1)*ST
      I3=IP*ST
      I2=(I1+I3)/2
C
      DO 10 L=I1,I2
      	XL=L-XC
      	DO 10 M=1,NSAM
      	  YM=M-YC
      	  DO 10 N=1,NSAM
      	    ZN=N-ZC
      	    RAD2=XL**2+YM**2+ZN**2
      	    IF(RAD2.LT.RI3) THEN
      	      X3=XC+TM(1,1)*XL+TM(1,2)*YM+TM(1,3)*ZN
      	      Y3=YC+TM(2,1)*XL+TM(2,2)*YM+TM(2,3)*ZN
      	      Z3=ZC+TM(3,1)*XL+TM(3,2)*YM+TM(3,3)*ZN
      	      AINPO=TRILINMAP(NSAM,A3DV,X3,Y3,Z3)
      	      ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
              B3DV(ID)=B3DV(ID)+AINPO
            ENDIF
10    CONTINUE
C
C$OMP BARRIER
C
      DO 11 L=I2+1,I3
      	XL=L-XC
      	DO 11 M=1,NSAM
      	  YM=M-YC
      	  DO 11 N=1,NSAM
      	    ZN=N-ZC
      	    RAD2=XL**2+YM**2+ZN**2
      	    IF(RAD2.LT.RI3) THEN
      	      X3=XC+TM(1,1)*XL+TM(1,2)*YM+TM(1,3)*ZN
      	      Y3=YC+TM(2,1)*XL+TM(2,2)*YM+TM(2,3)*ZN
      	      Z3=ZC+TM(3,1)*XL+TM(3,2)*YM+TM(3,3)*ZN
      	      AINPO=TRILINMAP(NSAM,A3DV,X3,Y3,Z3)
      	      ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
              B3DV(ID)=B3DV(ID)+AINPO
            ENDIF
11    CONTINUE
      RETURN
      END
