C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE BEAUTIFY(NSAM,NSYM,SYMOP,JSYM,A3DV,B3DV,RI,
     +                    HALFW,IMP)
C**************************************************************************
C Applies (again) the full symmetry, this time to the real space map.  
C It  also applies a 3D cosine bell mask, to beautify the final real space 
C structure before output.   May also help to improve the final phases 
C with a small amount of solvent flattening.
C Calls MATMUL and uses function TRILINMAP.
C Used in Frealign.
C**************************************************************************
      IMPLICIT NONE
      EXTERNAL MATMUL
C
      INTEGER JSYM(*),NSAM,NSYM,IMP,I1,I2
      INTEGER I,J,K,ISYM,ID,IS
      REAL SYMOP(3,3,*),A3DV(*),B3DV(*)
      REAL XC,YC,ZC,TM(3,3),TOTSYM,ST
      REAL RI,RII,RI3,HALFW
C**************************************************************************
C     XC=1.0 + 0.5*(FLOAT(NSAM)-1.0) ! in earlier versions 4.06 and 3.041
      XC=1.0 + 0.5*FLOAT(NSAM)
      YC=XC
      ZC=XC
      TOTSYM=1.0
      write(*,5) XC,YC,ZC,NSAM
5     format(/' Entering BEAUTIFY with centre at',3F9.3,
     .        '  for map with NSAM =',I5)
      RII=RI+HALFW/2.0
      RII=MIN(RII,FLOAT(NSAM/2)-1.0)
      RI3=RII**2

      DO 11 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 11 I=1,NSAM
          ID=IS+I
          B3DV(ID)=A3DV(ID)
11    CONTINUE

      ST=REAL(NSAM)/REAL(IMP)
      DO 20 K=1,NSYM
      	TOTSYM=TOTSYM*JSYM(K)
      	DO 15 I=1,3
      	  DO 15 J=1,3
      	    TM(I,J)=0.0
15	CONTINUE
      	TM(1,1)=1.0
      	TM(2,2)=1.0
      	TM(3,3)=1.0
      	DO 10 ISYM=1,JSYM(K)-1
          CALL MATMUL(TM,SYMOP(1,1,K),TM)
!$OMP PARALLEL DO
          DO 30 I=1,IMP
            CALL BEAUTIFY_S(NSAM,A3DV,B3DV,RI3,TM,XC,YC,ZC,
     +                      ST,I)
30	  CONTINUE
10	CONTINUE
        DO 13 J=0,NSAM*NSAM-1
          IS=J*(NSAM+2)
          DO 13 I=1,NSAM
            ID=IS+I
      	    A3DV(ID)=B3DV(ID)
13      CONTINUE
20    CONTINUE
      DO 14 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 14 I=1,NSAM
          ID=IS+I
      	  A3DV(ID)=A3DV(ID)/TOTSYM
14    CONTINUE
      write(*,21) INT(TOTSYM)
21    format(' TOTSYM =',I6/)
C
51    CONTINUE
      RETURN
      END
