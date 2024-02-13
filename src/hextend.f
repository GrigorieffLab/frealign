C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE HEXTEND(NSAM,ALPHA,RR,A3DV,B3DV,
     +                   RI,NS,NMIN,NMAX,LSUM)
C**************************************************************************
C Extends helical symmetry to the real space map.  
C also applies a 3D cosine bell mask to the final map.
C Uses function TRILINMAP.
C Used in Frealix.
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER I,L,M,N,NSAM,ID,ICNT,NS,K,IS
      INTEGER K1,K2,NZR,N1,N2,NMIN,NMAX,J
      REAL A3DV(*),B3DV(*),ALPHAS,RS,RR
      REAL X3,Y3,Z3,XC,YC,XL,YM,LSUM(*)
      REAL SUM,AS,SINA,COSA,TOTSYM
      REAL AINPO,TRILINMAP,A,R,ALPHA
      REAL RAD2,RI,RI2,RI3,MEAN,PI
      PARAMETER (PI=3.1415926535897)
C**************************************************************************
      XC=1.0+INT(NSAM/2)
      YC=XC
      WRITE(*,5) XC,YC,NSAM,NMIN,NMAX
5     FORMAT(/' Entering HEXTEND with centre at',2F9.3,
     .        '  for map with NSAM =',I5,
     .       /' Extending region between ZMIN,ZMAX =',2I9)
      TOTSYM=0.0
C
      R=-RR
      RI2=(RI-1)**2
      RI3=RI**2
      IF (RI.GT.NSAM/2-1) THEN
        RI2=(NSAM/2-2)**2
        RI3=(NSAM/2-1)**2
      ENDIF
      K1=-NINT(NS/2.0)
      K2=NINT(NS/2.0)
      IF (K2-K1+1.GT.NS) K2=K2-1
      IF (K2-K1+1.GT.NS) K1=K1+1
C
      DO 20 I=1,NSAM*NSAM*(NSAM+2)
        B3DV(I)=0.0
20    CONTINUE
      DO 21 I=1,NSAM
        LSUM(I)=0.0
21    CONTINUE
C
      I=NINT(2.0*PI/ALPHA)
      RS=I*R/NS
      ALPHAS=I*ALPHA-2.0*PI 
C
      NZR=NSAM/R+1
      DO 10 K=K1,K2
        AS=K*ALPHAS
        DO 10 I=1-NSAM/R-NS,NZR+NSAM/R+NS
          A=I*ALPHA+K*ALPHAS
          SINA=SIN(A)
          COSA=COS(A)
          N1=MAX(1,INT(NMIN+I*R+K*RS))
          IF (N1-I*R-K*RS.LT.1.0) N1=N1+1
          N2=MIN(NSAM,INT(NMAX+I*R+K*RS))
          IF ((N1.LE.NSAM/2).AND.(N2.GE.NSAM/2))
     +      TOTSYM=TOTSYM+1.0
!$OMP PARALLEL DO
          DO 11 N=N1,N2
            CALL HEXTEND_S(NSAM,N,I,R,K,RS,XC,YC,
     +              RI3,COSA,SINA,A3DV,B3DV,LSUM)
11        CONTINUE
10    CONTINUE
C
      MEAN=0.0
      ICNT=0
!$OMP PARALLEL DO REDUCTION(+:MEAN,ICNT)
      DO 30 N=1,NSAM
        CALL SUMA(NSAM,N,MEAN,ICNT,XC,YC,
     +            RI2,RI3,B3DV,LSUM)
30    CONTINUE
      IF (ICNT.NE.0) MEAN=MEAN/ICNT
C
!$OMP PARALLEL DO
      DO 40 N=1,NSAM
        CALL SETM(NSAM,N,MEAN,XC,YC,
     +                RI3,B3DV)
40    CONTINUE
C
      DO 50 J=0,NSAM*NSAM-1
        IS=J*(NSAM+2)
        DO 50 I=1,NSAM
          ID=IS+I
          A3DV(ID)=B3DV(ID)/TOTSYM
50    CONTINUE
      WRITE(*,22) INT(TOTSYM)
22    FORMAT(' TOTSYM =',I6/)
C
      RETURN
      END
C
C**************************************************************************
      SUBROUTINE HEXTEND_S(NSAM,N,I,R,K,RS,XC,YC,
     +              RI3,COSA,SINA,A3DV,B3DV,LSUM)
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER I,L,M,N,NSAM,ID,K
      REAL A3DV(*),B3DV(*),RS,LSUM(*)
      REAL X3,Y3,Z3,XC,YC,XL,YM
      REAL AINPO,TRILINMAP,R,SINA,COSA
      REAL RAD2,RI3
C**************************************************************************
      LSUM(N)=LSUM(N)+1.0
      Z3=N-I*R-K*RS
      DO 12 M=1,NSAM
        YM=M-YC
        DO 12 L=1,NSAM
          XL=L-XC
          RAD2=XL**2+YM**2
          ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
          IF (RAD2.LT.RI3) THEN
            X3=XC+COSA*XL+SINA*YM
            Y3=YC-SINA*XL+COSA*YM
            AINPO=TRILINMAP(NSAM,A3DV,X3,Y3,Z3)
            B3DV(ID)=B3DV(ID)+AINPO
          ENDIF
12    CONTINUE
C
      RETURN
      END
C
C**************************************************************************
      SUBROUTINE SUMA(NSAM,N,MEAN,ICNT,XC,YC,
     +                RI2,RI3,B3DV,LSUM)
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER L,M,N,NSAM,ID,ICNT
      REAL B3DV(*),YC,YM,XC,XL,RAD2
      REAL MEAN,RI2,RI3,LSUM(*)
C**************************************************************************
      DO 30 M=1,NSAM
        YM=M-YC
        DO 30 L=1,NSAM
          XL=L-XC
          ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
          B3DV(ID)=B3DV(ID)/LSUM(N)
          RAD2=XL**2+YM**2
          IF ((RAD2.GE.RI2).AND.(RAD2.LE.RI3)) THEN
            MEAN=MEAN+B3DV(ID)
            ICNT=ICNT+1
          ENDIF
30    CONTINUE
C
      RETURN
      END
C
C**************************************************************************
      SUBROUTINE SETM(NSAM,N,MEAN,XC,YC,
     +                RI3,B3DV)
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER L,M,N,NSAM,ID
      REAL B3DV(*),YC,YM,XC,XL,RAD2
      REAL MEAN,RI3
C**************************************************************************
      DO 40 M=1,NSAM
        YM=M-YC
      	DO 40 L=1,NSAM
          XL=L-XC
          RAD2=XL**2+YM**2
          IF (RAD2.GE.RI3) THEN
            ID=L+(NSAM+2)*((M-1)+NSAM*(N-1))
            B3DV(ID)=MEAN
          ENDIF
40    CONTINUE
C
      RETURN
      END
