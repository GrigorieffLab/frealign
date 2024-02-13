C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE GETSYMMAT(ASYM,NASYM,SYMOP,NSYM,JSYM,ISYMAX,IREDUN)
C**************************************************************************
C  Retrieves the number of rotational matrices required and the corresponding
C    symmetry operators for any requested point group symmetry from the
C    arrays NSYM_CDTOI, NSYMTYPE and SYMSTORE, and calculates the rotational
C    order in subroutine CHECKSYM.

C  Calls CHECKSYM.
C  Used in CARD5.
C**************************************************************************
      IMPLICIT NONE

      INTEGER ISYMAX
      REAL PI
      PARAMETER (PI=3.1415926535897)

      INTEGER J,K,L,NASYM,NSYM,NGET,N_CDTOI,IREDUN
      INTEGER JSYM(*),NSYM_CDTOI(6),NSYMTYPE(4,6)
      REAL ANGLE,SYMOP(3,3,*),SYMSTORE(3,3,8)
      CHARACTER*1 ASYM
      CHARACTER*16 ASYMTEST
      DATA  ASYMTEST/' CDTOI0123456789'/
      DATA  NSYM_CDTOI/1,2,3,3,4,4/
      DATA  NSYMTYPE/8,0,0,0,
     +                  8,5,0,0,
     +                  4,5,6,0,
     +                  3,7,4,0,
     +                  1,4,5,6,
     +                  2,4,5,6/
C      DATA  SYMSTORE/0.5,-0.309017,-0.809017,-0.309017,0.809017,-0.5,
C     +                                          0.809017,0.5,0.309017,
      DATA  SYMSTORE/0.309017,0.809017,-0.5,-0.809017,0.5,0.309017,
     +                                          0.5,0.309017,0.809017,
     2          0.809017,0.309017,0.5,0.309017,0.5,-0.809017,
     +                                          -0.5,0.809017,0.309017,
     3           0.,  1.,  0., -1.,  0.,  0.,  0.,  0.,  1.,
     4           0.,  0.,  1.,  1.,  0.,  0.,  0.,  1.,  0.,
     5           1.,  0.,  0.,  0., -1.,  0.,  0.,  0., -1.,
     6          -1.,  0.,  0.,  0.,  1.,  0.,  0.,  0., -1.,
     7           0.,  1.,  0.,  1.,  0.,  0.,  0.,  0., -1.,
     8           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1./
C**************************************************************************
C  SYMSTORE stores all the required symmetry operators in eight locations
C       1.  5-fold along 01t, where t=(1+SQRT(5))/2 - default Int Tables icosahedral 5-fold
C       2.  5-fold along 10t    alternative RAC icosahedral 5-fold
C       3.  4-fold along 001    used in 432 (O)
C       4.  3-fold along 111    used in 23,432 and 532 (T,O,I)
C       5.  2-fold along 100    used in D,T,I
C       6.  2-fold along 010    used in T,I
C       7.  2-fold along 110    used only in 432 (O)
C       8.  N-fold along 001    created below as required for C and D pointgroups
C**************************************************************************
      N_CDTOI = 0
      IREDUN = 1
      WRITE(*,*)' Entering GETSYMMAT with ASYM,NASYM  ',ASYM,NASYM
      DO 40 J = 2,5
        IF(ASYM.EQ.ASYMTEST(J:J)) N_CDTOI = J-1
40    CONTINUE
      IF(ASYM.EQ.ASYMTEST(6:6)) THEN
C               ! default Int Tables icosahedral
        N_CDTOI = 5

C               ! alternative RAC icosahedral
        IF(NASYM.EQ.2) N_CDTOI = 6
      ENDIF
C               ! calculate N-fold rotation matrix (C,D)
      IF(N_CDTOI.LE.2.AND.N_CDTOI.GE.1) THEN
        IF(NASYM.EQ.0) STOP ' pointgroups C and D must be followed by n'
        ANGLE = 2.0*PI/NASYM
        SYMSTORE(1,1,8) = COS(ANGLE)
        SYMSTORE(1,2,8) = SIN(ANGLE)
        SYMSTORE(2,1,8) =-SIN(ANGLE)
        SYMSTORE(2,2,8) = COS(ANGLE)
      ENDIF

      IF(N_CDTOI.EQ.0) STOP ' Invalid symmetry request'
      NSYM = NSYM_CDTOI(N_CDTOI)
      DO 50 J=1,NSYM
        NGET=NSYMTYPE(J,N_CDTOI)
        DO 55 K=1,3
        DO 55 L=1,3
55      SYMOP(K,L,J)=SYMSTORE(K,L,NGET)
        CALL CHECKSYM(SYMOP(1,1,J),JSYM(J))
        IREDUN=IREDUN*JSYM(J)
50    CONTINUE
      RETURN
      END
