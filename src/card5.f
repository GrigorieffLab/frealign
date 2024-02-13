C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      SUBROUTINE CARD5(ASYMTEST,ASYM,JSYM,IREDUN,SYMOP,ASYMSTORE,
     +        NSYM,NASYM,ISYMAX,ALPHA,RISE,NU,HSTART,STIF,NONSYM)
C**************************************************************************
C Read input card 5 (symmetry) global parameters
C Card 5:  ASYM - symmetry required Cn,Dn,T,O,I,I1,I2,N (can be zero)
C                                             or H (helical symmetry)
C Card 5a: only if N and N.ne.0, ((SYMOP(J,K,I), J=1,3), K=1,3), I=1,N
C          or H
C Calls GETSYMMAT and CHECKSYM.
C Called by Frealign.
C*************************************************************************
        IMPLICIT NONE
        CHARACTER*16 ASYMTEST
        CHARACTER*3 ASYM
        CHARACTER*5 UPC,LOC
        INTEGER ISYMAX,NSYM,NASYM,ISYM,J,I
        INTEGER KSYMN1,KSYMN2,KSYMN3,NU
        INTEGER JSYM(*),IREDUN,HSTART
        REAL SYMOP(3,3,*),ALPHA,RISE,STIF
        DATA UPC/'CDTOI'/
        DATA LOC/'cdtoi'/
        LOGICAL ASYMSTORE,NONSYM
C*************************************************************************
C  Read in or retrieve rotation symmetry matrices.
      	WRITE(*,*)' SYMMETRY CARD ?'
      	READ(*,7015) ASYM
      	WRITE(*,17015) ASYM
7015	FORMAT(A3)
17015	FORMAT(3X,A3)
C
        NONSYM=.FALSE.
        DO 10 I=1,3
          DO 10 J=1,5
            IF (ASYM(I:I).EQ.LOC(J:J)) THEN
              NONSYM=.TRUE.
              ASYM(I:I)=UPC(J:J)
            ENDIF
10      CONTINUE
        IF (NONSYM) WRITE(*,6005)
6005    FORMAT(/' ***WARNING*** Performing',
     +    ' nonsymmetric refinement...'/)
C
      	ASYMSTORE=.FALSE.
      	KSYMN1=-1
      	KSYMN2=-1
      	KSYMN3=-1
C
        IF (ASYM(1:1).EQ.'0') ASYM='C1 '
C
        IF (ASYM(1:1).EQ.'H') THEN
          NSYM=0
          WRITE(*,6004)
6004      FORMAT(' HELICAL PARAMETER INPUT ?',
     + /' (Alpha [deg], Rise [A], # Subunits to average,',
     +     ' # Starts, Stiffness)')
          READ(*,*)ALPHA,RISE,NU,HSTART,STIF
          WRITE(*,6003)ALPHA,RISE,NU,HSTART,STIF
6003      FORMAT(2F10.5,2I10,F10.5)
          IF (NU.LT.1) NU=1
          IF (HSTART.LT.1) HSTART=1
          IF (RISE.LT.0.0) THEN
            RISE=-RISE
            ALPHA=-ALPHA
          ENDIF
          IREDUN=HSTART*NU
        ELSE
          STIF=1.0
C
      		IF(ASYM(1:1).EQ.' ') STOP ' Funny 0???' 	
C							! only letters and digits allowed
      	DO 8500 J=2,16
      		IF(ASYM(3:3).EQ.ASYMTEST(J:J).AND.J.LE.6) STOP 'Very funny 1??'
C								! no letters allowed
      		IF(ASYM(3:3).EQ.ASYMTEST(J:J).AND.J.GT.6) KSYMN3=J-7
      		IF(ASYM(2:2).EQ.ASYMTEST(J:J).AND.J.LE.6) STOP 'Very funny 2??'
C								! no letters allowed
      		IF(ASYM(2:2).EQ.ASYMTEST(J:J).AND.J.GT.6) KSYMN2=J-7
      		IF(ASYM(1:1).EQ.ASYMTEST(J:J)) THEN
      			KSYMN1=0
      			IF(J.LE.6) ASYMSTORE=.TRUE.
      			IF(J.GT.6) KSYMN1=J-7
      		ENDIF
8500	CONTINUE
      		IF(KSYMN1.EQ.-1) STOP ' Symmetry request invalid'
8501	IF(ASYMSTORE) THEN
                    IF(KSYMN3.NE.-1) THEN
                        NASYM=KSYMN3+10*KSYMN2
                    ELSE
                        IF(KSYMN2.NE.-1) THEN
                                NASYM=KSYMN2
                        ELSE
                                NASYM=0
                        ENDIF
                    ENDIF
      		CALL GETSYMMAT(ASYM(1:1),NASYM,SYMOP,NSYM,JSYM,ISYMAX,IREDUN)
      	ELSE
      		    IF(KSYMN3.NE.-1) THEN
      			NSYM=KSYMN3+10*KSYMN2+100*KSYMN1
      		    ELSE
      			IF(KSYMN2.NE.-1) THEN
      				NSYM=KSYMN2+10*KSYMN1
      			ELSE
      				NSYM=KSYMN1
      			ENDIF
      		    ENDIF
      		    IF(NSYM.GT.ISYMAX)STOP' Maximum ISYMAX SYMOPs allowed'
      		WRITE(*,*)' SYMMETRY MATRIX INPUT ?'
      		IREDUN=1
      		DO 800 I=1,NSYM
      			READ(*,*) (SYMOP(J,1,I), J=1,3)
      			WRITE(*,6000) (SYMOP(J,1,I), J=1,3)
      			READ(*,*) (SYMOP(J,2,I), J=1,3)
      			WRITE(*,6000) (SYMOP(J,2,I), J=1,3)
      			READ(*,*) (SYMOP(J,3,I), J=1,3)
      			WRITE(*,6000) (SYMOP(J,3,I), J=1,3)
6000			FORMAT(3F12.4)
      			CALL CHECKSYM(SYMOP(1,1,I),ISYM)
      			IF (ISYM.LT.0) THEN
      				WRITE(*,6001) I,ISYM
6001				FORMAT(' Symmetry operator',I4,' wrong: ISYM =',I4)
      				STOP
      			ELSE
      				WRITE(*,6002) I,ISYM
6002				FORMAT(' Symmetry operator',I4,'  ISYM =',I4)
      			ENDIF
      			JSYM(I)=ISYM
      			IREDUN=IREDUN*ISYM
800		CONTINUE
      	ENDIF
      	WRITE(*,*) ' Symmetry point group  ',ASYM,'  redundancy',IREDUN
C
        ENDIF
C
	RETURN
	END
