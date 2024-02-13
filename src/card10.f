C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
        SUBROUTINE CARD10(NANG,IFLAG,FINPAR,NNPART,ILIST,PSI,THETA,
     +               PHI,SHX,SHY,FILM,DFMID1,DFMID2,ANGAST,OCC,
     +               PRESA,AMEANPHAS,NPHAS,NSAM,AMAGP,NSET,PSIZE,
     +		     MAXSET,RELMAG,NSS,IFIRST,ILAST,NNSET,NSETST,
     +               DSTEP,TNSET,FSCT,THRESH,CMIN,CMAX,XM,YM,SX,SY,
     +               NPB,ASYM,THETAM,STHETA,PSIM,SPSI,ALPHA,RISE,QM,
     +               ALGP,SIG,RI,PSSNR,FBOOST,IREDUN,NONSYM,IRAN,P)
C**************************************************************************
C Read input card 10 
C   - at least one of the following is required
C  Card 10a     FINPAR - INPUT PARAMETER FILE, required if mode key=0,1,2,3,4
C  Card 10b     NIN,ABSMAGPIN,IFILMIN,DFMID1IN,DFMID2IN,ANGASTIN,MORE if mode<0
C                 number of particles, magnification, film number, defocus
C                 parameters for this dataset.  If MORE='1', then more cards
C                 of this type follow, until MORE=0 terminates the parameter
C                 data.  The information on these cards is used to create a
C                 new parameter file.
C Called by Frealign.
C*************************************************************************
        IMPLICIT NONE
        REAL PI
	PARAMETER  (PI=3.1415926535897)
      	INTEGER NANG,NANG1,NANG2,IFLAG,NNPART,MAXSET
        INTEGER FILM(*),ILIST(*),TNSET(*),I,NPB,INPB,SLEN2
        INTEGER NNSET(*),NSAM,NSET,NPHAS,NIN,J,JN,IFIRST
        INTEGER NSS,NFILMS,NSETST,ILAST,IMORE,IFILMIN
        INTEGER OFILM,NF,IF1,LGP,ABSMAGP,ABSMAGPIN,JJ,CNT
        INTEGER IREDUN,ICNT,IRAN
        REAL PSIZE,DFMID1IN,DFMID2IN,ANGASTIN,PSSNR(*)
        REAL PSI(*),THETA(*),PHI(*),DFMID1(*),FSCT(*),R
        REAL DFMID2(*),ANGAST(*),SHX(*),SHY(*),THRESH(*)
        REAL PRESA(*),ABSMAGO,AMAGP(*),THETAM(*),ALPHA,D
        REAL DSTEP(*),RELMAG(*),AMEANPHAS,STHETA(*),DR
        REAL CMIN(*),CMAX(*),XM(*),YM(*),SX(*),SY(*),A
        REAL PSIM(*),SPSI(*),PS1,PS2,TS1,TS2,RISE,RX,RY
        REAL QM(*),OCC(*),PSX,ALGP(*),SIG(*),SIGN,RI,P(*)
	CHARACTER*200 FINPAR
	CHARACTER*200 LINE
        CHARACTER ASYM*3
        LOGICAL FOLD,FBOOST,NONSYM,MULTILINE
C*************************************************************************
        IF (NANG.GE.NNPART) THEN
          IF(IFLAG.GE.0) THEN
            WRITE(*,*)' INPUT PARAMETER FILE ?'
            READ(*,7006)FINPAR  
            WRITE(*,17006)FINPAR
          ELSE
7023        CONTINUE
            READ(*,*)NIN,ABSMAGPIN,IFILMIN,
     .               DFMID1IN,DFMID2IN,ANGASTIN,IMORE
            IF (IMORE.EQ.1) GOTO 7023
          ENDIF
        ELSE
C
        FOLD=.FALSE.
        MULTILINE=.FALSE.
	NANG1=NANG
        ICNT=1
        IF (NANG.EQ.0) INPB=0
      	NFILMS=0
        OFILM=-1
        NF=0
        CMAX(NSET)=0.0
        TNSET(NSET)=0
        XM(NSET)=0.0
        YM(NSET)=0.0
        SX(NSET)=0.0
        SY(NSET)=0.0
C        THETAM(NSET)=0.0
        STHETA(NSET)=0.0
        SPSI(NSET)=0.0
        QM(NSET)=0.0
        SIGN=0.0
        DO 33 I=2,NSAM/2
          PSSNR(I)=EXP(-(REAL(I-1)/NSAM*RI/NSAM)**2)
          IF ((PSIZE*NSAM/I.GT.30.0).OR.FBOOST) THEN
            FSCT(I)=1.0
          ELSE
            FSCT(I)=0.0
          ENDIF
33      CONTINUE
        PSSNR(1)=0.0
        FSCT(1)=1.0
C
      	IF(IFLAG.GE.0) THEN		
C				! needs input parameter file
                CMIN(NSET)=999.0
      		WRITE(*,*)' INPUT PARAMETER FILE ?'
      		READ(*,7006)FINPAR
      		WRITE(*,17006)FINPAR
7006		FORMAT(A200)
17006		FORMAT(3X,A200)
                OPEN(77,FILE=FINPAR,STATUS='OLD')
7009    	CONTINUE
      		READ(77,7007,END=7008)LINE
7007		FORMAT(A200)
C                IF (LINE(1:1).EQ.'C') GOTO 7009
      		IF (LINE(1:1).EQ.'C') THEN
                  IF ((NSET.NE.1).OR.(IFLAG.EQ.0)) GOTO 7009
                  DO 10 J=1,200-8
                    IF (LINE(J:J+8).EQ.'Part_SSNR') GOTO 20
10                CONTINUE
                  GOTO 7009
20                CONTINUE
         WRITE(*,*)' Reading Part_SSNR table from previous par file...'
                  DO 30 I=2,NSAM/2
C                    READ(77,7012,ERR=7013)J,FSCT(I),PSSNR(I)
                    READ(77,7012,ERR=7013)J,R,PSSNR(I)
7012                FORMAT(1X,I4,25X,F7.3,10X,F11.4)
                    IF (J.NE.I) STOP ' ERROR: Shell number mismatch'
C                    WRITE(*,7014)I,PSIZE*NSAM/REAL(I-1),FSCT(I),PSSNR(I)
                    WRITE(*,7014)I,PSIZE*NSAM/REAL(I-1),R,PSSNR(I)
7014                FORMAT( 'IRAD, RESOL, FSC, Part_SSNR = ',
     +                       I6,F7.1,2F11.4)
                    PSSNR(I)=PSSNR(I)**2
                    JJ=I
30                CONTINUE
7013              CONTINUE
                  PSSNR(1)=PSSNR(2)
C                  FSCT(1)=FSCT(2)
                  DO 31 I=JJ,NSAM/2
                    PSSNR(I)=0.0
C                    FSCT(I)=0.0
31                CONTINUE
                  WRITE(*,*)
                  GOTO 7009
                ENDIF
		NANG=NANG+1
7016            CONTINUE
                IF (FOLD) THEN
C      		  READ(LINE,*)ILIST(NANG),PSI(NANG),THETA(NANG),
      		  READ(LINE,7015)ILIST(NANG),PSI(NANG),THETA(NANG),
     +			PHI(NANG),SHX(NANG),SHY(NANG),ABSMAGO,FILM(NANG),
     +			DFMID1(NANG),DFMID2(NANG),ANGAST(NANG),PRESA(NANG)
7015              FORMAT(I7,5F8.2,F8.0,I6,2F9.1,F8.2,F7.2)
                  ABSMAGP=ABSMAGO
                  OCC(NANG)=100.0
                  SIG(NANG)=0.0
                ELSE
                  IF (SLEN2(LINE).LE.95) THEN
                    GOTO 99
                  ELSEIF (SLEN2(LINE).LE.128) THEN
C                    READ(LINE,*,ERR=99,IOSTAT=CNT)ILIST(NANG),
                    READ(LINE,70119,ERR=99,IOSTAT=CNT)ILIST(NANG),
     +                  PSI(NANG),THETA(NANG),PHI(NANG),
     +                  SHX(NANG),SHY(NANG),ABSMAGP,FILM(NANG),
     +                  DFMID1(NANG),DFMID2(NANG),ANGAST(NANG),
     +                  OCC(NANG),LGP,PRESA(NANG)
70119               FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I13,F8.2)
                    IF (CNT.NE.0) GOTO 99
                    SIG(NANG)=0.0
                    ALGP(NANG)=LGP
                    GOTO 98
                  ELSE
C                    READ(LINE,7011,ERR=99,IOSTAT=CNT)ILIST(NANG),
                    READ(LINE,*,ERR=99,IOSTAT=CNT)ILIST(NANG),
     +                  PSI(NANG),THETA(NANG),PHI(NANG),SHX(NANG),
     +                  SHY(NANG),ABSMAGP,FILM(NANG),DFMID1(NANG),
     +                  DFMID2(NANG),ANGAST(NANG),OCC(NANG),
     +                  LGP,SIG(NANG),PRESA(NANG)
7011                FORMAT(I7,3F8.2,2F10.2,I8,I6,2F9.1,2F8.2,I10,F11.4,
     +                     F8.2)
                    IF (CNT.NE.0) GOTO 99
                    ALGP(NANG)=LGP
                    GOTO 98
                  ENDIF
99                CONTINUE
                  IF (NANG.EQ.1) THEN
                    WRITE(*,*) ' Error reading par file.',
     +                         ' Trying pre v9 par file format...'
                    FOLD=.TRUE.
                    GOTO 7016
                  ELSE
                    WRITE(*,*)
     +                ' ERROR: something wrong in parameter file.'
                    WRITE(*,*) LINE
                    STOP 'Fatal error'
                  ENDIF
98                CONTINUE
                ENDIF
      		IF(ILIST(NANG).LE.0) THEN
      		  NANG=NANG-1
      		  WRITE(*,*) ' ***WARNING*** blank line in parameter file.'
      		  GO TO 7009
      		ENDIF
                IF (NANG.GT.1) THEN
                  IF((NANG.NE.NANG1+1).AND.
     +              (ABS(ILIST(NANG))-ABS(ILIST(NANG-1)).GT.1).AND.
     +              (NANG.LE.NNPART))
     +              WRITE(*,*) ' ***WARNING*** parameters missing for',
     +                       ' particle',ABS(ILIST(NANG-1))
                  IF (ABS(ILIST(NANG))-ABS(ILIST(NANG-1)).EQ.0) THEN
                    IF (NONSYM) THEN
                      ICNT=ICNT+1
                      MULTILINE=.TRUE.
                    ELSE
                      NANG=NANG-1
                      GO TO 7009
                    ENDIF
                  ENDIF
                  IF (ICNT.GT.IREDUN) THEN
                    WRITE(*,*) ' ***WARNING*** too many symmetry-',
     +              'related lines for particle',ABS(ILIST(NANG-1))
                    NANG=NANG-1
                    GO TO 7009
                  ENDIF
                  IF ((ABS(ILIST(NANG))-ABS(ILIST(NANG-1)).GE.1)
     +              .AND.NONSYM) THEN
                    IF (NANG.LE.NNPART) THEN
                      IF ((ICNT.NE.1).AND.(ICNT.LT.IREDUN))
     +                WRITE(*,*) ' ***WARNING*** too few symmetry-',
     +                'related lines for particle',ABS(ILIST(NANG-1))
                      IF (ICNT.LT.IREDUN) THEN
                        I=NANG+IREDUN-ICNT
                        ILIST(I)=ILIST(NANG)
                        PSI(I)=PSI(NANG)
                        THETA(I)=THETA(NANG)
                        PHI(I)=PHI(NANG)
                        SHX(I)=SHX(NANG)
                        SHY(I)=SHY(NANG)
                        FILM(I)=FILM(NANG)
                        AMAGP(I)=AMAGP(NANG)
                        DFMID1(I)=DFMID1(NANG)
                        DFMID2(I)=DFMID2(NANG)
                        ANGAST(I)=ANGAST(NANG)
                        OCC(I)=OCC(NANG)
                        SIG(I)=SIG(NANG)
                        PRESA(I)=PRESA(NANG)
                        CALL PARTITION(OCC(NANG-1),P,IREDUN,IRAN)
                        OCC(NANG-1)=P(1)
                        DO 60 I=NANG,NANG+IREDUN-ICNT-1
                          ILIST(I)=ILIST(NANG-1)
                          PSI(I)=PSI(NANG-1)
                          THETA(I)=THETA(NANG-1)
                          PHI(I)=PHI(NANG-1)
                          SHX(I)=SHX(NANG-1)
                          SHY(I)=SHY(NANG-1)
                          FILM(I)=FILM(NANG-1)
                          AMAGP(I)=AMAGP(NANG-1)
                          DFMID1(I)=DFMID1(NANG-1)
                          DFMID2(I)=DFMID2(NANG-1)
                          ANGAST(I)=ANGAST(NANG-1)
C                          OCC(I)=OCC(NANG-1)
                          OCC(I)=P(I-NANG+2)
                          SIG(I)=SIG(NANG-1)
                          PRESA(I)=PRESA(NANG-1)
                          PSIM(I)=PSI(NANG-1)
                          THETAM(I)=THETA(NANG-1)
60                      CONTINUE
                        NANG=NANG+IREDUN-ICNT
                        TNSET(NSET)=TNSET(NSET)+IREDUN-ICNT
                      ENDIF
                    ELSE
                      IF ((ICNT.LT.IREDUN).AND.(.NOT.MULTILINE)) THEN
                        TNSET(NSET)=TNSET(NSET)+IREDUN-ICNT
                      ENDIF
                    ENDIF
                    ICNT=1
                  ENDIF
                ENDIF
		ILIST(NANG)=ILIST(NANG)*NSS
      		IF (PSI(NANG).GT.180.0) PSI(NANG)=PSI(NANG)-360.0
      		IF (THETA(NANG).GT.180.0) THETA(NANG)=THETA(NANG)-360.0
      		IF (PHI(NANG).GT.180.0) PHI(NANG)=PHI(NANG)-360.0
      		IF (PSI(NANG).LT.-180.0) PSI(NANG)=PSI(NANG)+360.0
      		IF (THETA(NANG).LT.-180.0) THETA(NANG)=THETA(NANG)+360.0
      		IF (PHI(NANG).LT.-180.0) PHI(NANG)=PHI(NANG)+360.0
      		PSI(NANG)=PSI(NANG)/180.0*PI
                PSIM(NANG)=PSI(NANG)
      		THETA(NANG)=THETA(NANG)/180.0*PI
      		THETAM(NANG)=THETA(NANG)
                PSX=RELMAG(NSET)*DSTEP(NSET)*10000.0/ABS(ABSMAGP)
                IF (FOLD) THEN
                  PRESA(NANG)=COS(PRESA(NANG)/180.0*PI)*100.0
                  SHX(NANG)=SHX(NANG)*PSX
                  SHY(NANG)=SHY(NANG)*PSX
                ENDIF
      		IF (PRESA(NANG).GT.-100.0) THEN
      	  		AMEANPHAS=AMEANPHAS+PRESA(NANG)
      	  		NPHAS=NPHAS+1
      		ENDIF
      		PRESA(NANG)=PRESA(NANG)/100.0
      		SHX(NANG)=SHX(NANG)/PSX/NSAM*PI*2.0
      		SHY(NANG)=SHY(NANG)/PSX/NSAM*PI*2.0
C  Reset unrealistically large shifts
                IF (ABS(SHX(NANG)).GT.PI) SHX(NANG)=0.0
                IF (ABS(SHY(NANG)).GT.PI) SHY(NANG)=0.0
                OCC(NANG)=OCC(NANG)/100.0
C
C  The following IF/ENDIF block will not be executed if the helical
C  symmetry is given as 'HP' instead of 'H'. This allows the processing
C  of filaments that have a break in their helical lattice, such as
C  microtubules with a seam.
C
                IF ((ASYM(1:1).EQ.'H').AND.(THETA(NANG).NE.0.0)
     .            .AND.(ASYM(2:2).NE.'P')) THEN
      		  A=RELMAG(NSET)*DSTEP(NSET)*10000./
     .			(ABS(ABSMAGP)*PSIZE)
                  R=-RISE/PSIZE/NSAM*PI*2.0/A*SIN(THETA(NANG))
                  RX=-COS(PSI(NANG))
                  RY=SIN(PSI(NANG))
                  DR=(SHX(NANG)*RX+SHY(NANG)*RY)/R
                  SHX(NANG)=SHX(NANG)-NINT(DR)*RX*R
                  SHY(NANG)=SHY(NANG)-NINT(DR)*RY*R
                  PHI(NANG)=PHI(NANG)-NINT(DR)*ALPHA
                  IF (PHI(NANG).GT.180.0) PHI(NANG)=PHI(NANG)-360.0
                  IF (PHI(NANG).LT.-180.0) PHI(NANG)=PHI(NANG)+360.0
                ENDIF
C
      		PHI(NANG)=PHI(NANG)/180.0*PI
                IF (CMIN(NSET).GT.PRESA(NANG)) CMIN(NSET)=PRESA(NANG)
                IF (CMAX(NSET).LT.PRESA(NANG)) CMAX(NSET)=PRESA(NANG)
                TNSET(NSET)=TNSET(NSET)+1
                XM(NSET)=XM(NSET)+OCC(NANG)*SHX(NANG)
                YM(NSET)=YM(NSET)+OCC(NANG)*SHY(NANG)
                SX(NSET)=SX(NSET)+OCC(NANG)*SHX(NANG)**2
                SY(NSET)=SY(NSET)+OCC(NANG)*SHY(NANG)**2
                QM(NSET)=QM(NSET)+OCC(NANG)
                SIGN=SIGN+OCC(NANG)*SIG(NANG)
                IF (OFILM.NE.FILM(NANG)) THEN
                  IF (ASYM(1:1).EQ.'H') THEN
                    IF (NF.NE.0) THEN
                      PS1=PS1/NF
                      PS2=PS2/NF-PS1**2
                      IF (PS2.LT.0.0) PS2=0.0
                      TS1=TS1/NF
                      TS2=TS2/NF-TS1**2
                      IF (TS2.LT.0.0) TS2=0.0
                      DO 40 I=IF1,NANG-1
                        PSIM(I)=PS1
                        THETAM(I)=TS1
40                    CONTINUE
                      SPSI(NSET)=SPSI(NSET)+PS2
                      STHETA(NSET)=STHETA(NSET)+TS2
                    ENDIF
                    IF1=NANG
                    NF=1
                    PS1=PSI(NANG)
                    PS2=PSI(NANG)**2
                    TS1=THETA(NANG)
                    TS2=THETA(NANG)**2
                  ENDIF
                  OFILM=FILM(NANG)
                  NFILMS=NFILMS+1
                ELSE
                  IF (ASYM(1:1).EQ.'H') THEN
                    D=PSI(NANG)-PS1/NF
                    IF (D.GT.PI) D=D-2.0*PI
                    IF (D.LE.-PI) D=D+2.0*PI
                    PS2=PS2+(D+PS1/NF)**2
                    PS1=PS1+D+PS1/NF
                    D=THETA(NANG)-TS1/NF
                    IF (D.GT.PI) D=D-2.0*PI
                    IF (D.LE.-PI) D=D+2.0*PI
                    TS2=TS2+(D+TS1/NF)**2
                    TS1=TS1+D+TS1/NF
                    NF=NF+1
                  ENDIF
                ENDIF
C                THETAM(NSET)=THETAM(NSET)+THETA(NANG)
C                STHETA(NSET)=STHETA(NSET)+THETA(NANG)**2
                INPB=INPB+1
                IF (NANG.GT.1) THEN
                  IF ((FILM(NANG).NE.FILM(NANG-1)).OR.
     .              (NANG.EQ.ILAST)) THEN
                    IF (INPB.GT.NPB) NPB=INPB
                    INPB=0
                  ENDIF
                ENDIF
C      		IF (IFLAG.GE.3) THEN
C      			PSI(NANG)=0.0
C      			THETA(NANG)=0.0
C      			PHI(NANG)=0.0
C      			SHX(NANG)=0.0
C      			SHY(NANG)=0.0
C      			PRESA(NANG)=0.0
C      		ENDIF
      		ANGAST(NANG)=ANGAST(NANG)/180.0*PI
      		AMAGP(NANG)=RELMAG(NSET)*DSTEP(NSET)*10000./
     .			(ABS(ABSMAGP)*PSIZE)
      		IF(NANG.EQ.NANG1+1) THEN
      		   IF(AMAGP(NANG).LT.0.65.OR.AMAGP(NANG).GT.1.5) THEN
      			WRITE(*,7018) AMAGP(NANG),ABSMAGP,
     .				RELMAG(NSET),PSIZE,DSTEP(NSET)
7018			FORMAT(/' ***WARNING*** particle magnification and',
     .			' densitometer step size'/18X,
     .			' incompatible with selected pixel size'/
     .			60X,' AMAGP     = ',F10.4/
     .			60X,' ABSMAGP   = ',I10/
     .			60X,' RELMAG    = ',F10.4/
     .			60X,' PSIZE     = ',F10.4/
     .			60X,' DSTEP     = ',F10.4/)
      		   ELSE
c      			WRITE(*,*) AMAGP(NANG),ABSMAGP,
      			WRITE(*,7019) AMAGP(NANG),ABSMAGP,
     .			RELMAG(NSET),PSIZE,DSTEP(NSET)
7019			FORMAT(/' AMAGP     = ',F10.4/
     .                          ' ABSMAGP   = ',I10/
     .                          ' RELMAG    = ',F10.4/
     .                          ' PSIZE     = ',F10.4/
     .                          ' DSTEP     = ',F10.4/)
      		   ENDIF
      		ENDIF
		IF (NANG.GT.NNPART) NANG=NANG-1
                GOTO 7009
C
7008		CONTINUE
		CLOSE(77)
C
                IF (NONSYM.AND.(.NOT.MULTILINE)) THEN
                  IF ((ICNT.NE.1).AND.(ICNT.LT.IREDUN))
     +            WRITE(*,*) ' ***WARNING*** too few symmetry-',
     +            'related lines for particle',ABS(ILIST(NANG-1))
                  IF (ICNT.LT.IREDUN) THEN
                    DO 61 I=NANG+1,NANG+IREDUN-ICNT
                      ILIST(I)=ILIST(NANG)
                      PSI(I)=PSI(NANG)
                      THETA(I)=THETA(NANG)
                      PHI(I)=PHI(NANG)
                      SHX(I)=SHX(NANG)
                      SHY(I)=SHY(NANG)
                      FILM(I)=FILM(NANG)
                      AMAGP(I)=AMAGP(NANG)
                      DFMID1(I)=DFMID1(NANG)
                      DFMID2(I)=DFMID2(NANG)
                      ANGAST(I)=ANGAST(NANG)
                      OCC(I)=OCC(NANG)
                      SIG(I)=SIG(NANG)
                      PRESA(I)=PRESA(NANG)
                      PSIM(I)=PSI(NANG)
                      THETAM(I)=THETA(NANG)
61                  CONTINUE
                  ENDIF
                ENDIF
C
C                IF (PSSNR(1).EQ.1.0) THEN
C                  WRITE(*,*)
C                  WRITE(*,*)' ***WARNING*** No Part_SSNR table found!'
C                  WRITE(*,*)' Particle alignment will be less accurate.'
C                  WRITE(*,*)
C                ENDIF
                PS1=PS1/NF
                PS2=PS2/NF-PS1**2
                TS1=TS1/NF
                TS2=TS2/NF-TS1**2
                IF (ASYM(1:1).EQ.'H') THEN
                  DO 50 I=IF1,NANG
                    THETAM(I)=TS1
50                CONTINUE
                ENDIF
                STHETA(NSET)=STHETA(NSET)+TS2
C
      	ELSEIF(IFLAG.LT.0) THEN
                IF (NONSYM) THEN
                  WRITE(*,*)
     +              ' ERROR: startup not allowed with asymmetric',
     +              ' refinement.'
                  STOP 'Fatal error'
                ENDIF
                CMIN(NSET)=0.0
C				! create parameter file from scratch
      		FINPAR='    none                                          '
7022      	CONTINUE
		NANG2=NANG
      		READ(*,*)NIN,ABSMAGPIN,IFILMIN,DFMID1IN,
     .                   DFMID2IN,ANGASTIN,IMORE
      			NFILMS=NFILMS+1
      		NANG=NANG+NIN
                TNSET(NSET)=TNSET(NSET)+NIN
                IF (NANG.GT.NNPART) THEN
                  NIN=NIN-(NANG-NNPART)
                  NANG=NNPART
                ENDIF
      		WRITE(*,7021)NIN,ABSMAGPIN,IFILMIN,DFMID1IN,DFMID2IN,
     .			ANGASTIN,IFLAG,IMORE
7021	FORMAT(' Parameters created for this set of particles',/,
     .         '         number of particles ...........',I8,/,
     .         '         magnification .................',I10,/,
     .         '         film number ...................',I8,/,
     .         '         defocus 1 .....................',F10.1,/,
     .         '         defocus 2 .....................',F10.1,/,
     .         '         astigmatism angle..............',F10.1,/,
     .         '           IFLAG........................',I8,/,
     .         '           IMORE........................',I8)

      		DO 7020 J=1,NIN
      		      JN=NANG2+J
      			ILIST(JN)=JN*NSS
      			PSI(JN)=0.0
                        PSIM(JN)=0.0
      			THETA(JN)=0.0
                        IF (ASYM(1:1).EQ.'H') THETA(JN)=PI/2.0
      			THETAM(JN)=0.0
                        IF (ASYM(1:1).EQ.'H') THETAM(JN)=PI/2.0
      			PHI(JN)=0.0
      			SHX(JN)=0.0
      			SHY(JN)=0.0
      			AMAGP(JN)=RELMAG(NSET)*DSTEP(NSET)*10000./(ABSMAGPIN*PSIZE)
      			IF(JN.EQ.1) THEN
      			   IF(AMAGP(JN).LT.0.65.OR.AMAGP(JN).GT.1.5) THEN
      				WRITE(*,7018) AMAGP(JN),ABSMAGPIN,
     .				RELMAG(NSET),PSIZE,DSTEP(NSET)
      			   ELSE
      				WRITE(*,7019) AMAGP(JN),ABSMAGPIN,
     .				RELMAG(NSET),PSIZE,DSTEP(NSET)
      			   ENDIF
      			ENDIF
      			FILM(JN)=IFILMIN
      			DFMID1(JN)=DFMID1IN
      			DFMID2(JN)=DFMID2IN
      			ANGAST(JN)=ANGASTIN/180.0*PI
      			PRESA(JN)=0.0
                        OCC(JN)=1.0
                        SIG(JN)=0.0
7020		CONTINUE
      		IF ((NANG.LT.NNPART).AND.(IMORE.EQ.1)) THEN
                  GOTO 7022
                ELSEIF (IMORE.EQ.1) THEN
7024              CONTINUE
                  READ(*,*)NIN,ABSMAGPIN,IFILMIN,
     .                     DFMID1IN,DFMID2IN,ANGASTIN,IMORE
                  IF (IMORE.EQ.1) GOTO 7024
                ENDIF
      	ENDIF
C                TNSET(NSET)=(NANG-NANG1)
                IF (NONSYM) THEN
                  NANG=NANG+IREDUN-ICNT
                  TNSET(NSET)=TNSET(NSET)+IREDUN-ICNT
                  TNSET(NSET)=TNSET(NSET)/IREDUN
                  NANG=NANG/IREDUN
                ENDIF
                IF (TNSET(NSET).EQ.0)
     .            STOP ' ERROR: No particles in data set!'
                IF (QM(NSET).NE.0.0) THEN
                  XM(NSET)=XM(NSET)/QM(NSET)
                  YM(NSET)=YM(NSET)/QM(NSET)
                  SX(NSET)=SQRT(ABS(SX(NSET)/QM(NSET)-XM(NSET)**2))
                  SY(NSET)=SQRT(ABS(SY(NSET)/QM(NSET)-YM(NSET)**2))
                  SIGN=SIGN/QM(NSET)
C                  THETAM(NSET)=THETAM(NSET)/TNSET(NSET)
C                  STHETA(NSET)=SQRT(STHETA(NSET)
C     .                       /TNSET(NSET)-THETAM(NSET)**2)
                ENDIF
                IF(NFILMS.NE.0) SPSI(NSET)=SQRT(SPSI(NSET)/NFILMS)
                IF(NFILMS.NE.0) STHETA(NSET)=SQRT(STHETA(NSET)/NFILMS)
C
		IF(NFILMS.NE.0) WRITE(*,2001) NSET,TNSET(NSET),NFILMS,NANG
		IF(NFILMS.EQ.0) WRITE(*,2002) NSET,TNSET(NSET),NANG
2001		FORMAT(' Number of images in set',I3,' =',I8,' (on',
     .            I5,' new films)',', total particle parameters',
     .            ' so far =',I8)
2002		FORMAT(' Number of images in set',I3,' =',I5,
     .                 ', total particle parameters so far =',I6)
                IF ((SX(NSET).LT.1.0E-6).OR.(SY(NSET).LT.1.0E-6)) THEN
                  SX(NSET)=99999.0
                  SY(NSET)=99999.0
                  WRITE(*,2004)
2004              FORMAT(/' X,Y shifts distribution undetermined.',
     .                    ' X,Y ditribution function not used.',/)
                ELSE
		  WRITE(*,2003) CMIN(NSET)*100.0,CMAX(NSET)*100.0,
     .                        -XM(NSET)*NSAM/PI/2.0,
     .                        -YM(NSET)*NSAM/PI/2.0,
     .                        SX(NSET)*NSAM/PI/2.0,
     .                        SY(NSET)*NSAM/PI/2.0
2003		  FORMAT(/' Min, max score          =',2F12.3,/,
     .                  ' Using X,Y distribution function:',/,
     .                  ' Average X,Y shifts      =',2F12.3,/,
     .                  ' Sigma X,Y shifts        =',2F12.3,/)
                ENDIF
C
                IF (ASYM(1:1).EQ.'H') THEN
                  IF (STHETA(NSET).LT.1.0E-3) THEN
                    THETAM(NSET)=90.0/180.0*PI
                    STHETA(NSET)=2.0/180.0*PI
                    WRITE(*,2005)
2005                FORMAT(/' THETA angle distribution undetermined.',
     .                      ' THETA ditribution function not used.',/)
                  ELSE
		    WRITE(*,2006) STHETA(NSET)*180.0/PI
2006		    FORMAT(/' Using THETA distribution function:',/,
     .                      ' Sigma THETA angle        =',F8.3,/)
                  ENDIF
                  IF (SPSI(NSET).LT.1.0E-6) THEN
                    SPSI(NSET)=99999.0
                    WRITE(*,2007)
2007                FORMAT(/' PSI angle distribution undetermined.',
     .                      ' PSI ditribution function not used.',/)
                  ELSE
		    WRITE(*,2008) SPSI(NSET)*180.0/PI
2008		    FORMAT(/' Using PSI distribution function:',/,
     .                      ' Sigma PSI angle          =',F8.3,/)
                  ENDIF
                ENDIF
C
                IF (SIGN.EQ.0.0) THEN
                  WRITE(*,2016)
2016              FORMAT(/' Average noise sigma undetermined.'/)
                ELSE
                  WRITE(*,2015) SIGN
2015              FORMAT(/' Average noise sigma =',F8.3,/)
C  Reset unrealistically large or small sigmas
                  IF (NONSYM) THEN
                    DO 101 I=NANG1*IREDUN+1,NANG*IREDUN
C                      IF ((SIG(I)/SIGN).GT.10.0) SIG(I)=SIGN
                      IF (SIG(I).LE.0.0) THEN
                        SIG(I)=SIGN
C                      ELSEIF ((SIGN/SIG(I)).GT.10.0) THEN
C                        SIG(I)=SIGN
                      ENDIF
101                 CONTINUE
                  ELSE
                    DO 100 I=NANG1+1,NANG
C                      IF ((SIG(I)/SIGN).GT.10.0) SIG(I)=SIGN
                      IF (SIG(I).LE.0.0) THEN
                        SIG(I)=SIGN
C                      ELSEIF ((SIGN/SIG(I)).GT.10.0) THEN
C                        SIG(I)=SIGN
                      ENDIF
100                 CONTINUE
                  ENDIF
                ENDIF
C
C       	NANG - total number of images
C
        IF ((FOLD).AND.(ABS(THRESH(NSET)).GT.1.0)) THEN
          IF (ABS(THRESH(NSET)).EQ.90.0) THEN
            THRESH(NSET)=0.0
          ELSE
            THRESH(NSET)=COS(ABS(THRESH(NSET))/180.0*PI)*100.0
          ENDIF
        ELSEIF ((THRESH(NSET).LE.1.0).AND.
     .         (THRESH(NSET).GT.0.0)) THEN
          CALL PINC(NSET,THRESH,OCC,PRESA,NANG1,NANG)
        ENDIF
C
      		IF ((NANG.GE.IFIRST).AND.(NSETST.EQ.0)) THEN
      	  		NSETST=NSET
      	  		NNSET(NSET)=NANG-IFIRST+1
      	  		IF (NANG.GT.ILAST) NNSET(NSET)=ILAST-IFIRST+1
      		ELSEIF ((NANG1+1.GE.IFIRST).AND.(NANG.LE.ILAST)) THEN
      	  		NNSET(NSET)=NANG-NANG1
      		ELSEIF (NANG.GT.ILAST) THEN
      	  		NNSET(NSET)=ILAST-NANG1
      		ENDIF
C
      ENDIF
C
      RETURN
      END
C**************************************************************************
      SUBROUTINE PARTITION(X,P,N,IRAN)
C
      IMPLICIT NONE
C
      INTEGER N,IRAN,I
      REAL P(*),RANDOM,SUM,X
C
      SUM=0.0
      DO 10 I=1,N
        P(I)=RANDOM(IRAN)
        SUM=SUM+P(I)
10    CONTINUE
C
      DO 20 I=1,N
        P(I)=X*P(I)/SUM
20    CONTINUE
C
      RETURN
      END
