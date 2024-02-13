C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
c ---------------------------------------------------------------------
	SUBROUTINE OPRESSTATHALF(NDOC1,VSN,VN,RBINS,PF,PR,ISTP,
     .		NSTP,NSAM,RLIM,PSIZE,FSC,QF,QC,NPIX,NPIXT,
     . 		CDATE,CTIME,CZONE,DTVAL,VX,NN1,IC,VVSN,VVN,IALL,
     .          NA,PSSNR,FMASK,FPART,BFACT,MASK,WALL)
c ---------------------------------------------------------------------
c	O/P resolution statistics between two halves
c	Calls DATE_AND_TIME.
c	Used by Frealign.
c ---------------------------------------------------------------------
#ifdef _NAG
        USE F90_UNIX
#endif
        IMPLICIT NONE
        REAL PI
	PARAMETER  (PI=3.1415926535897)
	INTEGER DTVAL(8),I,NDOC1,NSTP,NSAM,NN1,ISTP,IC(*),IALL
        INTEGER ISUMV,ISUMP,NA
        CHARACTER CDATE*8,CTIME*10,CZONE*5
	CHARACTER*15 VX
        REAL PSIZE,RINGRAD,RLIM,SQRTSSNR,SSNR,FMASK,FPART
	REAL RBINS(*),RVSN,RVN,CC,RVVSN,RVVN,SCC,AVPFSC,WALL
	REAL QF(*),PF(*),PR(*),RL,ERC,PSSNR(*),PFSC,MASK
	REAL FSC(*),QC(*),AVFSPR,AVFSC,AVQ,AVQRAN,BFACT
        REAL AVSSNR,AVCFREE,AVEXPC,AVSIGC,AVLL,AVPSSNR
        DOUBLEPRECISION VSN(*),VN(*),VVSN(*),VVN(*)
	INTEGER NPIX(*),NPIXT(*)
c ---------------------------------------------------------------------
        AVFSPR=0.0
        AVFSC=0.0
        AVPFSC=0.0
        AVQ=0.0
        AVQRAN=0.0
        AVSSNR=0.0
        AVCFREE=0.0
        AVEXPC=0.0
        AVSIGC=0.0
        AVLL=0.0
        AVPSSNR=0.0
        ISUMV=0
        ISUMP=0
        WRITE(NDOC1+1,9995) FMASK,FPART
C        WRITE(*,9995),FMASK,FPART
9995    FORMAT('C'/
     +         'C  Fraction_mask, Fraction_particle = ',2F11.5/
     +         'C')
        IF (BFACT.NE.0.0) WRITE(NDOC1+1,9994) BFACT
9994    FORMAT('C  B-factor applied for sharpening (in A^2) = ',
     +         F14.5/'C')
        IF (MASK.NE.0.0) WRITE(NDOC1+1,9993) MASK
9993    FORMAT('C  Resolution limited by cosine filter to (in A) = ',
     +         F9.5/'C')
        WRITE(NDOC1+1,7013)
        WRITE(*,7013)
        CALL FLUSH(NDOC1+1)
        CALL FLUSH(6)
7012    FORMAT('C',49X,'sqrt',7X,'sqrt'/
     +         'C  NO.  RESOL  RING RAD   FSPR    FSC ',
     +         ' Part_FSC  Part_SSNR  Rec_SSNR       CC ',
     +         '  EXP. C    SIG C  ERFC  TOTVOX  APDIF')
7013    FORMAT('C',49X,'sqrt',7X,'sqrt'/
     +         'C  NO.  RESOL  RING RAD   FSPR    FSC ',
     +         ' Part_FSC  Part_SSNR  Rec_SSNR       CC ',
     +         '  EXP. C    SIG C  ERFC  TOTVOX')
        DO 62 I=2,NSTP
          IF(RBINS(2*I).GT.0.0) RBINS(2*I-1)=RBINS(2*I-1)/RBINS(2*I)
          RVSN=0.0
          RVN=0.0
          RVVSN=0.0
          RVVN=0.0
          IF(IALL.NE.0) THEN
            RVSN=REAL(VSN(I)/IALL)
            RVN=REAL(VN(I)/IALL)
            RVVSN=REAL(VVSN(I)/IALL)
            RVVN=REAL(VVN(I)/IALL)
            RVVSN=SQRT((RVVSN-RVSN**2)/IALL)
            RVVN=SQRT((RVVN-RVN**2)/IALL)
          ENDIF
          CC=0.0
          SCC=0.0
          RL=0.0
          IF(RVSN.NE.0.0) THEN
C            SSNR=(RVSN/RVN-1.0)*(REAL(NA)/NSAM/NSAM)
C            SSNR=(RVSN/RVN-1.0)
C            CC=SSNR/(1.0+SSNR)
C            CC=(RVSN-RVN)/RVSN
            IF (RVSN.GE.RVN) CC=SQRT((RVSN-RVN)/RVSN)
C            SCC=SQRT(RVVSN**2*(1.0-CC)**2+RVVN**2)/RVSN
            SCC=0.2*(ABS(RBINS(2*I-1))+ABS(CC))/2.0
            IF(SCC.NE.0.0)RL=ERC(ABS(RBINS(2*I-1)-CC)/SCC/SQRT(2.0))
          ENDIF
C          print *,I,RVSN,RVN,RVVSN,RVVN,CC,SCC,IC(I)
          SSNR=2.0*FSC(I)/(1.0-FSC(I)+0.0001)
          PFSC=0.0
          IF (FSC(I).GE.1.0) THEN
            PFSC=1.0
          ELSEIF (FPART.NE.0.0) THEN
            SSNR=FMASK/FPART*SSNR
            PFSC=SSNR/(ABS(SSNR)+2.0)
          ENDIF
          IF(PF(I)-1.0.GT.0.0) THEN
                SQRTSSNR=SQRT(PF(I)-1.0)
          ELSE
                SQRTSSNR=SQRT(ABS(SSNR))
          ENDIF
          RINGRAD=REAL(I-1)*ISTP/NSAM
          IF(RINGRAD.GT.RLIM) GO TO 62
            WRITE(NDOC1+1,7014)I,PSIZE/RINGRAD,RINGRAD,PR(I)*180.0/PI,
     +          FSC(I),PFSC,SQRT(PSSNR(I)*WALL),SQRTSSNR,
     +          RBINS(2*I-1),CC,SCC,RL,
     +          NPIXT(I)
            CALL FLUSH(NDOC1+1)
            WRITE(*,7014)I,PSIZE/RINGRAD,RINGRAD,PR(I)*180.0/PI,
     +          FSC(I),PFSC,SQRT(PSSNR(I)*WALL),SQRTSSNR,
     +          RBINS(2*I-1),CC,SCC,RL,
     +          NPIXT(I)
            CALL FLUSH(6)
C          END IF
7014      FORMAT('C',I4,F8.2,F10.4,F7.2,F7.3,F10.3,F11.4,F10.2,
     +           3F9.4,F6.2,I8)
          AVFSPR=AVFSPR+NPIX(I)*PR(I)*180.0/PI
          AVFSC=AVFSC+NPIX(I)*FSC(I)
          AVPFSC=AVPFSC+NPIX(I)*PFSC
          AVQ=AVQ+NPIX(I)*QF(I)
          AVQRAN=AVQRAN+NPIX(I)*QC(I)
          AVSSNR=AVSSNR+NPIX(I)*SQRTSSNR**2
          AVCFREE=AVCFREE+IC(I)*RBINS(2*I-1)
          AVEXPC=AVEXPC+IC(I)*CC
          AVSIGC=AVSIGC+IC(I)*SCC
          AVLL=AVLL+IC(I)*RL
          AVPSSNR=AVPSSNR+NPIX(I)*PSSNR(I)*WALL
          ISUMV=ISUMV+NPIX(I)
          ISUMP=ISUMP+IC(I)
62      CONTINUE

        IF (ISUMP.EQ.0) ISUMP=1
        WRITE(NDOC1+1,7015)AVFSPR/ISUMV,AVFSC/ISUMV,AVPFSC/ISUMV,
     +          SQRT(AVPSSNR/ISUMV),SQRT(AVSSNR/ISUMV),
     +          AVCFREE/ISUMP,AVEXPC/ISUMP,AVSIGC/ISUMP,AVLL/ISUMP
        CALL FLUSH(NDOC1+1)
        WRITE(*,7015)AVFSPR/ISUMV,AVFSC/ISUMV,AVPFSC/ISUMV,
     +          SQRT(AVPSSNR/ISUMV),SQRT(AVSSNR/ISUMV),
     +          AVCFREE/ISUMP,AVEXPC/ISUMP,AVSIGC/ISUMP,AVLL/ISUMP
        CALL FLUSH(6)
7015      FORMAT('C  Average             ',F7.2,F7.3,F10.3,F11.4,
     +            F10.2,3F9.4,F6.2)
        CALL DATE_AND_TIME(CDATE,CTIME,CZONE,DTVAL)
        WRITE(NDOC1+1,6710) CDATE(7:8),CDATE(5:6),CDATE(1:4),
     +          CTIME(1:2),CTIME(3:4),VX
        CALL FLUSH(NDOC1+1)
6710    FORMAT('C'/,
     +         'C Date and time      ',A2,'-',A2,'-',A4,',  ',A2,
     +         ':',A2,'    Frealign V',A15)

	RETURN
	END
