C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C ---------------------------------------------------------------------------------        
      SUBROUTINE DUMP(STD,VTD,NN1,NN2,NN3,IREDUN,FCREF,
     +          RAWBINS,WC2,RLIM,RI,RIC,HALFW,ALPHA,RISE,
     +          HSTART,PSIZE,MW,DALT,CFORM,ASYM,VX,F3D,FWEIGH,
     +          F3D1,F3D2,FPHA,FPOI,FOUTPAR1,QCP,ASUM,KSUM,
     +		PSUM,VSUM,A3DF,D3DF,S3DF,V3DF,
     +		IPAD,NNSTAT,NKSUM,NSYM,JSYM,SYMOP,
     +          INTERP,IFSC,INPIC,FBEAUT,FBFACT,IALL,NA,WALL,
     +          IWALL,VSN,VN,VVSN,VVN,NDOC1,NDOC2,NSET,IC,
     +          IEXCL,PALL,OALL,ICCOUNT,CCAVE)
C ---------------------------------------------------------------------------------        
C     Dump intermediate arrays to disk for parallelized reconstruction
C     Used in Frealign.
C ---------------------------------------------------------------------------------        
      IMPLICIT NONE
C
      INTEGER IREDUN,ISTAT,INTERP,NN1,NN2,NN3,IEXCL,NNSTAT
      INTEGER KSUM(NN3*NNSTAT+1),IPAD,NA,INPIC,IFSC,SLEN2
      INTEGER HSTART,NSYM,ISYMAX,NDOC1,NDOC2,NSET,I,IALL
      PARAMETER (ISYMAX=100)
      INTEGER JSYM(ISYMAX),IC(NN1),IWALL,ICCOUNT,IREC
      INTEGER*8 NKSUM
      DOUBLEPRECISION STD,VTD
      DOUBLEPRECISION VSN(NN1),VN(NN1),VVSN(NN1),VVN(NN1)
      REAL WC2,RLIM,RI,RIC,HALFW,ALPHA,RISE,PSIZE,MW
      REAL ASUM(NN3*NNSTAT+1),VSUM(NN3*NNSTAT+1),DALT
      REAL PSUM(NN3*NNSTAT+1),SYMOP(3,3,ISYMAX),WALL
      REAL S3DF(NN3),RAWBINS(NN1*2)
      REAL QCP(NN1/2*NN1*NN1*NNSTAT+1)
      REAL V3DF(NN3),PALL,OALL,CCAVE
      COMPLEX A3DF(NN3),D3DF(NN3)
      CHARACTER*200 F3D,FWEIGH,F3D1,F3D2,FPHA,FPOI
      CHARACTER*200 FOUTPAR1
      CHARACTER CFORM*1,ASYM*3,VX*15
      LOGICAL FCREF,FBEAUT,FBFACT
C ---------------------------------------------------------------------------------        
C
      WRITE(*,*)
      WRITE(*,*) ' Dumping intermediate results...'
      WRITE(*,1000)F3D(1:SLEN2(F3D))
1000  FORMAT(' Dump file: ',A)
      WRITE(*,*)
C
      DO 59 I=1,NSET
        CLOSE(NDOC1+NSET)
        CLOSE(NDOC2+NSET)
59    CONTINUE
      CALL ICLOSE(INPIC)
C
      OPEN(UNIT=INPIC,ACCESS='STREAM',FILE=F3D,ACTION='WRITE',
     .FORM='UNFORMATTED',STATUS='REPLACE')
C
      IREC=1
      IF (IFSC.EQ.0) THEN
        WRITE(INPIC,POS=IREC,IOSTAT=ISTAT)NN1,NN2,NN3,
     .    NNSTAT,IPAD,IREDUN,NKSUM,INTERP,IFSC,NSET,IWALL,
     .    WC2,RLIM,RI,RIC,
     .    HALFW,ALPHA,RISE,HSTART,PSIZE,MW,DALT,NSYM,JSYM,
     .    SYMOP,IALL,NA,WALL,CFORM,ASYM,VX,FCREF,FBEAUT,
     .    FBFACT,F3D,FWEIGH,F3D1,F3D2,FPHA,FPOI,FOUTPAR1,
     .    STD,VTD,VSN,VN,VVSN,VVN,ASUM,KSUM,PSUM,VSUM,QCP,
     .    A3DF,S3DF,D3DF,V3DF,RAWBINS,
     .    IC,IEXCL,PALL,OALL,ICCOUNT,CCAVE
      ELSE
        WRITE(INPIC,POS=IREC,IOSTAT=ISTAT)NN1,NN2,NN3,
     .    NNSTAT,IPAD,IREDUN,NKSUM,INTERP,IFSC,NSET,IWALL,
     .    WC2,RLIM,RI,RIC,
     .    HALFW,ALPHA,RISE,HSTART,PSIZE,MW,DALT,NSYM,JSYM,
     .    SYMOP,IALL,NA,WALL,CFORM,ASYM,VX,FCREF,FBEAUT,
     .    FBFACT,F3D,FWEIGH,F3D1,F3D2,FPHA,FPOI,FOUTPAR1,
     .    STD,VTD,VSN,VN,VVSN,VVN,ASUM,KSUM,PSUM,VSUM,QCP,
     .    A3DF,S3DF,RAWBINS,IC,IEXCL,PALL,OALL,
     .    ICCOUNT,CCAVE
      ENDIF
      IF (ISTAT.NE.0) THEN
        STOP ' ERROR: Dump file could not be written'
      ENDIF
C
      CLOSE(INPIC)
C
      RETURN
      END
