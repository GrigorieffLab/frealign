C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C  -----------------------------------------------------------------------------
      SUBROUTINE OPMAPS(FPOI,FPHA,
     .			INSTAT,IPSTAT,
     .			IPOINT,INPIC,
     .			PSUM,
     .			NN1,JC,
     .			NSAM,NSAMH,PSIZE,
     .			A3DV,D3DV,S3DF,
     .			CFORM,ASYM,VX,SBUF,NNSTAT,IFSC)
C -----------------------------------------------------------------------------
C	O/P all maps
C       IOPEN, IWRITE, ICLOSE.
C       Used in Frealign.
C -----------------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER INSTAT,IPSTAT,IPOINT,INPIC,II,KK
	INTEGER NN1,JC,NSAM,NSAMH,N1,N2,N3,I,J,K,ID,IS
	INTEGER NNSTAT,IFSC
        REAL PSIZE,A3DV(*),D3DV(*),SBUF(*),S3DF(*),PSUM(*)
	CHARACTER*200 FPHA,FPOI
        CHARACTER*1 CFORM
        CHARACTER*3 ASYM
        CHARACTER*15 VX
C ----------------------------------------------------------------------------------
      	N1=NSAM
      	N2=NSAM
      	N3=NSAM
        IF ((NNSTAT.NE.0).AND.(IFSC.EQ.0)) THEN
      	  CALL IOPEN(FPOI,IPOINT,CFORM,N1,N2,N3,'NEW',ASYM,PSIZE,
     +		     VX)
        ENDIF
      	K=0
	DO 50 I=1,NSAM
	  DO 50 J=1,NSAM
      	    K=K+1
            DO 55 II=1,NSAM
      	      ID=II+(NSAM+2)*((J-1)+NSAM*(I-1))
              A3DV(ID)=A3DV(ID)/NSAM/NSAM/NSAM
              IF ((NNSTAT.NE.0).AND.(IFSC.EQ.0))
     +          D3DV(ID)=D3DV(ID)/NSAM/NSAM/NSAM
55          CONTINUE
      	    ID=1+(NSAM+2)*((J-1)+NSAM*(I-1))
      	    CALL IWRITE(INPIC,A3DV(ID),K)
            IF ((NNSTAT.NE.0).AND.(IFSC.EQ.0))
     +	      CALL IWRITE(IPOINT,D3DV(ID),K)
50	CONTINUE

	CALL ICLOSE(INPIC)
        IF ((NNSTAT.NE.0).AND.(IFSC.EQ.0))
     +    CALL ICLOSE(IPOINT)

      	K=0
	DO 51 I=1,NSAM
	  DO 51 J=1,NSAM
      	    K=K+1
      	    DO 52 KK=1,JC
      	      ID=KK+JC*((J-1)+NSAM*(I-1))
      	      SBUF(KK)=S3DF(ID)
52	    CONTINUE
      	    CALL IWRITE(INSTAT,SBUF,K)
51	CONTINUE

        IF (NNSTAT.NE.0) THEN
C
      	K=0
      	N1=JC
      	N2=NSAM
      	N3=NSAM
      	CALL IOPEN(FPHA,IPSTAT,CFORM,N1,N2,N3,'NEW',ASYM,PSIZE,VX)
      	SBUF(JC)=0.0
	DO 551 I=1,NSAM
	  DO 551 J=1,NSAM
      	    K=K+1
      	    DO 554 KK=1,JC
      	      ID=KK+JC*((J-1)+NSAM*(I-1))
      	      SBUF(KK)=PSUM(ID)
554	    CONTINUE
      	    CALL IWRITE(IPSTAT,SBUF,K)
551	CONTINUE

	CALL ICLOSE(IPSTAT)
C
        ENDIF
C
	CALL ICLOSE(INSTAT)
	END
