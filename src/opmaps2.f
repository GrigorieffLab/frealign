C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C  -----------------------------------------------------------------------------
	SUBROUTINE OPMAPS2(F3D1,F3D2,I3D1,I3D2,
     .			JC,NSAM,NSAMH,PSIZE,A3DV,B3DV,
     .			CFORM,ASYM,VX,SBUF)
C -----------------------------------------------------------------------------
C	O/P halfset maps (need to be shifted into center of box)
C       IOPEN, IWRITE, ICLOSE.
C       Used in Frealign.
C -----------------------------------------------------------------------------
        IMPLICIT NONE
        INTEGER I3D1,I3D2,II,KK
	INTEGER JC,NSAM,NSAMH,N1,N2,N3,I,J,K,ID
	INTEGER L,M,N
        REAL PSIZE,A3DV(*),B3DV(*),SBUF(*)
	CHARACTER*200 F3D1,F3D2
        CHARACTER*1 CFORM
        CHARACTER*3 ASYM
        CHARACTER*15 VX
C ----------------------------------------------------------------------------------
      	N1=NSAM
      	N2=NSAM
      	N3=NSAM
      	CALL IOPEN(F3D1,I3D1,CFORM,N1,N2,N3,'NEW',ASYM,PSIZE,
     +		   VX)
      	CALL IOPEN(F3D2,I3D2,CFORM,N1,N2,N3,'NEW',ASYM,PSIZE,
     +		   VX)
C
      	K=0
	DO 50 I=1,NSAM
C          L=I-NSAMH
C          IF (L.LT.1) L=L+NSAM
	  DO 50 J=1,NSAM
C            M=J-NSAMH
C            IF (M.LT.1) M=M+NSAM
      	    K=K+1
            DO 55 II=1,NSAM
C              N=II-NSAMH
C              IF (N.LT.1) N=N+NSAM
C      	      ID=N+(NSAM+2)*((M-1)+NSAM*(L-1))
      	      ID=II+(NSAM+2)*((J-1)+NSAM*(I-1))
              SBUF(II)=2.0*A3DV(ID)
55          CONTINUE
      	    CALL IWRITE(I3D1,SBUF,K)
            DO 56 II=1,NSAM
C              N=II-NSAMH
C              IF (N.LT.1) N=N+NSAM
C      	      ID=N+(NSAM+2)*((M-1)+NSAM*(L-1))
      	      ID=II+(NSAM+2)*((J-1)+NSAM*(I-1))
              SBUF(II)=2.0*B3DV(ID)
56          CONTINUE
     	    CALL IWRITE(I3D2,SBUF,K)
50	CONTINUE
C
	CALL ICLOSE(I3D1)
	CALL ICLOSE(I3D2)
C
	END
