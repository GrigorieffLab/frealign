C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C -----------------------------------------------------------
	SUBROUTINE SHIFT2D(JC,NSAM,NSAMH,SHX,SHY,
     . 			DSHX,DSHY,OUTC)
C -----------------------------------------------------------
C           SHIFT PARTICLE INTO CENTRE
C	    Used in subroutine MATCH
C -----------------------------------------------------------
        REAL SHX,SHY
        COMPLEX OUTC(*)
	COMPLEX PSHFT
C -----------------------------------------------------------
            DO 109 L=1,JC
              LL=L-1
              DO 109 M=1,NSAM
                MM=M-1
                IF (MM.GE.JC) MM=MM-NSAM
                ISUM=(LL+MM)
                PSHFTR=1.0
                IF (MOD(ISUM,2).NE.0) PSHFTR=-1.0
                PHASE=-(SHX+DSHX)*LL-(SHY+DSHY)*MM
                PSHFT=CMPLX(COS(PHASE),SIN(PHASE))
C                IF (L.NE.JC) THEN
                  ID=L+JC*(M-1)
                  OUTC(ID)=OUTC(ID)*PSHFTR*PSHFT
C                ELSE
C                  OUTQ(M)=OUTQ(M)*PSHFTR*PSHFT
C                ENDIF
109         CONTINUE
	RETURN
	END
