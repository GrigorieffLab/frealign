      SUBROUTINE CROP_FT(NIN,NOUT,DINC,DINQ,DOUTC,DOUTQ)
C
      IMPLICIT NONE
C
      INTEGER NIN,NOUT,I,J,ID,IS,NINH,NOUTH,ND
      COMPLEX DINC(*),DINQ(*),DOUTC(*),DOUTQ(*)
C
      NINH=NIN/2
      NOUTH=NOUT/2
      ND=NIN-NOUT
      IF (ND.LT.0) RETURN
C
      DO 30 J=1,NOUTH
        IF (NOUT.EQ.NIN) THEN
          DOUTQ(J)=DINQ(J)
        ELSE
          DOUTQ(J)=DINC(J)
        ENDIF
        DOUTQ(J)=0.0
30    CONTINUE
C
      DO 40 J=NOUTH+1,NOUT
        IF (NOUT.EQ.NIN) THEN
          DOUTQ(J)=DINQ(J)
        ELSE
          DOUTQ(J)=DINC(ND+J)
        ENDIF
        DOUTQ(J)=0.0
40    CONTINUE
C
      DO 10 J=1,NOUTH
        DO 10 I=1,NOUTH
          ID=I+NOUTH*(J-1)
          IS=I+NINH*(J-1)
          DOUTC(ID)=DINC(IS)
10    CONTINUE
C
      DO 20 J=NOUTH+1,NOUT
        DO 20 I=1,NOUTH
          ID=I+NOUTH*(J-1)
          IS=I+NINH*(ND+J-1)
          DOUTC(ID)=DINC(IS)
20    CONTINUE
C
      RETURN
      END
