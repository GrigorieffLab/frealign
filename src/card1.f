C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
	SUBROUTINE CARD1(VX,CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,
     .                FMATCH,FSTAT,FBEAUT,FCREF,FBFACT,IFSC,IMP,
     .                INTERP,IMEM,FDUMP,FBOOST)
C**************************************************************************
C Card 1: describes the overall program flow control logic cards
C         CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FMATCH,FSTAT,FBEAUT
C Called by Frealign.
C*************************************************************************
        IMPLICIT NONE
        INTEGER IFLAG,IEWALD,IFSC,IMP,INTERP,IMEM
        CHARACTER*15 VX
	CHARACTER*1 CFORM
	CHARACTER*200 INLINE
        LOGICAL FMAG,FDEF,FPART,FMATCH,FSTAT,FBOOST
        LOGICAL FASTIG,FBEAUT,FCREF,FBFACT,FDUMP
C*************************************************************************

        IF (IMP.GT.1) THEN
C
        WRITE(*,7000) VX,IMP
7000    FORMAT(/' FREALIGN - Fourier REconstruction and ALIGNment'/,
     *         /' V',A15,
     +         /'  Copyright 2013 Howard Hughes Medical Institute.',
     +         /'  All rights reserved.',
     +         /'  Use is subject to Janelia Farm Research Campus',
     +          ' Software Copyright 1.1',
     +         /'  license terms ( http://license.janelia.org/license',
     +          '/jfrc_copyright_1_1.html )'/,
     *         /' Parallel processing: NCPUS =   ',I8,/)
C
        ELSE
C
        WRITE(*,7001) VX
7001    FORMAT(/' FREALIGN - Fourier REconstruction and ALIGNment',
     *         /' V',A15,
     +         /'  Copyright 2013 Howard Hughes Medical Institute.',
     +         /'  All rights reserved.',
     +         /'  Use is subject to Janelia Farm Research Campus',
     +          ' Software Copyright 1.1',
     +         /'  license terms ( http://license.janelia.org/license',
     +          '/jfrc_copyright_1_1.html )'//)
C
        ENDIF
C
        INTERP=0
        FBFACT=.FALSE.
        FCREF=.FALSE.
        FSTAT=.FALSE.
        IFSC=0
C
        WRITE(*,*)' CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,',
     *          'FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FDUMP,IMEM,INTERP?'
        READ(*,'(A)')INLINE
        CFORM=INLINE(1:1)
        READ(INLINE(3:),*,ERR=96)IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,
     .           FBEAUT,FCREF,FBFACT,FMATCH,IFSC,FDUMP,IMEM,INTERP
        GOTO 98
96      CONTINUE
        WRITE(*,*) 'Card 1 error. Trying old CARD 1 input...'
        READ(INLINE(3:),*,ERR=97)IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,
     .                      FBEAUT,FCREF,FMATCH,IFSC,FDUMP,IMEM
        GOTO 98
97      CONTINUE
        WRITE(*,*) 'Card 1 error. Trying older CARD 1 input...'
        READ(INLINE(3:),*,ERR=99)IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,
     .                      FMATCH,FDUMP,FBEAUT,FCREF,IFSC
        GOTO 98
99      CONTINUE
        WRITE(*,*) 'Card 1 error. Trying oldest CARD 1 input...'
        READ(INLINE(3:),*)IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,
     .                      FMATCH,FDUMP,FBEAUT
98      CONTINUE
        WRITE(*,17001)CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,
     .                FBEAUT,FCREF,FBFACT,FMATCH,IFSC,FDUMP,IMEM,
     .                INTERP
17001   FORMAT(3X,A1,I3,4L2,I5,4L2,I3,L2,I3,I3)
      	IF(FDEF.AND.IFLAG.LT.0) THEN
      	WRITE(*,*)'Cannot refine defocus without first determining',
     .	' particle parameters at least roughly'
      	STOP 'FDEF true therefore not allowed'
      	ENDIF
      	IF(FMAG.AND.IFLAG.LT.0) THEN
          WRITE(*,*)
      	  WRITE(*,6999)
6999	  FORMAT('Cannot refine magnification of all particle',
     .	  ' images from same film'/' without first determining',
     .	  ' particle parameters at least roughly')
      	  STOP 'FMAG true therefore not allowed'
        ENDIF
        IF(IFLAG.EQ.0.AND.FMAG)THEN
          WRITE(*,*)
          WRITE(*,*) 'Cannot refine magnification with MODE=0;',
     +               ' FMAG set to F'
          WRITE(*,*)
          FMAG=.FALSE.
        ENDIF
        IF(IFLAG.EQ.0.AND.FDEF)THEN
          WRITE(*,*)
          WRITE(*,*) 'Cannot refine defocus with MODE=0;',
     +               ' FDEF set to F'
          WRITE(*,*)
          FDEF=.FALSE.
        ENDIF
        IF(IFLAG.EQ.0.AND.FASTIG)THEN
          WRITE(*,*)
          WRITE(*,*) 'Cannot refine astigmatism with MODE=0;',
     +               ' FASTIG set to F'
          WRITE(*,*)
          FASTIG=.FALSE.
      	ENDIF
      	IF(FCREF.AND.(IFSC.GT.0)) THEN
          WRITE(*,*)
      	  WRITE(*,*)'Cannot apply Wiener filter if IFSC not',
     .	  ' equal to 0'
      	  STOP 'FFILT true therefore not allowed'
      	ENDIF
        IF((IMEM.LT.0).OR.(IMEM.GT.3)) THEN
          WRITE(*,*)
          WRITE(*,*) 'IMEM must be betweeen 0 and 3;',
     +               ' IMEM set to 0'
          WRITE(*,*)
          IMEM=0
      	ENDIF
        FBOOST=.FALSE.
        IF (IFSC.LT.0) FBOOST=.TRUE.
C
	RETURN
	END
