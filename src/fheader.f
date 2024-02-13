C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
      PROGRAM FHEADER
C**************************************************************************
C Program to get header info of image files
C**************************************************************************
      IMPLICIT NONE
C
      INTEGER NSAM,COMMAND_ARGUMENT_COUNT,NARG,SLEN2
      INTEGER IERR,LEN
      REAL PSIZE
      CHARACTER*200 FNAME
      CHARACTER VX*15,CFORM*1,ASYM*3,VXX*15
      CHARACTER CDATE*8,CTIME*10,CZONE*5
      LOGICAL EX
C**************************************************************************
C
      DATA VXX/'1.00 - 22.05.14'/
C     15 chars 'X.XX - XX.XX.XX' <--<--<--<--<--<--<--<--
C
      WRITE(*,1010)VXX
1010  FORMAT(/'  FHEADER ',A15,/,
     +       /'  Copyright 2013 Howard Hughes Medical Institute.',
     +       /'  All rights reserved.',
     +       /'  Use is subject to Janelia Farm Research Campus',
     +        ' Software Copyright 1.1',
     +       /'  license terms ( http://license.janelia.org/license',
     +        '/jfrc_copyright_1_1.html )'/)
C
C  Read input file
C
      NARG=COMMAND_ARGUMENT_COUNT()
C
      IF (NARG.EQ.0) THEN
        WRITE(*,*)' Image file for intput?'
        READ(*,1020)FNAME
1020    FORMAT(A200)
C        WRITE(*,1000)FNAME(1:SLEN2(FNAME))
1000    FORMAT(A)
      ELSE
        CALL GET_COMMAND_ARGUMENT(1,FNAME,LEN,IERR)
        IF (IERR.NE.0) THEN
          WRITE(*,*)' ERROR reading command line argument'
          GOTO 999
        ENDIF
C        WRITE(*,1000)FNAME(1:SLEN2(FNAME))
      ENDIF
C
      CALL GUESSF(FNAME,CFORM,EX)
      IF (EX) THEN
        CALL IOPEN(FNAME,10,CFORM,NSAM,NSAM,NSAM,'OLD',ASYM,PSIZE,VX)
        CALL ICLOSE(10)
      ELSE
        WRITE(*,*)' ERROR: File does not exist'
      ENDIF
C
999   CONTINUE
C
      END
