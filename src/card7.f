C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
        SUBROUTINE CARD7(NSET,MAXSET,RREC,RMAX1,RMAX2,
     +                   RCLAS,DFSTD,RBFAC)
C**************************************************************************
C Read input card 7 
C RREC,RMAX1,RMAX2,RCLAS,DFSTD,RBFACT - map resoln, refine low/high,
C                        class. resoln, defocus uncertainty, B-factor
C Called by Frealign.
C*************************************************************************
        IMPLICIT NONE
        INTEGER NSET,MAXSET
	REAL RREC(*),RMAX1(*),RMAX2(*),DFSTD(*),RBFAC(*),RCLAS(*)
        CHARACTER*80 INLINE
C*************************************************************************
        WRITE(*,*)' RESOLUTION RECONST., REFINE LOW/HIGH,',
     +            ' CLASSIFY, DFSTD, RBFACT ?'
        READ(*,'(A)')INLINE
        READ(INLINE,*,ERR=99) RREC(NSET),RMAX1(NSET),RMAX2(NSET),
     +            RCLAS(NSET),DFSTD(NSET),RBFAC(NSET)
        GOTO 98 
99      CONTINUE
        WRITE(*,*)
        WRITE(*,*) 'Card 7 error. Trying old CARD 7 input...'
        WRITE(*,*)
        READ(INLINE,*,ERR=97) RREC(NSET),RMAX1(NSET),RMAX2(NSET),
     +            RBFAC(NSET)
        DFSTD(NSET)=100.0
        RCLAS(NSET)=RMAX2(NSET)
        GOTO 98
97      CONTINUE
        READ(INLINE,*) RREC(NSET),RMAX1(NSET),RMAX2(NSET),
     +            DFSTD(NSET),RBFAC(NSET)
        RCLAS(NSET)=RMAX2(NSET)
98      CONTINUE
        DFSTD(NSET)=ABS(DFSTD(NSET))
        WRITE(*,17005) RREC(NSET),RMAX1(NSET),RMAX2(NSET),
     +            RCLAS(NSET),DFSTD(NSET),RBFAC(NSET)
17005   FORMAT(6F12.1)
	RETURN
	END
