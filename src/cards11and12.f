C*
C* Copyright 2013 Howard Hughes Medical Institute.
C* All rights reserved.
C* Use is subject to Janelia Farm Research Campus Software Copyright 1.1
C* license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
C*
C**************************************************************************
	SUBROUTINE CARDS11AND12(NDOC1,NDOC2,NSET,CDATE,CTIME,
     .	CZONE,DTVAL,VX,CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,
     .  FMATCH,FDUMP,FBEAUT,FCREF,FBFACT,RI,RIC,PSIZE,WGH,XSTD,
     .  PBC,BOFF,ASYM,NSYM,SYMOP,IPMAX,ITMAX,DANGIN,IFIRST,ILAST,
     .  NSS,MAXSET,RELMAG,DSTEP,TARGET,THRESH,CS,AKV,RREC,RMAX1,
     .	RMAX2,RBFAC,FINPAT1,FINPAT2,FINPAR,ISYMAX,TX,TY,IFSC,
     .  DFSTD,IMEM,ALPHA,RISE,NU,HSTART,STIF,PMASK,DMASK,MW,
     .  INTERP,RCLAS,FOUTPAR1)
C**************************************************************************
C Read input cards 11 and 12
C Card 11       FOUTPAR - OUTPUT PARAMETER FILE
C Card 12       FOUTSH  - OUTPUT PARAMETER SHIFTS FILE, File name for file
C                       containing the parameter shifts for this refinement
C                       cycle (should become smaller if refinement converges).
C Calls DATE_AND_TIME.
C Called by Frealign.
C MW added Eald flag IEWALD
C*************************************************************************
#ifdef _NAG
        USE F90_UNIX
#endif
        IMPLICIT NONE
        INTEGER MAXSET,ISYMAX,NDOC1,NDOC2,NSET,PMASK(*)
        INTEGER ITMAX,IFIRST,ILAST,NU,IFSC,IMEM,SLEN2
        INTEGER IFLAG,I,J,K,NSYM,NSS,IEWALD,IPMAX,INTERP
        REAL PI,PBC,TX(*),TY(*),ALPHA,RISE,MW,DMASK(*)
	PARAMETER (PI=3.1415926535897)
	CHARACTER*15 VX
        CHARACTER*1 CFORM
	CHARACTER*200 FOUTPAR,FOUTPAR1,FOUTSH,FINPAR
        CHARACTER CDATE*8,CTIME*10,CZONE*5
	CHARACTER*3 ASYM
        CHARACTER*200 FINPAT1(*),FINPAT2(*)
        LOGICAL FMAG,FDEF,FPART,FMATCH,FDUMP,FCREF
        LOGICAL FASTIG,FBEAUT,FBFACT
	INTEGER DTVAL(8),HSTART
	REAL SYMOP(3,3,*),XSTD,RI,WGH,PSIZE,DANGIN,AKV
        REAL RELMAG(*),DSTEP(*),TARGET(*),BOFF
        REAL THRESH(*),CS(*),RIC,DFSTD(*),STIF
	REAL RREC(*),RMAX1(*),RMAX2(*),RBFAC(*),RCLAS(*)
C**************************************************************************
	WRITE(*,*)' OUTPUT PARAMETER FILE ?'
        READ(*,7006)FOUTPAR
        WRITE(*,17006)FOUTPAR
        IF (NSET.EQ.1) FOUTPAR1=FOUTPAR
        OPEN(NDOC1+NSET,FILE=FOUTPAR,STATUS='UNKNOWN')

        WRITE(*,*)' OUTPUT PARAMETER SHIFTS FILE ?'
        READ(*,7006)FOUTSH
        WRITE(*,17006)FOUTSH
        WRITE(*,*)
        OPEN(NDOC2+NSET,FILE=FOUTSH,STATUS='UNKNOWN')
7006            FORMAT(A200)
17006           FORMAT(3X,A200)

      	CALL DATE_AND_TIME(CDATE,CTIME,CZONE,DTVAL)

        IF (ASYM(1:1).EQ.'H') THEN
C
      	IF(NSET.EQ.1)WRITE(*,6705) CDATE(7:8),CDATE(5:6),CDATE(1:4),
     +		CTIME(1:2),CTIME(3:4),VX,CFORM,IFLAG,PMASK(3),
     +		PMASK(2),PMASK(1),PMASK(4),PMASK(5),
     +          (DMASK(I)*PSIZE,I=1,4),FMAG,FDEF,
     +		FASTIG,FPART,IEWALD,FBEAUT,FCREF,FBFACT,FMATCH,
     +		IFSC,FDUMP,IMEM,INTERP,RI*PSIZE,RIC*PSIZE,PSIZE,
     +          MW,WGH,XSTD,PBC,BOFF,ASYM,ALPHA,RISE,NU,HSTART,STIF
      	WRITE(NDOC1+NSET,6705) CDATE(7:8),CDATE(5:6),CDATE(1:4),
     +		CTIME(1:2),CTIME(3:4),VX,CFORM,IFLAG,PMASK(3),
     +		PMASK(2),PMASK(1),PMASK(4),PMASK(5),
     +          (DMASK(I)*PSIZE,I=1,4),FMAG,FDEF,
     +		FASTIG,FPART,IEWALD,FBEAUT,FCREF,FBFACT,FMATCH,
     +		IFSC,FDUMP,IMEM,INTERP,RI*PSIZE,RIC*PSIZE,PSIZE,
     +          MW,WGH,XSTD,PBC,BOFF,ASYM,ALPHA,RISE,NU,HSTART,STIF
	CALL FLUSH(NDOC1+NSET)
6705	FORMAT(  'C Date and time      ',A2,'-',A2,'-',A4,
     +           ',  ',A2,':',A2,'    Frealign V',A15,
     +          /'C Image format . . . . . . . . . . . . ',11X,A,
     +          /'C Mode . . . . . . . . . . . . . . . . ',I12,
     +          /'C PMASK for parameter refinement . . . ',5I2,
     +          /'C DMASK for parameter refinement . . . ',4F8.2,
     +          /'C Magnification refinement . . . . . . ',10X,L1,
     +          /'C Defocus refinement . . . . . . . . . ',10X,L1,
     +          /'C Astigmatism refinement . . . . . . . ',10X,L1,
     +          /'C Defocus ref. of individual particles ',10X,L1,
     +          /'C Ewald sphere correction. . . . . . . ',I12,
     +          /'C Beautify the final real space map. . ',10X,L1,
     +          /'C Apply Wiener filter to final map . . ',10X,L1,
     +          /'C B-factor correction of final map . . ',10X,L1,
     +          /'C Write out matching projections . . . ',10X,L1,
     +          /'C Calculate FSPR and FSC curves. . . . ',I12,
     +          /'C Dump intermediate 3D files . . . . . ',10X,L1,
     +          /'C Memory/speed optimization. . . . . . ',I12,
     +          /'C 3D interpolation . . . . . . . . . . ',I12,
     +          /'C Outer Radius of object [Angstroms] . ',F12.2,
     +          /'C Inner Radius of object [Angstroms] . ',F12.2,
     +          /'C Pixel size [Angstroms] . . . . . . . ',F12.5,
     +          /'C Molecular mass [kDa] . . . . . . . . ',F12.3,
     +          /'C % Amplitude contrast . . . . . . . . ',F12.2,
     +          /'C STD level for 3D mask. . . . . . . . ',F12.2,
     +          /'C Score / B factor constant. . . . . . ',F12.2,
     +          /'C Average score for weighting. . . . . ',F12.2,
     +          /'C Symmetry card as input . . . . . . . ',11X,A3,
     +          /'C   Helical rotation per subunit . . . ',F12.2,
     +          /'C   Helical rise per subunit . . . . . ',F12.2,
     +          /'C   Number of subunits to average. . . ',I12,
     +          /'C   Number of starts . . . . . . . . . ',I12,
     +          /'C   Stiffness parameter. . . . . . . . ',F12.2)
C
        ELSE
C
      	IF(NSET.EQ.1)WRITE(*,6700) CDATE(7:8),CDATE(5:6),CDATE(1:4),
     +		CTIME(1:2),CTIME(3:4),VX,CFORM,IFLAG,PMASK(3),
     +		PMASK(2),PMASK(1),PMASK(4),PMASK(5),
     +          (DMASK(I)*PSIZE,I=1,4),FMAG,FDEF,
     +		FASTIG,FPART,IEWALD,FBEAUT,FCREF,FBFACT,FMATCH,
     +		IFSC,FDUMP,IMEM,INTERP,RI*PSIZE,RIC*PSIZE,PSIZE,
     +          MW,WGH,XSTD,PBC,BOFF,ASYM,NSYM
      	WRITE(NDOC1+NSET,6700) CDATE(7:8),CDATE(5:6),CDATE(1:4),
     +		CTIME(1:2),CTIME(3:4),VX,CFORM,IFLAG,PMASK(3),
     +		PMASK(2),PMASK(1),PMASK(4),PMASK(5),
     +          (DMASK(I)*PSIZE,I=1,4),FMAG,FDEF,
     +		FASTIG,FPART,IEWALD,FBEAUT,FCREF,FBFACT,FMATCH,
     +		IFSC,FDUMP,IMEM,INTERP,RI*PSIZE,RIC*PSIZE,PSIZE,
     +          MW,WGH,XSTD,PBC,BOFF,ASYM,NSYM
	CALL FLUSH(NDOC1+NSET)
6700	FORMAT(  'C Date and time      ',A2,'-',A2,'-',A4,
     +           ',  ',A2,':',A2,'    Frealign V',A15,
     +          /'C Image format . . . . . . . . . . . . ',9X,A,
     +          /'C Mode . . . . . . . . . . . . . . . . ',I10,
     +          /'C PMASK for parameter refinement . . . ',5I2,
     +          /'C DMASK for parameter refinement . . . ',4F8.2,
     +          /'C Magnification refinement . . . . . . ',8X,L1,
     +          /'C Defocus refinement . . . . . . . . . ',8X,L1,
     +          /'C Astigmatism refinement . . . . . . . ',8X,L1,
     +          /'C Defocus ref. of individual particles ',8X,L1,
     +          /'C Ewald sphere correction. . . . . . . ',I10,
     +          /'C Beautify the final real space map. . ',8X,L1,
     +          /'C Apply Wiener filter to final map . . ',8X,L1,
     +          /'C B-factor correction of final map . . ',8X,L1,
     +          /'C Write out matching projections . . . ',8X,L1,
     +          /'C Calculate FSPR and FSC curves. . . . ',I10,
     +          /'C Dump intermediate 3D files . . . . . ',8X,L1,
     +          /'C Memory/speed optimization. . . . . . ',I10,
     +          /'C 3D interpolation . . . . . . . . . . ',I10,
     +          /'C Outer Radius of object [Angstroms] . ',F10.2,
     +          /'C Inner Radius of object [Angstroms] . ',F10.2,
     +          /'C Pixel size [Angstroms] . . . . . . . ',F10.5,
     +          /'C Molecular mass [kDa] . . . . . . . . ',F10.3,
     +          /'C % Amplitude contrast . . . . . . . . ',F10.2,
     +          /'C STD level for 3D mask. . . . . . . . ',F10.2,
     +          /'C Score / B factor constant. . . . . . ',F10.2,
     +          /'C Average score for weighting. . . . . ',F10.2,
     +          /'C Symmetry card as input . . . . . . . ',9X,A3,
     +          /'C Number of symmetry operators . . . . ',I10)
      	DO 801 I=1,NSYM
      	  DO 802 K=1,3
      	  IF(NSET.EQ.1) WRITE(*,6702) I,(SYMOP(J,K,I),J=1,3)
802       WRITE(NDOC1+NSET,6702) I,(SYMOP(J,K,I),J=1,3)
	  CALL FLUSH(NDOC1+NSET)
6702	  FORMAT('C',I3,3F12.4)
801	CONTINUE
C
        ENDIF
C
      		IF(IFLAG.EQ.2.OR.IABS(IFLAG).EQ.4) THEN
      		 WRITE(*,6717) ITMAX
      		 WRITE(NDOC1+NSET,6717) ITMAX
6717		 FORMAT('C No. cycles for randomised search . . ',I10)
      		 WRITE(*,6719) IPMAX
      		 WRITE(NDOC1+NSET,6719) IPMAX
6719		 FORMAT('C No. of search peaks to refine. . . . ',I10)
	         CALL FLUSH(NDOC1+NSET)
      		ENDIF
6718		 FORMAT('C Angular step size for search . . . . ',F10.1)
      		IF(IFLAG.EQ.3.OR.IABS(IFLAG).EQ.4)
     +            WRITE(NDOC1+NSET,6718) DANGIN
	        CALL FLUSH(NDOC1+NSET)
      	WRITE(NDOC1+NSET,6701)IFIRST,ILAST,
     +	      RELMAG(NSET),DSTEP(NSET),TARGET(NSET),THRESH(NSET),
     +	      CS(NSET),AKV/1000.0,TX(NSET),TY(NSET),
     +	      PSIZE/RREC(NSET),PSIZE/RMAX1(NSET),PSIZE/RMAX2(NSET),
     +	      PSIZE/RCLAS(NSET),DFSTD(NSET),RBFAC(NSET),
     +        FINPAT1(NSET)(1:SLEN2(FINPAT1(NSET)))
	CALL FLUSH(NDOC1+NSET)
6701	FORMAT(  'C First particle . . . . . . . . . . . ',I10,
     +          /'C Last particle. . . . . . . . . . . . ',I10,
     +          /'C Relative magnification . . . . . . . ',F12.4,
     +          /'C Densitometer step size (microns) . . ',F9.1,
     +          /'C Score target . . . . . . . . . . . . ',F10.2,
     +          /'C Score threshold. . . . . . . . . . . ',F10.2,
     +          /'C Cs [mm]. . . . . . . . . . . . . . . ',F10.2,
     +          /'C Voltage [kV] . . . . . . . . . . . . ',F10.2,
     +          /'C Beam tilt Tx, Ty [mrad]. . . . . . . ',2F10.2,
     +          /'C Resolution of reconstruction . . . . ',F11.3,
     +          /'C Low resol. limit refinement. . . . . ',F11.3,
     +          /'C High resol. limit refinement . . . . ',F11.3,
     +          /'C High resol. limit classification . . ',F11.3,
     +          /'C Defocus uncertainty. . . . . . . . . ',F11.3,
     +          /'C B-factor for parameter refinement. . ',F11.3,
     +          /'C Input image stack           ',A)
      	IF (FMATCH) WRITE(NDOC1+NSET,6703)
     +    FINPAT2(NSET)(1:SLEN2(FINPAT2(NSET)))
	CALL FLUSH(NDOC1+NSET)
6703	FORMAT('C Matching projections stack ',A)
      	WRITE(NDOC1+NSET,6704)FINPAR(1:SLEN2(FINPAR)),
     +                        FOUTPAR(1:SLEN2(FOUTPAR)),
     +                        FOUTSH(1:SLEN2(FOUTSH))
	CALL FLUSH(NDOC1+NSET)
6704	FORMAT(  'C Input parameter file        ',A,
     +          /'C Output parameter file       ',A,
     +          /'C Output shifts file          ',A)
      	CS(NSET)=CS(NSET)*(10.0**7.0)
      	TARGET(NSET)=TARGET(NSET)/100.0
      	THRESH(NSET)=THRESH(NSET)/100.0
      	NSET=NSET+1
      	NSS=-NSS
	RETURN
	END
