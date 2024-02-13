# 
# ------ Lists source files for use by Makefile  ------
#  
.SUFFIXES :  .f .o .a .f90

ELEMENTS =  \
 $(LIB)(a3d3.o) \
 $(LIB)(a3d3_s.o) \
 $(LIB)(ainterpo3ds.o) \
 $(LIB)(ainterpo3dbig.o) \
 $(LIB)(applycref.o) \
 $(LIB)(beautify.o) \
 $(LIB)(beautify_s.o) \
 $(LIB)(bfactorsub.o) \
 $(LIB)(boxft_lut.o) \
 $(LIB)(calcfx.o) \
 $(LIB)(calcsinc.o) \
 $(LIB)(card10.o) \
 $(LIB)(cards11and12.o) \
 $(LIB)(cards13and14.o) \
 $(LIB)(cards15to18.o) \
 $(LIB)(card1.o) \
 $(LIB)(card2.o) \
 $(LIB)(card3.o) \
 $(LIB)(card4.o) \
 $(LIB)(card5.o) \
 $(LIB)(card6.o) \
 $(LIB)(card7.o) \
 $(LIB)(cards8and9.o) \
 $(LIB)(ccoef.o) \
 $(LIB)(cc3.o) \
 $(LIB)(cc3m.o) \
 $(LIB)(cc3_c.o) \
 $(LIB)(cc3m_c.o) \
 $(LIB)(ccp.o) \
 $(LIB)(checksym.o) \
 $(LIB)(crop_ft.o) \
 $(LIB)(ctfapply.o) \
 $(LIB)(ctfapply_phase_only.o) \
 $(LIB)(ctf.o) \
 $(LIB)(ctfrefine.o) \
 $(LIB)(combine_arrays.o) \
 $(LIB)(correct3d.o) \
 $(LIB)(correct3d_c.o) \
 $(LIB)(d2mask.o) \
 $(LIB)(d3mask.o) \
 $(LIB)(dump.o) \
 $(LIB)(erc.o) \
 $(LIB)(ewaldex.o) \
 $(LIB)(ewaldin.o) \
 $(LIB)(fd.o) \
 $(LIB)(fsh.o) \
 $(LIB)(fangle.o) \
 $(LIB)(find_fpart.o) \
 $(LIB)(find_rlim.o) \
 $(LIB)(fpart_cc.o) \
 $(LIB)(fftw_plans.o) \
 $(LIB)(getsymmat.o) \
 $(LIB)(hextend.o) \
 $(LIB)(iof.o) \
 $(LIB)(limitsymm.o) \
 $(LIB)(magrefine.o) \
 $(LIB)(lmain.o) \
 $(LIB)(maskcos.o) \
 $(LIB)(maskcos2.o) \
 $(LIB)(maskcos_c.o) \
 $(LIB)(maskcos2_c.o) \
 $(LIB)(maskenv.o) \
 $(LIB)(mask3d.o) \
 $(LIB)(mask3d_c.o) \
 $(LIB)(match.o) \
 $(LIB)(matmul.o) \
 $(LIB)(matmul_t.o) \
 $(LIB)(opmaps.o) \
 $(LIB)(opmaps2.o) \
 $(LIB)(opresstathalf.o) \
 $(LIB)(opresstatmaps.o) \
 $(LIB)(pad.o) \
 $(LIB)(pdiff.o) \
 $(LIB)(pextract.o) \
 $(LIB)(phaseflip.o) \
 $(LIB)(pinc.o) \
 $(LIB)(pinsert.o) \
 $(LIB)(pinsert_c.o) \
 $(LIB)(pinsert_nearest.o) \
 $(LIB)(pinsert_trilin.o) \
 $(LIB)(pinsert_nearestb.o) \
 $(LIB)(pinsert_trilinb.o) \
 $(LIB)(pinsert_nearest_se.o) \
 $(LIB)(pinsert_trilin_se.o) \
 $(LIB)(pinsert_nearestb_se.o) \
 $(LIB)(pinsert_trilinb_se.o) \
 $(LIB)(prefine.o) \
 $(LIB)(presb.o) \
 $(LIB)(psearch.o) \
 $(LIB)(pweight.o) \
 $(LIB)(random.o) \
 $(LIB)(remap_theta.o) \
 $(LIB)(rotate2d.o) \
 $(LIB)(rotmat.o) \
 $(LIB)(searchang.o) \
 $(LIB)(sheltest.o) \
 $(LIB)(shift2D.o) \
 $(LIB)(shift.o) \
 $(LIB)(shiftvol.o) \
 $(LIB)(sigma.o) \
 $(LIB)(slen2.o) \
 $(LIB)(stamp.o) \
 $(LIB)(storeshift.o) \
 $(LIB)(trilinmap.o) \
 $(LIB)(va04.o) \
 $(LIB)(variance.o) \
 $(LIB)(volmeasure.o) \
 $(LIB)(window.o)
ifeq ($(COMP),ifort)
ELEMENTS += \
 $(LIB)(gcvspl.o) \
 $(LIB)(sort2.o) \
 $(LIB)(dates_and_times.o) \
 $(LIB)(program_instances.o) \
 $(LIB)(globals.o) \
 $(LIB)(look_up_tables.o) \
 $(LIB)(useful_functions.o) \
 $(LIB)(string_manipulations.o) \
 $(LIB)(user_supplied_parameters.o) \
 $(LIB)(user_input.o) \
 $(LIB)(histograms.o) \
 $(LIB)(empirical_distributions.o) \
 $(LIB)(peaks.o) \
 $(LIB)(image_header_records.o) \
 $(LIB)(image_headers.o) \
 $(LIB)(image_files.o) \
 $(LIB)(units.o) \
 $(LIB)(progress_bars.o) \
 $(LIB)(numeric_text_files.o) \
 $(LIB)(curves.o) \
 $(LIB)(contrast_transfer_functions.o) \
 $(LIB)(images_core.o)
endif

ifeq ($(MKLROOT),"")
 FFTLIB  = libfftw3f.a
else
 FFTLIB  =
endif

ifeq ($(LIB),$(PGM)_mp.a)
$(PGM) : $(PGM).o $(LIB) resample_mp merge_3d_mp
	$(COMP) $(PGM).o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/$(PGM)_mp.exe
	@\rm $(PGM).o
	@\rm $(PGM).a
	@\rm fftw3.3.o
	@\rm fftw33.mod
else
$(PGM) : $(PGM).o $(LIB) resample fheader merge_classes select_classes merge_3d calc_occ calc_occ_helical set_polarity rsample bootvar bfactor apply_mask
	$(COMP) $(PGM).o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/$(PGM).exe
	@\rm $(PGM).o
	@\rm $(PGM).a
	@\rm fftw3.3.o
	@\rm fftw33.mod
endif

ifeq ($(COMP),ifort)
resample : resample.o $(LIB)
	$(COMP) resample.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/resample.exe
	@\rm resample.o

resample_mp : resample.o $(LIB)
	$(COMP) resample.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/resample_mp.exe
	@\rm resample.o

resample.o  : resample.f90
	$(COMP) $(FFLAGS) -o resample.o  resample.f90
else
resample :
resample_mp :
endif

fheader : fheader.o $(LIB)
	$(COMP) fheader.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/fheader.exe
	@\rm fheader.o

merge_classes : merge_classes.o $(LIB)
	$(COMP) merge_classes.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/merge_classes.exe
	@\rm merge_classes.o

select_classes : select_classes.o $(LIB)
	$(COMP) select_classes.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/select_classes.exe
	@\rm select_classes.o

merge_3d : merge_3d.o $(LIB)
	$(COMP) merge_3d.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/merge_3d.exe
	@\rm merge_3d.o

merge_3d_mp : merge_3d.o $(LIB)
	$(COMP) merge_3d.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/merge_3d_mp.exe
	@\rm merge_3d.o

calc_occ : calc_occ.o
	$(COMP) calc_occ.o $(FFTLIB) $(LF) $(OSX) -o ../bin/calc_occ.exe
	@\rm calc_occ.o

calc_occ_helical : calc_occ_helical.o
	$(COMP) calc_occ_helical.o $(FFTLIB) $(LF) $(OSX) -o ../bin/calc_occ_helical.exe
	@\rm calc_occ_helical.o

set_polarity : set_polarity.o
	$(COMP) set_polarity.o $(FFTLIB) $(LF) $(OSX) -o ../bin/set_polarity.exe
	@\rm set_polarity.o

rsample : rsample.o
	$(COMP) rsample.o $(FFTLIB) $(LF) $(OSX) -o ../bin/rsample.exe
	@\rm rsample.o

bootvar : bootvar.o $(LIB)
	$(COMP) bootvar.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/bootvar.exe
	@\rm bootvar.o

bfactor : bfactor.o $(LIB)
	$(COMP) bfactor.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/bfactor.exe
	@\rm bfactor.o

apply_mask : apply_mask.o $(LIB)
	$(COMP) apply_mask.o $(LIB) $(FFTLIB) $(LF) $(OSX) -o ../bin/apply_mask.exe
	@\rm apply_mask.o

$(PGM).o  : $(PGM).f
	$(COMP) $(FFLAGS) fftw3.3.f90
	$(COMP) $(FFLAGS) -o $(PGM).o  $(PGM).f
	$(AR) r $(LIB) fftw3.3.o
	

fheader.o  : fheader.f
	$(COMP) $(FFLAGS) -o fheader.o  fheader.f

merge_classes.o  : merge_classes.f
	$(COMP) $(FFLAGS) -o merge_classes.o  merge_classes.f

select_classes.o  : select_classes.f
	$(COMP) $(FFLAGS) -o select_classes.o  select_classes.f

merge_3d.o  : merge_3d.f
	$(COMP) $(FFLAGS) -o merge_3d.o  merge_3d.f

calc_occ.o  : calc_occ.f
	$(COMP) $(FFLAGS) -o calc_occ.o  calc_occ.f

calc_occ_helical.o  : calc_occ_helical.f
	$(COMP) $(FFLAGS) -o calc_occ_helical.o  calc_occ_helical.f

set_polarity.o  : set_polarity.f
	$(COMP) $(FFLAGS) -o set_polarity.o  set_polarity.f

rsample.o  : rsample.f
	$(COMP) $(FFLAGS) -o rsample.o  rsample.f

bootvar.o  : bootvar.f
	$(COMP) $(FFLAGS) -o bootvar.o  bootvar.f

bfactor.o  : bfactor.f
	$(COMP) $(FFLAGS) -o bfactor.o  bfactor.f

apply_mask.o  : apply_mask.f
	$(COMP) $(FFLAGS) -o apply_mask.o  apply_mask.f
