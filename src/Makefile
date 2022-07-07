include ../make.inc

#----------------------------------------------------------------------------------------
# General
#----------------------------------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .f .F .f90 .F90 .o

.f.o:
	$(FCOMPL) -c $(FFLAGC) $(FFDBG) $(FFLAGE) -o $*.o $<

.F.o:
	$(FCOMPL) -c $(FFLAGC) $(FFDBG) $(FFLAGE) -o $*.o $<

.f90.o:
	$(FCOMPL) -c $(FFLAGC) $(FFDBG) $(FFLAGE) -o $*.o $<

.F90.o:
	$(FCOMPL) -c $(FFLAGC) $(FFDBG) $(FFLAGE) -o $*.o $<

.s.o:
	as $<


#----------------------------------------------------------------------------------------
# Source Files
#----------------------------------------------------------------------------------------
OBJS = fvs2d.o
MODOBJS =	mainparam.o input.o kdtree2.o \
					data_grid.o data_solution.o \
					grid_procs.o \
					interpolation.o \
					gradient_ggcb.o gradient_ggnb.o gradient_lsq.o gradient.o gradient_limiter.o \
					runge_kutta.o \
					mms.o \
					initialize.o \
					flux_invscid.o residual.o \
					tecplot.o \
					io.o \
					test.o


#----------------------------------------------------------------------------------------
# FLAGS for Code
#----------------------------------------------------------------------------------------
FT =
ext=
CODE="fvs2d_nlin"


#----------------------------------------------------------------------------------------
# Debug option
#----------------------------------------------------------------------------------------
ifneq (,$(findstring DBG,$(DEBUG)))
FFDBG = $(FFDEBG)
CODEMPI=$(CODE)"_debug.x"
else
FFDBG =
CODEMPI=$(CODE)".x"
endif


#----------------------------------------------------------------------------------------
# Generate code
#----------------------------------------------------------------------------------------
FF = $"FFLAGE=$(FT)
solver:
	$(MAKE) $(FF) \
	$(MODOBJS) $(OBJS)
	@echo " "
	@echo "======================================"
	@echo " Linking the target at run/: $(CODEMPI)"
	@echo " host: $(host)"
	@echo "======================================"
	@echo " "
	$(FCOMPL) $(FFLAGC) $(FFDBG) $(FOMP) -o $(RUN)/$(CODEMPI) \
	$(OBJS) $(MODOBJS) $(LIBS)
	chmod 755 $(RUN)/$(CODEMPI)


#----------------------------------------------------------------------------------------
# Clean
#----------------------------------------------------------------------------------------
cleanob:
	@rm -f *.o *.mod

clean:
	@rm -f *.o *.mod $(RUN)/*.x
