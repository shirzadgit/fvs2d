include make.inc

fvs2d:  setup lib nlin

clean: cleancode cleanlib


######################################
setup:
	(./make.sh)


######################################
nlin:
	(sed -i '/SOLVER_TYPE/c\SOLVER_TYPE=NLIN' "make.inc";	cd src; $(MAKE) solver)


######################################
lib:
	( cd libs; $(MAKE) lib )


######################################
cleancode:
	( cd src; $(MAKE) clean )


cleanlib:
	( cd libs; $(MAKE) clean )
