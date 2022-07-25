include ../../make.inc
include ../make.inc

OBJ =	dgemm.o dger.o dgesv.o dgetf2.o dgetrf.o dgetrs.o \
	dlamch.o dlaswp.o dscal.o dswap.o dtrsm.o idamax.o \
	ieeeck.o ilaenv.o iparmq.o lsame.o xerbla.o

.f.o:
	$(FCOMPL) -c $(FFLAGL) $<

$(BLASLIB): $(OBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(OBJ)
	$(RANLIB) $@

all: $(BLASLIB)

clean:
	rm -f $(OBJ) $(BLASLIB)
