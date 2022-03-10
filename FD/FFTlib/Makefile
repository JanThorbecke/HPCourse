# Makefile for ...FFT

ROOT = ../

B = $(ROOT)/bin
I = $(ROOT)/include
L = $(ROOT)/lib
M = $(ROOT)/man

include ../Make_include

LIB = $L/libgenfft.a

CSRC =			\
	cc1fft.c	\
	ccmfft.c	\
	cc2dfft.c	\
	rc1fft.c	\
	cr1fft.c	\
	rcmfft.c	\
	crmfft.c	\
	rc2dfft.c	\
	cr2dfft.c	\
	xt2wx.c	\
	xt2wkx.c	\
	wkx2xt.c	\
	wx2xt.c	\
	pfafft.c \
	fft_mayer.c \
	wallclock_time.c \
	optnumber.c \
	lib_fft.c

OBJC    = $(CSRC:%.c=%.o)
OBJF    = $(FSRC:%.f=%.o)

CARCH = $(patsubst %.c, $(LIB)(%.o), $(CSRC))
FARCH = $(patsubst %.f, $(LIB)(%.o), $(FSRC))

ARCH := $(CARCH) $(FARCH)

INSTALL: mkdirs $(LIB)

$(LIB)	: $(ARCH)
	$(RANLIB) $(LIB)

mkdirs:
	-mkdir $(ROOT)/lib
	-mkdir $(ROOT)/include
	-ln -sf $(PWD)/genfft.h $I/genfft.h 

bins: 
	cd test ; $(MAKE)

remake	:
	@-rm -f $(LIB) INSTALL
	@$(MAKE)

list	:
	$(AR) tv $(LIB)

clean:
	-rm -f core a.out $(LIB) *.vect $I/genfft.h

#	@cd test ; $(MAKE) clean

realclean:
	$(MAKE) clean
	-rm -f $(LIB) 
#	@cd test ; $(MAKE) realclean 

.PRECIOUS:	 $(LIB)


