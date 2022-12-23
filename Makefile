PNAME         = imthreshold
PROGNAME      = $(PNAME)-deskew $(PNAME)-rotate $(PNAME)-filter $(PNAME)-fautoinv $(PNAME)-fdespeckle $(PNAME)-finfo $(PNAME)-fpmean $(PNAME)-math $(PNAME)-separate $(PNAME)-size $(PNAME)-sbwmag2 $(PNAME)-shris $(PNAME)-tglobal $(PNAME)-tlocal $(PNAME)-tlayer $(PNAME)-tdjvul $(PNAME)-tgatos $(PNAME)-thalftone2 $(PNAME)-ttext $(PNAME)-twhiterohrer
SRCS          = src
CC            = gcc
CPP           = g++
CFLAGS        = -DUNIX -I$(SRCS) -Wall -s
VER           = 0
VERB          = 20221223
ifeq ($(OS),Windows_NT)
LIBS          = $(SRCS)/FreeImage.lib
PLIBF         = $(PNAME).$(VER).dll
PLIBFI        = $(PNAME)freeimage.$(VER).dll
RM            = del /Q
else
LIBS          = -lfreeimage
PLIBF         = lib$(PNAME).so.$(VER)
PLIBFI        = lib$(PNAME)freeimage.so.$(VER)
RM            = rm -f
endif
PLIB          = $(PLIBF) $(PLIBFI)
PREFIX        = /usr/local
INSTALL       = install
LN            = ln -fs

.PHONY: all clean install

all: $(PROGNAME)

clean:
	$(RM) $(PROGNAME) $(PLIB) *.exe

$(PLIBF): $(SRCS)/lib/color.c $(SRCS)/lib/commom.c $(SRCS)/lib/denoise.c $(SRCS)/lib/despeckle.c $(SRCS)/lib/filter.c $(SRCS)/lib/math.c $(SRCS)/lib/pmean.c $(SRCS)/lib/rotate.c $(SRCS)/lib/separate.c $(SRCS)/lib/size.c $(SRCS)/lib/threshold.c
	$(CC) $(CFLAGS) -shared -Wl,-soname,$@ $^ -o $@

$(PLIBFI): $(SRCS)/imthresholdfreeimage.cpp $(PLIBF)
	$(CPP) $(CFLAGS) -shared -Wl,-soname,$@ $^ $(LIBS) -o $@

$(PNAME)-deskew: $(SRCS)/deskew.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-rotate: $(SRCS)/rotate.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-filter: $(SRCS)/filter.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-fautoinv: $(SRCS)/fautoinv.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-fdespeckle: $(SRCS)/fdespeckle.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-finfo: $(SRCS)/finfo.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-fpmean: $(SRCS)/fpmean.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-math: $(SRCS)/math.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-separate: $(SRCS)/separate.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-size: $(SRCS)/size.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-sbwmag2: $(SRCS)/sbwmag2.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-shris: $(SRCS)/shris.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-tglobal: $(SRCS)/tglobal.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-tlocal: $(SRCS)/tlocal.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-tlayer: $(SRCS)/tlayer.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-tdjvul: $(SRCS)/tdjvul.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-tgatos: $(SRCS)/tgatos.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-thalftone2: $(SRCS)/thalftone2.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-ttext: $(SRCS)/ttext.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

$(PNAME)-twhiterohrer: $(SRCS)/twhiterohrer.cpp $(PLIB)
	$(CPP) $(CFLAGS) $^ $(LIBS) -o $@

install: $(PROGNAME)
	$(INSTALL) -d $(PREFIX)/lib
	$(INSTALL) -m 0644 lib$(PNAME).so.$(VER) $(PREFIX)/lib/lib$(PNAME).so.$(VER).$(VERB)
	$(LN)  lib$(PNAME).so.$(VER).$(VERB) $(PREFIX)/lib/lib$(PNAME).so.$(VER)
	$(LN)  lib$(PNAME).so.$(VER) $(PREFIX)/lib/lib$(PNAME).so
	$(INSTALL) -m 0644 lib$(PNAME)freeimage.so.$(VER) $(PREFIX)/lib/lib$(PNAME)freeimage.so.$(VER).$(VERB)
	$(LN)  lib$(PNAME)freeimage.so.$(VER).$(VERB)  $(PREFIX)/lib/lib$(PNAME)freeimage.so.$(VER)
	$(LN)  lib$(PNAME)freeimage.so.$(VER)  $(PREFIX)/lib/lib$(PNAME)freeimage.so
	$(INSTALL) -d $(PREFIX)/include
	$(INSTALL) -m 0644 src/imthreshold*.h $(PREFIX)/include
	$(INSTALL) -d $(PREFIX)/bin
	$(INSTALL) -m 0755 $(PROGNAME) $(PREFIX)/bin/
	$(INSTALL) -d $(PREFIX)/share/man/man1
	$(INSTALL) -m 0644 man/man1/imthreshold*.1 $(PREFIX)/share/man/man1
	$(INSTALL) -d $(PREFIX)/share/doc/$(PNAME)
	$(INSTALL) -m 0644 README.md doc/imthreshold/* $(PREFIX)/share/doc/$(PNAME)
