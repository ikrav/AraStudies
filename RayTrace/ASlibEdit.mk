#**** Shared Library Generation Makefile
#**** Home Version
#**** AJS Original Made on December 27, 2014 (CsvCreator)

#**** Control Variables

DEPENDANCIES = ASlibEdit
#**** ABSOLUTELY no spaces after LIBRARY_CLASS equivalence
LIBRARY_CLASS = ASlibEdit

#ROOTCXXFLAGS = -O2 -fPIC -pthread -m64 -I/home/andrew/root/include
#ROOTLIBFLAGS = -L/home/andrew/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic

ROOTCXXFLAGS = -O2 -fPIC -pthread -m64 -I/net/ara/repositories/root_5.34/include
ROOTLIBFLAGS = -L/net/ara/repositories/root_5.34/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic

#**** End Control Variables

OBJ_SUFF = o
SO_SUFF = so
HEAD_SUFF = h
DEF_SUFF = cxx

SHAREDLIBNAME = lib$(LIBRARY_CLASS)
SHAREDLIB = $(SHAREDLIBNAME).so

OBJS = $(addsuffix .$(OBJ_SUFF), $(DEPENDANCIES)) $(SHAREDLIBNAME)Dict.o
DEFINITIONS = $(addsuffix .$(DEF_SUFF), $(DEPENDANCIES))
HEADERS = $(addsuffix .$(HEAD_SUFF), $(DEPENDANCIES))

COMPILE = g++ # -std=c++11
MAKE_SHARED_LIBRARY = -shared
MAKE_DICTIONARY = rootcint -f

all : $(SHAREDLIB)

$(SHAREDLIBNAME)Dict.$(DEF_SUFF) : $(HEADERS)
	@echo "**Generating " $(SHAREDLIBNAME)"Dict."$(DEF_SUFF) " dictionary ..."
	$(MAKE_DICTIONARY) $(SHAREDLIBNAME)Dict.$(DEF_SUFF) -c $(ROOTCXXFLAGS) -p $(HEADERS) LinkDef.h

%.$(OBJ_SUFF) : %.$(DEF_SUFF)
	@echo "**Compiling " $^ " to generate " $@ " ..."
	$(COMPILE) $(ROOTCXXFLAGS) -c $< -o $@

$(SHAREDLIB) : $(OBJS)
	@echo "**Finishing up... " $@
	$(COMPILE) $(MAKE_SHARED_LIBRARY) -o $(SHAREDLIB) $(ROOTLIBFLAGS) $(OBJS)
	@echo "**Done!!!"

clean:
	rm *.$(OBJ_SUFF)
	rm *Dict*

