EXEDIR := run
OBJDIR := bin
SRCDIR := src
INCDIR := inc
MAKEDIR := bin
LIBFILE := $(OBJDIR)/libStatObj.a

CXX := $(shell root-config --cxx)
EXTRA_WARNINGS := -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Winit-self -Winvalid-pch -Wlong-long -Wmissing-format-attribute -Wmissing-include-dirs -Wmissing-noreturn -Wpacked -Wpointer-arith -Wredundant-decls -Wstack-protector -Wswitch-default -Wswitch-enum -Wundef -Wunused -Wvariadic-macros -Wwrite-strings -Wabi -Wctor-dtor-privacy -Wnon-virtual-dtor -Wsign-promo -Wsign-compare #-Wunsafe-loop-optimizations -Wfloat-equal -Wsign-conversion -Wunreachable-code
CXXFLAGS := -isystem $(shell root-config --incdir) -Wall -Wextra -pedantic -Werror -Wshadow -Woverloaded-virtual -Wold-style-cast $(EXTRA_WARNINGS) $(shell root-config --cflags) -O2 -I $(INCDIR)
LD := $(shell root-config --ld)
LDFLAGS := $(shell root-config --ldflags)
LDLIBS := $(shell root-config --libs) -lMinuit -lRooStats -lTreePlayer 

EXECUTABLES := $(addprefix $(EXEDIR)/, $(addsuffix .exe, $(notdir $(basename $(wildcard $(SRCDIR)/*.cxx))))) 
OBJECTS := $(addprefix $(OBJDIR)/, $(addsuffix .o, $(notdir $(basename $(wildcard $(SRCDIR)/*.cc)))))

FIND_DEPS = $(CXX) $(CXXFLAGS) -MM -MG -MF $@ $<
EXPAND_DEPS = perl -pi -e 's|$*.o|$(OBJDIR)/$*.o $(MAKEDIR)/$*.d|g' $@
GET_DEPS = $(FIND_DEPS) && $(EXPAND_DEPS)
COMPILE = $(CXX) $(CXXFLAGS) -o $@ -c $<
LINK = $(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

vpath %.cc $(SRCDIR)
vpath %.cxx $(SRCDIR)
vpath %.hh $(INCDIR)
vpath %.o $(OBJDIR)
vpath %.exe $(EXEDIR)
vpath %.d $(MAKEDIR)

all: $(EXECUTABLES) $(CMSSW_BASE)/src/babymaker/bmaker/src/baby_base.cc $(CMSSW_BASE)/src/babymaker/bmaker/interface/baby_base.hh

-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cc)))))
-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cxx)))))
-include $(MAKEDIR)/baby.d 

$(LIBFILE): $(OBJECTS)

$(MAKEDIR)/%.d: $(SRCDIR)/%.cc
	$(GET_DEPS)

$(MAKEDIR)/%.d: $(SRCDIR)/%.cxx
	$(GET_DEPS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(COMPILE)

$(OBJDIR)/%.o: $(SRCDIR)/%.cxx
	$(COMPILE)

$(OBJDIR)/%.a:
	ar rcsv $@ $^

$(EXEDIR)/generate_baby.exe: $(OBJDIR)/generate_baby.o
	$(LINK)

$(EXEDIR)/%.exe: $(OBJDIR)/%.o $(LIBFILE)
	$(LINK)

# Auto-generated code
.SECONDARY: dummy_baby.all 
.PRECIOUS: generate_baby.o 

$(CMSSW_BASE)/src/babymaker/bmaker/src/baby_base.cc $(CMSSW_BASE)/src/babymaker/bmaker/interface/baby_base.hh: dummy_baby.all
dummy_baby.all: $(EXEDIR)/generate_baby.exe $(shell find $(CMSSW_BASE)/src/babymaker/variables -type f)
	./$<

.DELETE_ON_ERROR:
