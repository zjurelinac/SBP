CXX = g++#clang++-4.0 # g++

PROFILE = true
DEBUG = true
OPT_LVL = 3

DFLAGS =
DFLAGS_ = $(foreach flag,$(DFLAGS),-D$(flag))
CFLAGS = -std=c++14 -O$(OPT_LVL) -msse4.2 -Wall -Wextra -Wno-char-subscripts -I./src/include $(DFLAGS_)

ifeq ($(DEBUG), true)
	CFLAGS += -g -DDEBUG
endif

ifeq ($(PROFILE), true)
	CFLAGS += -pg
	LDFLAGS += -pg
endif

SDIR = src
IDIR = src/include
ODIR = build

EXEC = build/sbr
SRCS = $(wildcard $(SDIR)/*.cpp)
DEPS = $(wildcard $(IDIR)/*.h $(IDIR)/*.hpp)
OBJS = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRCS))

all: .build $(EXEC)

.build:
	@mkdir -p build

$(EXEC): $(OBJS) $(DEPS)
	@echo $(CXX): building $@...
	@$(CXX) -o $@ $(OBJS) $(LDFLAGS)
	@echo Completed.

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	@echo $(CXX): compiling $<...
	@$(CXX) -c -o $@ $< $(CFLAGS)

clean:
	@\rm -r build
