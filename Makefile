CXX = clang++-4.0

PROFILE = false
DEBUG = true
OPT_LVL = 3

DFLAGS = POPULATION=50 ELITENESS=2 TERMINATE_TIME=60 # DYNAMIC_FITNESS
DFLAGS_ = $(foreach flag,$(DFLAGS),-D$(flag))
CFLAGS = -std=c++1z -O$(OPT_LVL) -Wall -Wextra -Wno-char-subscripts -Wno-unused-result -I./src/include $(DFLAGS_) -march=native -mfpmath=sse
LDFLAGS = -pthread

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
