DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

# Source files
ALL_SRC := $(wildcard ${DIR_SRC}/*.cpp) $(wildcard ${DIR_SRC}/bmdscore/*.cpp)

# Compiler and flags
CXX ?= g++

# Release flags
CXXFLAGS_RELEASE := -std=c++11 -pthread -g -O3 -MD -MP -I${DIR_INC} \
		-I/usr/local/include/nlopt \
		-I/usr/include/eigen3 \
		-I./src/zlib \
		-I./src/include \
		-I./src/bmdscore \
		$(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))

# Debug flags
CXXFLAGS_DEBUG := -std=c++11 -pthread -g -O0 -MD -MP -fno-omit-frame-pointer -fsanitize=address -I${DIR_INC} \
		-I/usr/local/include/nlopt \
		-I/usr/include/eigen3 \
		-I./src/zlib \
		-I./src/include \
		-I./src/bmdscore \
		$(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))

LIBS := -lnlopt -lgsl -lgslcblas -lisal -ldeflate -lpthread -lz
STATIC_FLAGS := -static -Wl,--no-as-needed -pthread

LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)
STATIC_LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(STATIC_FLAGS) $(LIBS)

# Select files for each executable
BMDSEQ_SRC := $(filter-out ${DIR_SRC}/kmer2abunt.cpp, ${ALL_SRC})
KMER2ABUNT_SRC := ${DIR_SRC}/options.cpp \
				  ${DIR_SRC}/kmer2abunt.cpp \
				  ${DIR_SRC}/bmdseeker.cpp \
				  ${DIR_SRC}/fastareader.cpp \
				  ${DIR_SRC}/kmer.cpp \
				  ${DIR_SRC}/read.cpp \
				  ${DIR_SRC}/writer.cpp \
				  ${DIR_SRC}/sequence.cpp \
				  $(wildcard ${DIR_SRC}/bmdscore/*.cpp)

# Object files per target
BMDSEQ_OBJ := $(patsubst ${DIR_SRC}/%.cpp,${DIR_OBJ}/%.o,${BMDSEQ_SRC})
KMER2ABUNT_OBJ := $(patsubst ${DIR_SRC}/%.cpp,${DIR_OBJ}/%.o,${KMER2ABUNT_SRC})

# Targets
TARGETS := bmdseq kmer2abunt

.PHONY: all debug clean static install

all: CXXFLAGS := $(CXXFLAGS_RELEASE)
all: LD_FLAGS := $(LD_FLAGS)
all: ${TARGETS}

debug: CXXFLAGS := $(CXXFLAGS_DEBUG)
debug: LD_FLAGS := $(LD_FLAGS) -fsanitize=address
debug: clean
debug: ${TARGETS}

bmdseq: ${BMDSEQ_OBJ}
	$(CXX) $^ -o $@ $(LD_FLAGS)

kmer2abunt: ${KMER2ABUNT_OBJ}
	$(CXX) $^ -o $@ $(LD_FLAGS)

static:
	$(CXX) $(BMDSEQ_OBJ) -o bmdseq $(STATIC_LD_FLAGS)
	$(CXX) $(KMER2ABUNT_OBJ) -o kmer2abunt $(STATIC_LD_FLAGS)

# Compile rules
${DIR_OBJ}/%.o: ${DIR_SRC}/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(DIR_OBJ)/bmdscore/%.o: ${DIR_SRC}/bmdscore/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	@rm -rf $(DIR_OBJ)
	@rm -f ${TARGETS}

install:
	install bmdseq $(BINDIR)/bmdseq
	install kmer2abunt $(BINDIR)/kmer2abunt
	@echo "Installed."

-include $(BMDSEQ_OBJ:.o=.d) $(KMER2ABUNT_OBJ:.o=.d)
