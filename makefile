CLANG = clang
GCC = gcc
ICC = icc
CC = $(CLANG)

# IMPORTANT steps for choosing and configuring the RNG used!
# Please provide EXACTLY ONE of
# -DuseSFMT -DDSFMT_MEXP=19937
# -DuseMRG
# -DuseDRAND
PARAMS = -DDSFMT_MEXP=19937 -DuseSFMT

SOURCES = phi_quad.c dSFMT/dSFMT.c

# flags
DEBUG = -g -fsanitize=address,undefined,bounds,null -fno-omit-frame-pointer
DEBUG_GCC = -g -fsanitize=address -fno-omit-frame-pointer
LIBS = -lm
LIBS_ICC = -lm
RELEASE = -Ofast -march="native" -Rpass-missed -Rpass-analysis -ffast-math
RELEASE_GCC = -Ofast -march="native" #-fopt-info-missed
PROFILE =  -Ofast -march="native" -fprofile-instr-generate
PROFILE_GCC =  -O1 -march="native" -pg 
WARN = -Weverything -Wno-sign-conversion -Wno-pedantic -Wno-padded
WARN_GCC = -Wall -Wextra
MODE = 
MODE_GCC = -std=c11
MODE_ICC = -std=gnu11

all: clang
clang: clang-debug clang-release
gcc: gcc-debug gcc-release
icc: icc-release

CLANG_COMMON = $(CC) $(PARAMS) $(MODE) $(WARN) $(LIBS) 

clang-release: $(SOURCES)
	$(CLANG_COMMON) $(RELEASE) -o clang_release $(SOURCES)

clang-profile: $(SOURCES)
	$(CLANG_COMMON) $(PROFILE) -o clang_profile $(SOURCES)

clang-debug: $(SOURCES)
	$(CLANG_COMMON) $(DEBUG) -o clang_debug $(SOURCES)


GCC_COMMON = $(GCC) $(PARAMS) $(MODE_GCC) $(WARN_GCC)

gcc-release: $(SOURCES)
	$(GCC_COMMON) $(RELEASE_GCC) -o gcc_release $(SOURCES) $(LIBS)

gcc-profile: $(SOURCES)
	$(GCC_COMMON) $(PROFILE_GCC) -o gcc_profile $(SOURCES) $(LIBS)

gcc-debug: $(SOURCES)
	$(GCC_COMMON) $(DEBUG_GCC) -o gcc_debug $(SOURCES) $(LIBS)


icc-release: $(SOURCES)
	$(ICC) $(PARAMS) $(MODE_ICC) $(WARN_GCC) $(RELEASE_GCC) -o icc_release $(SOURCES) $(LIBS_ICC) 