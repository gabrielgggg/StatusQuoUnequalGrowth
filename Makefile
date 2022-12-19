#
# Set compiler. Valid options are:
#   gnu    for gfortran
#   cray   for Cray ftn
#   intel  for Intel ifort
#   arm    for ARM LLVM (untested)
#
tools = intel

ifeq ($(tools),gnu)
  compiler = gfortran
  flags = -Ofast -march=native -mtune=native -pipe -fopenmp -flto
  fWarn = -Wall -pedantic
  libs = -llapack -lblas -lm
else ifeq ($(tools),cray)
  compiler = ftn
  flags = -O3 -h omp,scalar3,vector3,ipa3,contiguous,nobounds
  fWarn = 
  libs =
else ifeq ($(tools),intel)
  compiler = ifort
  flags = -O3 -parallel -ipo -xHost -align -qopenmp -std18
  fWarn = -warn all,noexternal
  libs = -qmkl
else ifeq ($(tools),arm)
  compiler = armflang
  flags = -Ofast -fopenmp -mcpu=a64fx -armpl -flto
  fWarn = -Wall
  libs = 
else
  $(error Unknown toolchain?)
endif

.PHONY: all clean

all:	unequal.f90
	$(compiler) $(flags) $(fWarn) -o unequal unequal.f90 ${libs}

clean:
	rm -f ./unequal ./*.mod ./*.o

cleanresults:
	rm -f results/*.tab results/*.bin

cleanmlab:
	rm -f results/*.mat results/*.fig

ccc:	clean cleanresults

