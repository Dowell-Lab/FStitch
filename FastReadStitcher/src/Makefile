#***************************************************************************/
#**    Author: Joey Azofeifa                                               */
#**    Date:   1/2/2015	                                                   */
#**    Description: Makefile for FastReadStitcher                          */
#***************************************************************************/

CXX             = g++
CXXFLAGS        = -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11  -fopenmp -Wno-write-strings -g
EXEC            = ${PWD}/FStitch
UNITTEST        = ${PWD}/unittest
ARCH  			= getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 		= $(shell gcc --version | grep ^gcc | sed 's/^.* //g')
FStitch: main.o BaumWelch.o NewtonsMethod.o read_in_parameters.o read.o main_train.o grabTrainingExamples.o interval_tree.o \
		 split.o write.o main_segment.o viterbi.o validate.o ParamWrapper.o unittests.o
	@${CXX} ${CXXFLAGS} ${PWD}/main.o  ${PWD}/BaumWelch.o \
		${PWD}/NewtonsMethod.o ${PWD}/read_in_parameters.o ${PWD}/read.o \
		${PWD}/interval_tree.o ${PWD}/grabTrainingExamples.o ${PWD}/main_train.o \
		${PWD}/split.o ${PWD}/write.o ${PWD}/main_segment.o ${PWD}/viterbi.o ${PWD}/validate.o ${PWD}/ParamWrapper.o -o ${EXEC}
	@echo "========================================="
	@echo Sucessfully Compiled FStitch executable.
	@echo Building unit test executable...
	@${CXX} ${CXXFLAGS} ${PWD}/unittests.o  ${PWD}/BaumWelch.o \
		${PWD}/NewtonsMethod.o ${PWD}/read_in_parameters.o ${PWD}/read.o \
		${PWD}/interval_tree.o ${PWD}/grabTrainingExamples.o ${PWD}/main_train.o \
		${PWD}/split.o ${PWD}/write.o ${PWD}/main_segment.o ${PWD}/viterbi.o ${PWD}/validate.o ${PWD}/ParamWrapper.o -o ${UNITTEST}
	@echo Done.
	
		
main.o: main.cpp
		@echo "========================="
		
		@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
		@echo Compiled main
read_in_parameters.o: read_in_parameters.cpp
		@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp 
		@echo Compiled read_in_parameters
read.o: read.cpp		
		@@${CXX} -c ${CXXFLAGS} ${PWD}/read.cpp 
		@echo Compiled read
main_train.o: main_train.cpp		
		@${CXX} -c ${CXXFLAGS} ${PWD}/main_train.cpp 
		@echo Compiled main_train
split.o: split.cpp
		@@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp
		@echo Compiled split
grabTrainingExamples.o: grabTrainingExamples.cpp 
		@${CXX} -c ${CXXFLAGS} ${PWD}/grabTrainingExamples.cpp 
		@echo Compiled grabTrainingExamples
interval_tree.o: interval_tree.cpp
		@${CXX} -c ${CXXFLAGS} ${PWD}/interval_tree.cpp 
		@echo Compiled interval_tree
NewtonsMethod.o: NewtonsMethod.cpp
		@${CXX} -c ${CXXFLAGS} ${PWD}/NewtonsMethod.cpp 
		@echo Compiled NewtonsMethod
BaumWelch.o: BaumWelch.cpp
		@${CXX} -c ${CXXFLAGS} ${PWD}/BaumWelch.cpp 
		@echo Compiled BaumWelch
main_segment.o: main_segment.cpp
		@${CXX} -c ${CXXFLAGS} ${PWD}/main_segment.cpp 
		@echo Compiled main_segment
viterbi.o: viterbi.cpp
		@${CXX} -c ${CXXFLAGS} ${PWD}/viterbi.cpp 
		@echo Compiled viterbi
validate.o: validate.cpp
		@${CXX} -c ${CXXFLAGS} ${PWD}/validate.cpp 
		@echo Compiled validate
		
write.o: write.cpp
	@${CXX} -c ${CXXFLAGS} ${PWD}/write.cpp 
	@echo Compiled write

ParamWrapper.o: ParamWrapper.cpp
	@${CXX} -c ${CXXFLAGS} ${PWD}/ParamWrapper.cpp 
	@echo Compiled ParamWrapper

unittests.o: unittests.cpp
	@${CXX} -c ${CXXFLAGS} ${PWD}/unittests.cpp
	@echo Compiled unittests

unittests: unittests.o BaumWelch.o NewtonsMethod.o read_in_parameters.o read.o main_train.o grabTrainingExamples.o interval_tree.o \
		 split.o write.o main_segment.o viterbi.o validate.o ParamWrapper.o
	@${CXX} ${CXXFLAGS} ${PWD}/unittests.o  ${PWD}/BaumWelch.o \
		${PWD}/NewtonsMethod.o ${PWD}/read_in_parameters.o ${PWD}/read.o \
		${PWD}/interval_tree.o ${PWD}/grabTrainingExamples.o ${PWD}/main_train.o \
		${PWD}/split.o ${PWD}/write.o ${PWD}/main_segment.o ${PWD}/viterbi.o ${PWD}/validate.o ${PWD}/ParamWrapper.o -o unittests
	
	@echo Unit test executable built. Run ./unittests <training bed4 file> <bed4 dataset with positive and negative reads>
	

clean:
	@rm -f *.o

fullclean:
	@rm -f *.o FStitch unittest
