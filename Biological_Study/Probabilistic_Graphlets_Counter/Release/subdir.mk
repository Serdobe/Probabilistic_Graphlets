################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../CmdLine.cpp \
../Conditions.cpp \
../Error.cpp \
../Esu.cpp \
../GTrie.cpp \
../GraphMatrix.cpp \
../GraphTree.cpp \
../GraphUtils.cpp \
../Isomorphism.cpp \
../Random.cpp \
../Timer.cpp \
../main.cpp 

OBJS += \
./CmdLine.o \
./Conditions.o \
./Error.o \
./Esu.o \
./GTrie.o \
./GraphMatrix.o \
./GraphTree.o \
./GraphUtils.o \
./Isomorphism.o \
./Random.o \
./Timer.o \
./main.o 

CPP_DEPS += \
./CmdLine.d \
./Conditions.d \
./Error.d \
./Esu.d \
./GTrie.d \
./GraphMatrix.d \
./GraphTree.d \
./GraphUtils.d \
./Isomorphism.d \
./Random.d \
./Timer.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


