################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../nauty/naugraph.c \
../nauty/nautil.c \
../nauty/nauty.c 

OBJS += \
./nauty/naugraph.o \
./nauty/nautil.o \
./nauty/nauty.o 

C_DEPS += \
./nauty/naugraph.d \
./nauty/nautil.d \
./nauty/nauty.d 


# Each subdirectory must supply rules for building sources it contributes
nauty/%.o: ../nauty/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


