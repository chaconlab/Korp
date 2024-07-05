################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../korpe.cpp 

OBJS += \
./korpe.o 

CPP_DEPS += \
./korpe.d 


# Each subdirectory must supply rules for building sources it contributes
korpe.o: ../korpe.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../../src -O3 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"korpe.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


