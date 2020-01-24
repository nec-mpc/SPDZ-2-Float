################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Auth/MAC_Check.cpp \
../Auth/Subroutines.cpp \
../Auth/Summer.cpp \
../Auth/fake-stuff.cpp 

O_SRCS += \
../Auth/MAC_Check.o \
../Auth/Subroutines.o \
../Auth/Summer.o \
../Auth/fake-stuff.o 

OBJS += \
./Auth/MAC_Check.o \
./Auth/Subroutines.o \
./Auth/Summer.o \
./Auth/fake-stuff.o 

CPP_DEPS += \
./Auth/MAC_Check.d \
./Auth/Subroutines.d \
./Auth/Summer.d \
./Auth/fake-stuff.d 


# Each subdirectory must supply rules for building sources it contributes
Auth/%.o: ../Auth/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


