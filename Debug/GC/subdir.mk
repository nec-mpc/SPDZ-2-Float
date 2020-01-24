################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../GC/FakeSecret.cpp \
../GC/Instruction.cpp \
../GC/Machine.cpp \
../GC/Memory.cpp \
../GC/Processor.cpp \
../GC/Program.cpp \
../GC/Secret.cpp 

OBJS += \
./GC/FakeSecret.o \
./GC/Instruction.o \
./GC/Machine.o \
./GC/Memory.o \
./GC/Processor.o \
./GC/Program.o \
./GC/Secret.o 

CPP_DEPS += \
./GC/FakeSecret.d \
./GC/Instruction.d \
./GC/Machine.d \
./GC/Memory.d \
./GC/Processor.d \
./GC/Program.d \
./GC/Secret.d 


# Each subdirectory must supply rules for building sources it contributes
GC/%.o: ../GC/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


