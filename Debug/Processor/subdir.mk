################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Processor/Binary_File_IO.cpp \
../Processor/Buffer.cpp \
../Processor/Data_Files.cpp \
../Processor/ExternalClients.cpp \
../Processor/Input.cpp \
../Processor/Instruction.cpp \
../Processor/Machine.cpp \
../Processor/Memory.cpp \
../Processor/Online-Thread.cpp \
../Processor/PrivateOutput.cpp \
../Processor/Processor.cpp \
../Processor/Program.cpp 

O_SRCS += \
../Processor/Binary_File_IO.o \
../Processor/Buffer.o \
../Processor/Data_Files.o \
../Processor/ExternalClients.o \
../Processor/Input.o \
../Processor/Instruction.o \
../Processor/Machine.o \
../Processor/Memory.o \
../Processor/Online-Thread.o \
../Processor/PrivateOutput.o \
../Processor/Processor.o \
../Processor/Program.o 

OBJS += \
./Processor/Binary_File_IO.o \
./Processor/Buffer.o \
./Processor/Data_Files.o \
./Processor/ExternalClients.o \
./Processor/Input.o \
./Processor/Instruction.o \
./Processor/Machine.o \
./Processor/Memory.o \
./Processor/Online-Thread.o \
./Processor/PrivateOutput.o \
./Processor/Processor.o \
./Processor/Program.o 

CPP_DEPS += \
./Processor/Binary_File_IO.d \
./Processor/Buffer.d \
./Processor/Data_Files.d \
./Processor/ExternalClients.d \
./Processor/Input.d \
./Processor/Instruction.d \
./Processor/Machine.d \
./Processor/Memory.d \
./Processor/Online-Thread.d \
./Processor/PrivateOutput.d \
./Processor/Processor.d \
./Processor/Program.d 


# Each subdirectory must supply rules for building sources it contributes
Processor/%.o: ../Processor/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


