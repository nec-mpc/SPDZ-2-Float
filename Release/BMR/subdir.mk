################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BMR/AndJob.cpp \
../BMR/BooleanCircuit.cpp \
../BMR/CommonParty.cpp \
../BMR/GarbledGate.cpp \
../BMR/Key.cpp \
../BMR/Party.cpp \
../BMR/Register.cpp \
../BMR/SpdzWire.cpp \
../BMR/TrustedParty.cpp \
../BMR/aes.cpp \
../BMR/msg_types.cpp \
../BMR/prf.cpp \
../BMR/proto_utils.cpp 

OBJS += \
./BMR/AndJob.o \
./BMR/BooleanCircuit.o \
./BMR/CommonParty.o \
./BMR/GarbledGate.o \
./BMR/Key.o \
./BMR/Party.o \
./BMR/Register.o \
./BMR/SpdzWire.o \
./BMR/TrustedParty.o \
./BMR/aes.o \
./BMR/msg_types.o \
./BMR/prf.o \
./BMR/proto_utils.o 

CPP_DEPS += \
./BMR/AndJob.d \
./BMR/BooleanCircuit.d \
./BMR/CommonParty.d \
./BMR/GarbledGate.d \
./BMR/Key.d \
./BMR/Party.d \
./BMR/Register.d \
./BMR/SpdzWire.d \
./BMR/TrustedParty.d \
./BMR/aes.d \
./BMR/msg_types.d \
./BMR/prf.d \
./BMR/proto_utils.d 


# Each subdirectory must supply rules for building sources it contributes
BMR/%.o: ../BMR/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


