################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BMR/network/Client.cpp \
../BMR/network/Node.cpp \
../BMR/network/Server.cpp \
../BMR/network/utils.cpp 

OBJS += \
./BMR/network/Client.o \
./BMR/network/Node.o \
./BMR/network/Server.o \
./BMR/network/utils.o 

CPP_DEPS += \
./BMR/network/Client.d \
./BMR/network/Node.d \
./BMR/network/Server.d \
./BMR/network/utils.d 


# Each subdirectory must supply rules for building sources it contributes
BMR/network/%.o: ../BMR/network/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


