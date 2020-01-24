################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ExternalIO/bankers-bonus-client.cpp \
../ExternalIO/bankers-bonus-commsec-client.cpp 

OBJS += \
./ExternalIO/bankers-bonus-client.o \
./ExternalIO/bankers-bonus-commsec-client.o 

CPP_DEPS += \
./ExternalIO/bankers-bonus-client.d \
./ExternalIO/bankers-bonus-commsec-client.d 


# Each subdirectory must supply rules for building sources it contributes
ExternalIO/%.o: ../ExternalIO/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


