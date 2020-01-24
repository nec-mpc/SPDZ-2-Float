################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Networking/Player.cpp \
../Networking/Receiver.cpp \
../Networking/STS.cpp \
../Networking/Sender.cpp \
../Networking/Server.cpp \
../Networking/ServerSocket.cpp \
../Networking/sockets.cpp 

O_SRCS += \
../Networking/Player.o \
../Networking/Receiver.o \
../Networking/STS.o \
../Networking/Sender.o \
../Networking/Server.o \
../Networking/ServerSocket.o \
../Networking/sockets.o 

OBJS += \
./Networking/Player.o \
./Networking/Receiver.o \
./Networking/STS.o \
./Networking/Sender.o \
./Networking/Server.o \
./Networking/ServerSocket.o \
./Networking/sockets.o 

CPP_DEPS += \
./Networking/Player.d \
./Networking/Receiver.d \
./Networking/STS.d \
./Networking/Sender.d \
./Networking/Server.d \
./Networking/ServerSocket.d \
./Networking/sockets.d 


# Each subdirectory must supply rules for building sources it contributes
Networking/%.o: ../Networking/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


