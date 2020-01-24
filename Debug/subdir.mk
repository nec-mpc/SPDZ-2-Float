################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Check-Offline.cpp \
../Fake-Offline.cpp \
../Player-Online.cpp \
../Server.cpp \
../bmr-program-party.cpp \
../bmr-program-tparty.cpp \
../check-passive.cpp \
../client-setup.cpp \
../cnc-offline.cpp \
../ot-offline.cpp \
../pairwise-offline.cpp \
../simple-offline.cpp \
../spdz2-offline.cpp 

OBJS += \
./Check-Offline.o \
./Fake-Offline.o \
./Player-Online.o \
./Server.o \
./bmr-program-party.o \
./bmr-program-tparty.o \
./check-passive.o \
./client-setup.o \
./cnc-offline.o \
./ot-offline.o \
./pairwise-offline.o \
./simple-offline.o \
./spdz2-offline.o 

CPP_DEPS += \
./Check-Offline.d \
./Fake-Offline.d \
./Player-Online.d \
./Server.d \
./bmr-program-party.d \
./bmr-program-tparty.d \
./check-passive.d \
./client-setup.d \
./cnc-offline.d \
./ot-offline.d \
./pairwise-offline.d \
./simple-offline.d \
./spdz2-offline.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


