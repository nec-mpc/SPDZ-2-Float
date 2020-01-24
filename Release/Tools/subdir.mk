################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tools/Commit.cpp \
../Tools/Config.cpp \
../Tools/FlexBuffer.cpp \
../Tools/Lock.cpp \
../Tools/MMO.cpp \
../Tools/OfflineMachineBase.cpp \
../Tools/Signal.cpp \
../Tools/aes-ni.cpp \
../Tools/aes.cpp \
../Tools/mkpath.cpp \
../Tools/names.cpp \
../Tools/octetStream.cpp \
../Tools/random.cpp \
../Tools/sha1.cpp \
../Tools/time-func.cpp 

O_SRCS += \
../Tools/Commit.o \
../Tools/Config.o \
../Tools/FlexBuffer.o \
../Tools/Lock.o \
../Tools/MMO.o \
../Tools/OfflineMachineBase.o \
../Tools/Signal.o \
../Tools/aes-ni.o \
../Tools/aes.o \
../Tools/mkpath.o \
../Tools/names.o \
../Tools/octetStream.o \
../Tools/random.o \
../Tools/sha1.o \
../Tools/time-func.o 

OBJS += \
./Tools/Commit.o \
./Tools/Config.o \
./Tools/FlexBuffer.o \
./Tools/Lock.o \
./Tools/MMO.o \
./Tools/OfflineMachineBase.o \
./Tools/Signal.o \
./Tools/aes-ni.o \
./Tools/aes.o \
./Tools/mkpath.o \
./Tools/names.o \
./Tools/octetStream.o \
./Tools/random.o \
./Tools/sha1.o \
./Tools/time-func.o 

CPP_DEPS += \
./Tools/Commit.d \
./Tools/Config.d \
./Tools/FlexBuffer.d \
./Tools/Lock.d \
./Tools/MMO.d \
./Tools/OfflineMachineBase.d \
./Tools/Signal.d \
./Tools/aes-ni.d \
./Tools/aes.d \
./Tools/mkpath.d \
./Tools/names.d \
./Tools/octetStream.d \
./Tools/random.d \
./Tools/sha1.d \
./Tools/time-func.d 


# Each subdirectory must supply rules for building sources it contributes
Tools/%.o: ../Tools/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


