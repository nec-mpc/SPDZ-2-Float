################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Math/Integer.cpp \
../Math/Setup.cpp \
../Math/Share.cpp \
../Math/Subroutines.cpp \
../Math/Zp_Data.cpp \
../Math/bigint.cpp \
../Math/gf2n.cpp \
../Math/gf2nlong.cpp \
../Math/gfp.cpp \
../Math/modp.cpp 

O_SRCS += \
../Math/Integer.o \
../Math/Setup.o \
../Math/Share.o \
../Math/Subroutines.o \
../Math/Zp_Data.o \
../Math/bigint.o \
../Math/gf2n.o \
../Math/gf2nlong.o \
../Math/gfp.o \
../Math/modp.o 

OBJS += \
./Math/Integer.o \
./Math/Setup.o \
./Math/Share.o \
./Math/Subroutines.o \
./Math/Zp_Data.o \
./Math/bigint.o \
./Math/gf2n.o \
./Math/gf2nlong.o \
./Math/gfp.o \
./Math/modp.o 

CPP_DEPS += \
./Math/Integer.d \
./Math/Setup.d \
./Math/Share.d \
./Math/Subroutines.d \
./Math/Zp_Data.d \
./Math/bigint.d \
./Math/gf2n.d \
./Math/gf2nlong.d \
./Math/gfp.d \
./Math/modp.d 


# Each subdirectory must supply rules for building sources it contributes
Math/%.o: ../Math/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


