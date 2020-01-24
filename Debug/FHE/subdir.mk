################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../FHE/Ciphertext.cpp \
../FHE/DiscreteGauss.cpp \
../FHE/FFT.cpp \
../FHE/FFT_Data.cpp \
../FHE/FHE_Keys.cpp \
../FHE/FHE_Params.cpp \
../FHE/Matrix.cpp \
../FHE/NTL-Subs.cpp \
../FHE/NoiseBounds.cpp \
../FHE/P2Data.cpp \
../FHE/PPData.cpp \
../FHE/Plaintext.cpp \
../FHE/QGroup.cpp \
../FHE/Random_Coins.cpp \
../FHE/Ring.cpp \
../FHE/Ring_Element.cpp \
../FHE/Rq_Element.cpp 

OBJS += \
./FHE/Ciphertext.o \
./FHE/DiscreteGauss.o \
./FHE/FFT.o \
./FHE/FFT_Data.o \
./FHE/FHE_Keys.o \
./FHE/FHE_Params.o \
./FHE/Matrix.o \
./FHE/NTL-Subs.o \
./FHE/NoiseBounds.o \
./FHE/P2Data.o \
./FHE/PPData.o \
./FHE/Plaintext.o \
./FHE/QGroup.o \
./FHE/Random_Coins.o \
./FHE/Ring.o \
./FHE/Ring_Element.o \
./FHE/Rq_Element.o 

CPP_DEPS += \
./FHE/Ciphertext.d \
./FHE/DiscreteGauss.d \
./FHE/FFT.d \
./FHE/FFT_Data.d \
./FHE/FHE_Keys.d \
./FHE/FHE_Params.d \
./FHE/Matrix.d \
./FHE/NTL-Subs.d \
./FHE/NoiseBounds.d \
./FHE/P2Data.d \
./FHE/PPData.d \
./FHE/Plaintext.d \
./FHE/QGroup.d \
./FHE/Random_Coins.d \
./FHE/Ring.d \
./FHE/Ring_Element.d \
./FHE/Rq_Element.d 


# Each subdirectory must supply rules for building sources it contributes
FHE/%.o: ../FHE/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


