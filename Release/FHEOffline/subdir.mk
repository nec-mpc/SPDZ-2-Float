################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../FHEOffline/CutAndChooseMachine.cpp \
../FHEOffline/DataSetup.cpp \
../FHEOffline/DistDecrypt.cpp \
../FHEOffline/DistKeyGen.cpp \
../FHEOffline/EncCommit.cpp \
../FHEOffline/FHE-Subroutines.cpp \
../FHEOffline/FullSetup.cpp \
../FHEOffline/Multiplier.cpp \
../FHEOffline/PairwiseGenerator.cpp \
../FHEOffline/PairwiseMachine.cpp \
../FHEOffline/PairwiseSetup.cpp \
../FHEOffline/Producer.cpp \
../FHEOffline/Proof.cpp \
../FHEOffline/Prover.cpp \
../FHEOffline/Reshare.cpp \
../FHEOffline/Sacrificing.cpp \
../FHEOffline/SimpleDistDecrypt.cpp \
../FHEOffline/SimpleEncCommit.cpp \
../FHEOffline/SimpleGenerator.cpp \
../FHEOffline/SimpleMachine.cpp \
../FHEOffline/Verifier.cpp 

OBJS += \
./FHEOffline/CutAndChooseMachine.o \
./FHEOffline/DataSetup.o \
./FHEOffline/DistDecrypt.o \
./FHEOffline/DistKeyGen.o \
./FHEOffline/EncCommit.o \
./FHEOffline/FHE-Subroutines.o \
./FHEOffline/FullSetup.o \
./FHEOffline/Multiplier.o \
./FHEOffline/PairwiseGenerator.o \
./FHEOffline/PairwiseMachine.o \
./FHEOffline/PairwiseSetup.o \
./FHEOffline/Producer.o \
./FHEOffline/Proof.o \
./FHEOffline/Prover.o \
./FHEOffline/Reshare.o \
./FHEOffline/Sacrificing.o \
./FHEOffline/SimpleDistDecrypt.o \
./FHEOffline/SimpleEncCommit.o \
./FHEOffline/SimpleGenerator.o \
./FHEOffline/SimpleMachine.o \
./FHEOffline/Verifier.o 

CPP_DEPS += \
./FHEOffline/CutAndChooseMachine.d \
./FHEOffline/DataSetup.d \
./FHEOffline/DistDecrypt.d \
./FHEOffline/DistKeyGen.d \
./FHEOffline/EncCommit.d \
./FHEOffline/FHE-Subroutines.d \
./FHEOffline/FullSetup.d \
./FHEOffline/Multiplier.d \
./FHEOffline/PairwiseGenerator.d \
./FHEOffline/PairwiseMachine.d \
./FHEOffline/PairwiseSetup.d \
./FHEOffline/Producer.d \
./FHEOffline/Proof.d \
./FHEOffline/Prover.d \
./FHEOffline/Reshare.d \
./FHEOffline/Sacrificing.d \
./FHEOffline/SimpleDistDecrypt.d \
./FHEOffline/SimpleEncCommit.d \
./FHEOffline/SimpleGenerator.d \
./FHEOffline/SimpleMachine.d \
./FHEOffline/Verifier.d 


# Each subdirectory must supply rules for building sources it contributes
FHEOffline/%.o: ../FHEOffline/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


