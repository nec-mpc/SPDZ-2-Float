################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Scripts/gen_input_f2n.cpp \
../Scripts/gen_input_fp.cpp 

OBJS += \
./Scripts/gen_input_f2n.o \
./Scripts/gen_input_fp.o 

CPP_DEPS += \
./Scripts/gen_input_f2n.d \
./Scripts/gen_input_fp.d 


# Each subdirectory must supply rules for building sources it contributes
Scripts/%.o: ../Scripts/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


