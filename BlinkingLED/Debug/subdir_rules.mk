################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Each subdirectory must supply rules for building sources it contributes
blink_led.obj: ../blink_led.c $(GEN_OPTS) $(GEN_SRCS)
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	"E:/Software/Code Composer Studio/ccsv4/tools/compiler/c5500/bin/cl55" -v5515 -g --include_path="E:/Software/Code Composer Studio/ccsv4/tools/compiler/c5500/include" --include_path="E:/Software/Code Composer Studio/ccsv4/emulation/boards/ezdsp5535_v1/c55xx_csl/inc" --include_path="E:/Software/Code Composer Studio/ccsv4/emulation/boards/ezdsp5535_v1/include" --diag_warning=225 --ptrdiff_size=16 --memory_model=large --preproc_with_compile --preproc_dependency="blink_led.pp" $(GEN_OPTS_QUOTED) $(subst #,$(wildcard $(subst $(SPACE),\$(SPACE),$<)),"#")
	@echo 'Finished building: $<'
	@echo ' '

main.obj: ../main.c $(GEN_OPTS) $(GEN_SRCS)
	@echo 'Building file: $<'
	@echo 'Invoking: Compiler'
	"E:/Software/Code Composer Studio/ccsv4/tools/compiler/c5500/bin/cl55" -v5515 -g --include_path="E:/Software/Code Composer Studio/ccsv4/tools/compiler/c5500/include" --include_path="E:/Software/Code Composer Studio/ccsv4/emulation/boards/ezdsp5535_v1/c55xx_csl/inc" --include_path="E:/Software/Code Composer Studio/ccsv4/emulation/boards/ezdsp5535_v1/include" --diag_warning=225 --ptrdiff_size=16 --memory_model=large --preproc_with_compile --preproc_dependency="main.pp" $(GEN_OPTS_QUOTED) $(subst #,$(wildcard $(subst $(SPACE),\$(SPACE),$<)),"#")
	@echo 'Finished building: $<'
	@echo ' '


