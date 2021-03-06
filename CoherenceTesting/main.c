/*****************************************************************************/
/*                                                                           */
/* FILENAME                                                                  */
/* 	 main.c                                                               */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   main function for codec on the EZDSP 5535 USB Stick.          */
/*                                                                           */
/* REVISION                                                                  */
/*   Revision: 1.00	                                                         */
/*   Author  : R.M.N. Upeksha                                                  */
/*****************************************************************************/

#include "stdio.h"
#include "ezdsp5535.h"
#include "ezdsp5535_i2c.h"
#include "aic3204.h"
#include "ezdsp5535_aic3204_dma.h"
#include "ezdsp5535_i2s.h"
#include "ezdsp5535_sar.h"
#include "math.h"

#define SAMPLE_RATE 8000L

#define PI 3.14159265

#pragma DATA_ALIGN(sampleBufferL, 4)
Int16 sampleBufferL[AUDIO_IO_SIZE];
#pragma DATA_ALIGN(sampleBufferR, 4)
Int16 sampleBufferR[AUDIO_IO_SIZE];
short tempBuff[128];

void main(void)
{
	int i = 0;

	EZDSP5535_init();

	EZDSP5535_SAR_init();

	printf("\n Getting audio signals \n");

	aic3204_hardware_init();

	aic3204_init();

	aic3204_dma_init();

	set_sampling_frequency_and_gain(SAMPLE_RATE, 0);

	while (1)
	{
		aic3204_read_block(sampleBufferL, sampleBufferR);
		printf("%d", ++i);
		/* Your code here */
		if (i == 100)
		{
			break;
		}

		aic3204_write_block(sampleBufferR, sampleBufferR);
	}

	/* Prekid veze sa AIC3204 kodekom */
	aic3204_disable();

	printf("\n***End Program***\n");
	SW_BREAKPOINT;
}
