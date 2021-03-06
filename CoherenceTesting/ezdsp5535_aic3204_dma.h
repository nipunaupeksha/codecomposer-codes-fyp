//////////////////////////////////////////////////////////////////////////////
// *
// * FILENAME
// *  ezdsp5535_aic3204_dma.h
// *
// * DESCRIPTION
// *   Header file for ezdsp5535_aic3204_dma.c
// *
// * REVISION                                                             
// *   Revision: 1.00	                                                       
// *   Author  : R.M.N. Upeksha 
// *
//////////////////////////////////////////////////////////////////////////////
#ifndef MY_DMA_PING_PONG_REGISTER_SETUP_H_
#define MY_DMA_PING_PONG_REGISTER_SETUP_H_

#define AUDIO_IO_SIZE  128		// DMA transfer size
#define PING_PONG_SIZE (2 * AUDIO_IO_SIZE)

/* ============ Some Data Types Macros ======================= */
#define Uint32  unsigned long
#define Uint16  unsigned short
#define Uint8   unsigned char
#define Int32   int
#define Int16   short
#define Int8    char


#define GPIO20  0x14 //for 20 pins its 0x14 and for 26 pins its GPIO26 0x1A

// My_DMA
/* ==================== System Registers ====================== */

#define EBSR      		    *(volatile ioport Uint16*)(0x1c00)
#define PCGCR1         		*(volatile ioport Uint16*)(0x1c02)
#define PCGCR2         		*(volatile ioport Uint16*)(0x1c03)
#define PSRCR         		*(volatile ioport Uint16*)(0x1c04)
#define PRCR       		    *(volatile ioport Uint16*)(0x1c05)
#define ODSCR		       	*(volatile ioport Uint16*)(0x1c16)
#define PDINHIBR1        	*(volatile ioport Uint16*)(0x1c17)

/* =================== GPIO Registers ========================= */

#define IODIR1      		*(volatile ioport Uint16*)(0x1c06)
#define IODIR2          	*(volatile ioport Uint16*)(0x1c07)
#define IOINDATA1		   	*(volatile ioport Uint16*)(0x1c08)
#define IOINDATA2		   	*(volatile ioport Uint16*)(0x1c09)
#define IODATAOUT1		  	*(volatile ioport Uint16*)(0x1c0a)
#define IODATAOUT2  		*(volatile ioport Uint16*)(0x1c0b)

/* ================== I2C Registers =========================== */

#define IER    	       		*(volatile ioport Uint16*)(0x1A04) //IMR
#define STR    	       		*(volatile ioport Uint16*)(0x1A08)
#define CLKL           		*(volatile ioport Uint16*)(0x1A0C)
#define CLKH           		*(volatile ioport Uint16*)(0x1A10)
#define CNT    		   		*(volatile ioport Uint16*)(0x1A14)
#define DRR    		   		*(volatile ioport Uint16*)(0x1A18)
#define SAR    	       		*(volatile ioport Uint16*)(0x1A1C)
#define DXR    	       		*(volatile ioport Uint16*)(0x1A20)
#define MDR            		*(volatile ioport Uint16*)(0x1A24)
#define EDR    	       		*(volatile ioport Uint16*)(0x1A2C)
#define PSC    	       		*(volatile ioport Uint16*)(0x1A30)

/* ================== I2S Registers =========================== */

#define CR            		*(volatile ioport Uint16*)(0x2800)
#define SRGR          		*(volatile ioport Uint16*)(0x2804)
#define W0_LSW_W      		*(volatile ioport Uint16*)(0x2808)
#define W0_MSW_W      		*(volatile ioport Uint16*)(0x2809)
#define W1_LSW_W      		*(volatile ioport Uint16*)(0x280C)
#define W1_MSW_W      		*(volatile ioport Uint16*)(0x280D)
#define IR            		*(volatile ioport Uint16*)(0x2810)
#define ICMR          		*(volatile ioport Uint16*)(0x2814)
#define W0_LSW_R      		*(volatile ioport Uint16*)(0x2828)
#define W0_MSW_R      		*(volatile ioport Uint16*)(0x2829)
#define W1_LSW_R      		*(volatile ioport Uint16*)(0x282C)
#define W1_MSW_R      		*(volatile ioport Uint16*)(0x282D)

/*================= FFT ======================================= */
#define FFTPOINTS 512

/* ==================== Some CPU Registers ==================== */
#define IER0        			 *(volatile unsigned *)0x0000
#define IFR0        			 *(volatile unsigned *)0x0001
#define IER1        			 *(volatile unsigned *)0x0045
#define IFR1        			 *(volatile unsigned *)0x0046
//#define PRCR       				 *(volatile ioport Uint16*)(0x1C05)

/* ==================== General DMA Registers ==================== */
#define DMA_IFR     			 *(ioport volatile unsigned *)0x1C30    // DMA Interrupt Flag Register
#define DMA_IER     			 *(ioport volatile unsigned *)0x1C31    // DMA Interrupt Enable Register

/* ==================== DMA Controller 1 ==================== */
#define DMA1_CESR1 *(ioport volatile unsigned *)0x1C1C //DMA1 Channel Event Source Register 1
#define DMA1_CESR2 *(ioport volatile unsigned *)0x1C1D //DMA1 Channel Event Source Register 2
/* ------------------------ Channel 0 ----------------------- */
#define DMA1_CH0_SSAL *(ioport volatile unsigned *)0x0D00 //Channel 0 Source Start Address (Lower Part) Register
#define DMA1_CH0_SSAU *(ioport volatile unsigned *)0x0D01 //Channel 0 Source Start Address (Upper Part) Register
#define DMA1_CH0_DSAL *(ioport volatile unsigned *)0x0D02 //Channel 0 Destination Start Address (Lower Part) Register
#define DMA1_CH0_DSAU *(ioport volatile unsigned *)0x0D03 //Channel 0 Destination Start Address (Upper Part) Register
#define DMA1_CH0_TCR1 *(ioport volatile unsigned *)0x0D04 //Channel 0 Transfer Control Register 1
#define DMA1_CH0_TCR2 *(ioport volatile unsigned *)0x0D05 //Channel 0 Transfer Control Register 2
/* ---------------------------------------------------------- */

/* ------------------------ Channel 1 ----------------------- */
#define DMA1_CH1_SSAL *(ioport volatile unsigned *)0x0D20 //Channel 1 Source Start Address (Lower Part) Register
#define DMA1_CH1_SSAU *(ioport volatile unsigned *)0x0D21 //Channel 1 Source Start Address (Upper Part) Register
#define DMA1_CH1_DSAL *(ioport volatile unsigned *)0x0D22 //Channel 1 Destination Start Address (Lower Part) Register
#define DMA1_CH1_DSAU *(ioport volatile unsigned *)0x0D23 //Channel 1 Destination Start Address (Upper Part) Register
#define DMA1_CH1_TCR1 *(ioport volatile unsigned *)0x0D24 //Channel 1 Transfer Control Register 1
#define DMA1_CH1_TCR2 *(ioport volatile unsigned *)0x0D25 //Channel 1 Transfer Control Register 2
/* ---------------------------------------------------------- */

/* ------------------------ Channel 2 ----------------------- */
#define DMA1_CH2_SSAL *(ioport volatile unsigned *)0x0D40 //Channel 2 Source Start Address (Lower Part) Register
#define DMA1_CH2_SSAU *(ioport volatile unsigned *)0x0D41 //Channel 2 Source Start Address (Upper Part) Register
#define DMA1_CH2_DSAL *(ioport volatile unsigned *)0x0D42 //Channel 2 Destination Start Address (Lower Part) Register
#define DMA1_CH2_DSAU *(ioport volatile unsigned *)0x0D43 //Channel 2 Destination Start Address (Upper Part) Register
#define DMA1_CH2_TCR1 *(ioport volatile unsigned *)0x0D44 //Channel 2 Transfer Control Register 1
#define DMA1_CH2_TCR2 *(ioport volatile unsigned *)0x0D45 //Channel 2 Transfer Control Register 2
/* ---------------------------------------------------------- */

/* ------------------------ Channel 3 ----------------------- */
#define DMA1_CH3_SSAL *(ioport volatile unsigned *)0x0D60 //Channel 3 Source Start Address (Lower Part) Register
#define DMA1_CH3_SSAU *(ioport volatile unsigned *)0x0D61 //Channel 3 Source Start Address (Upper Part) Register
#define DMA1_CH3_DSAL *(ioport volatile unsigned *)0x0D62 //Channel 3 Destination Start Address (Lower Part) Register
#define DMA1_CH3_DSAU *(ioport volatile unsigned *)0x0D63 //Channel 3 Destination Start Address (Upper Part) Register
#define DMA1_CH3_TCR1 *(ioport volatile unsigned *)0x0D64 //Channel 3 Transfer Control Register 1
#define DMA1_CH3_TCR2 *(ioport volatile unsigned *)0x0D65 //Channel 3 Transfer Control Register 2
/* ---------------------------------------------------------- */

/* =========================== FUNCTION PROTOTYPES ============================= */

typedef struct{
	float real; //real part of the complex number
	float imag;// imaginary part of the complex number
}complex;

void aic3204_dma_init(void);
void aic3204_read_block(Int16* buffer_left, Int16* buffer_right);
void aic3204_write_block(Int16* buffer_left, Int16* buffer_right);


void GPIO_direction( Uint16 number, Uint16 direction );
void GPIO_dataout( Uint16 number, Uint16 output );

void wait(Uint32 delay);
void waitusec(Uint32 usec);


#endif /*MY_DMA_PING_PONG_REGISTER_SETUP_H_*/
