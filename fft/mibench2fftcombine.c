#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

void fft_double (
    unsigned  NumSamples,          /* must be a power of 2 */
    int       InverseTransform,    /* 0=forward FFT, 1=inverse FFT */
    double   *RealIn,              /* array of input's real samples */
    double   *ImaginaryIn,         /* array of input's imag samples */
    double   *RealOut,             /* array of output's reals */
    double   *ImaginaryOut );      /* array of output's imaginaries */


void fft_float (
    unsigned  NumSamples,          /* must be a power of 2 */
    int       InverseTransform,    /* 0=forward FFT, 1=inverse FFT */
    float    *RealIn,              /* array of input's real samples */
    float    *ImaginaryIn,         /* array of input's imag samples */
    float    *RealOut,             /* array of output's reals */
    float    *ImaginaryOut );      /* array of output's imaginaries */


int IsPowerOfTwo ( unsigned x );
unsigned NumberOfBitsNeeded ( unsigned PowerOfTwo );
unsigned ReverseBits ( unsigned index, unsigned NumBits );

#define fprintf(...)
    
double Index_to_frequency ( unsigned NumSamples, unsigned Index );

#ifdef __cplusplus
}
#endif


/*--- end of file fourier.h ---*/

#ifndef __DDC_DDC_H
#define __DDC_DDC_H

// If you add something to DDCRET, please add the appropriate string
// to the function DDCRET_String() in the file 'source\ddcret.cpp'.

enum DDCRET
{
   DDC_SUCCESS,           // The operation succeded
   DDC_FAILURE,           // The operation failed for unspecified reasons
   DDC_OUT_OF_MEMORY,     // Operation failed due to running out of memory
   DDC_FILE_ERROR,        // Operation encountered file I/O error
   DDC_INVALID_CALL,      // Operation was called with invalid parameters
   DDC_USER_ABORT,        // Operation was aborted by the user
   DDC_INVALID_FILE       // File format does not match
};


char *DDCRET_String ( DDCRET );   // See source\ddcret.cpp


#define  TRUE     1
#define  FALSE    0

typedef int dBOOLEAN;

typedef unsigned char BYTE;

typedef unsigned char        UINT8;
typedef signed   char        INT8;

typedef unsigned short int   UINT16;
typedef signed   short int   INT16;
typedef unsigned long  int   UINT32;
typedef signed   long  int   INT32;

#ifdef __BORLANDC__
	#if sizeof(UINT16) != 2
	  #error Need to fix UINT16 and INT16
	#endif

	#if sizeof(UINT32) != 4
	  #error Need to fix UINT32 and INT32
	#endif
#endif

#endif /* __DDC_DDC_H */

/*--- end of file ddc.h ---*/

#ifndef __ddcmath_h
#define __ddcmath_h

#define  DDC_PI  (3.14159265358979323846)

#endif /* __ddcmath_h */

/*--- end of file ddcmath.h ---*/

#define TRUE  1
#define FALSE 0

#define BITS_PER_WORD   (sizeof(unsigned) * 8)


int IsPowerOfTwo ( unsigned x )
{
    if ( x < 2 )
        return FALSE;

    if ( x & (x-1) )        // Thanks to 'byang' for this cute trick!
        return FALSE;

    return TRUE;
}


unsigned NumberOfBitsNeeded ( unsigned PowerOfTwo )
{
    unsigned i;

    if ( PowerOfTwo < 2 )
    {
        fprintf (
            stderr,
            ">>> Error in fftmisc.c: argument %d to NumberOfBitsNeeded is too small.\n",
            PowerOfTwo );

        exit(1);
    }

    for ( i=0; ; i++ )
    {
        if ( PowerOfTwo & (1 << i) )
            return i;
    }
}



unsigned ReverseBits ( unsigned index, unsigned NumBits )
{
    unsigned i, rev;

    for ( i=rev=0; i < NumBits; i++ )
    {
        rev = (rev << 1) | (index & 1);
        index >>= 1;
    }

    return rev;
}


double Index_to_frequency ( unsigned NumSamples, unsigned Index )
{
    if ( Index >= NumSamples )
        return 0.0;
    else if ( Index <= NumSamples/2 )
        return (double)Index / (double)NumSamples;

    return -(double)(NumSamples-Index) / (double)NumSamples;
}


/*--- end of file fftmisc.c---*/

#define CHECKPOINTER(p)  CheckPointer(p,#p)

static void CheckPointer ( void *p, char *name )
{
    if ( p == NULL )
    {
        fprintf ( stderr, "Error in fft_float():  %s == NULL\n", name );
        exit(1);
    }
}


void fft_float (
    unsigned  NumSamples,
    int       InverseTransform,
    float    *RealIn,
    float    *ImagIn,
    float    *RealOut,
    float    *ImagOut )
{
    unsigned NumBits;    /* Number of bits needed to store indices */
    unsigned i, j, k, n;
    unsigned BlockSize, BlockEnd;

    double angle_numerator = 2.0 * DDC_PI;
    double tr, ti;     /* temp real, temp imaginary */

    if ( !IsPowerOfTwo(NumSamples) )
    {
        fprintf (
            stderr,
            "Error in fft():  NumSamples=%u is not power of two\n",
            NumSamples );

        exit(1);
    }

    if ( InverseTransform )
        angle_numerator = -angle_numerator;

    CHECKPOINTER ( RealIn );
    CHECKPOINTER ( RealOut );
    CHECKPOINTER ( ImagOut );

    NumBits = NumberOfBitsNeeded ( NumSamples );

    /*
    **   Do simultaneous data copy and bit-reversal ordering into outputs...
    */

    for ( i=0; i < NumSamples; i++ )
    {
        j = ReverseBits ( i, NumBits );
        RealOut[j] = RealIn[i];
        ImagOut[j] = (ImagIn == NULL) ? 0.0 : ImagIn[i];
    }

    /*
    **   Do the FFT itself...
    */

    BlockEnd = 1;
    for ( BlockSize = 2; BlockSize <= NumSamples; BlockSize <<= 1 )
    {
        double delta_angle = angle_numerator / (double)BlockSize;
        double sm2 = sin ( -2 * delta_angle );
        double sm1 = sin ( -delta_angle );
        double cm2 = cos ( -2 * delta_angle );
        double cm1 = cos ( -delta_angle );
        double w = 2 * cm1;
        double ar[3], ai[3];

        for ( i=0; i < NumSamples; i += BlockSize )
        {
            ar[2] = cm2;
            ar[1] = cm1;

            ai[2] = sm2;
            ai[1] = sm1;

            for ( j=i, n=0; n < BlockEnd; j++, n++ )
            {
                ar[0] = w*ar[1] - ar[2];
                ar[2] = ar[1];
                ar[1] = ar[0];

                ai[0] = w*ai[1] - ai[2];
                ai[2] = ai[1];
                ai[1] = ai[0];

                k = j + BlockEnd;
                tr = ar[0]*RealOut[k] - ai[0]*ImagOut[k];
                ti = ar[0]*ImagOut[k] + ai[0]*RealOut[k];

                RealOut[k] = RealOut[j] - tr;
                ImagOut[k] = ImagOut[j] - ti;

                RealOut[j] += tr;
                ImagOut[j] += ti;
            }
        }

        BlockEnd = BlockSize;
    }

    /*
    **   Need to normalize if inverse transform...
    */

    if ( InverseTransform )
    {
        double denom = (double)NumSamples;

        for ( i=0; i < NumSamples; i++ )
        {
            RealOut[i] /= denom;
            ImagOut[i] /= denom;
        }
    }
}



#include "../bareBench.h"

int invfft=0;
unsigned MAXSIZE; // small 4096, 8192 inverse, 512 for memory-limited systems
unsigned MAXWAVES=4; //large has 8
static float realin[1024];
static float imagin[1024];
static float realout[1024];
static float imagout[1024];
static float Coeff[16];
static float Amp[16];

int old_main();
    
// main for benchmark purposes that does fft and inverse fft
int main() {
    MAXSIZE = 128;
    old_main();
    invfft = 1;
    MAXSIZE = 256;
    old_main();
    return 0;
}

int old_main() {
	unsigned i,j;
	float *RealIn;
	float *ImagIn;
	float *RealOut;
	float *ImagOut;
	float *coeff;
	float *amp;

		
 srand(1);

/* RealIn=(float*)malloc(sizeof(float)*MAXSIZE);
 ImagIn=(float*)malloc(sizeof(float)*MAXSIZE);
 RealOut=(float*)malloc(sizeof(float)*MAXSIZE);
 ImagOut=(float*)malloc(sizeof(float)*MAXSIZE);
 coeff=(float*)malloc(sizeof(float)*MAXWAVES);
 amp=(float*)malloc(sizeof(float)*MAXWAVES);
*/

 //Statically allocate
	RealIn = realin;
	ImagIn = imagin;
	RealOut = realout;
	ImagOut = imagout;
	coeff = Coeff;
	amp = Amp;

 /* Makes MAXWAVES waves of random amplitude and period */
	for(i=0;i<MAXWAVES;i++) 
	{
		coeff[i] = rand()%1000;
		amp[i] = rand()%1000;
	}
 for(i=0;i<MAXSIZE;i++) 
 {
   /*   RealIn[i]=rand();*/
	 RealIn[i]=0;
	 for(j=0;j<MAXWAVES;j++) 
	 {
		 /* randomly select sin or cos */
		 if (rand()%2)
		 {
		 		RealIn[i]+=coeff[j]*cos(amp[j]*i);
			}
		 else
		 {
		 	RealIn[i]+=coeff[j]*sin(amp[j]*i);
		 }
  	 ImagIn[i]=0;
	 }
 }

 /* regular*/
 fft_float (MAXSIZE,invfft,RealIn,ImagIn,RealOut,ImagOut);
 
 printf("RealOut:\n");
 for (i=0;i<MAXSIZE;i++)
   printf("%f \t", RealOut[i]);
 printf("\n");

printf("ImagOut:\n");
 for (i=0;i<MAXSIZE;i++)
   printf("%f \t", ImagOut[i]);
   printf("\n");

/* free(RealIn);
 free(ImagIn);
 free(RealOut);
 free(ImagOut);
 free(coeff);
 free(amp);
*/ return 0;


}
