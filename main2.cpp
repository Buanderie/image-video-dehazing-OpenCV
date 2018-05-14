/*
 * Copyright 2011 Juan Gabriel Gomila Salas <joan-by@hotmail.com>
 *
 * All rights reserved.
 *
 *
 * Patent warning:
 *
 * This file implement one algorithm possibly linked to the patent
 *
 * # N. Moroney, R.G. Beausoleil and I. Sobel, Local Color Correction, 
 * US Patent 6,822,762, (November 23, 2004)
 *
 * This file is made available for the exclusive aim of serving as
 * scientific tool to verify the soundness and completeness of the
 * algorithm description. Compilation, execution and redistribution
 * of this file may violate patents rights in certain countries.
 * The situation being different for every country and changing
 * over time, it is your responsibility to determine which patent
 * rights restrictions apply to you before you compile, use,
 * modify, or redistribute this file. A patent lawyer is qualified
 * to make this determination.
 * If and only if they don't conflict with any patent terms, you
 * can benefit from the following license terms attached to this
 * file.
 *
 * License:
 *
 * This program is provided for scientific and educational only:
 * you can use and/or modify it for these purposes, but you are
 * not allowed to redistribute this work or derivative works in
 * source or executable form. A license must be obtained from the
 * patent right holders for any other use.
 *
 *
 */

/**
 * @file localcolorcorrection.c
 * @brief Local Color Correction algorithm
 *
 * Parameters:
 *
 * input image
 * output image
 * radius (default: 40), radius of Gaussian kernel
 * option (default: 1), processing option
 *                      (1, separate processing of R, G and B channels;
 *                       2, processing of intensity channel;
 *                       3, processing of luma channel in YPbPr color space;
 *                       4, processing of luma channel in HSL color space )
 *
 * Read/write operations (png format) make use
 * of io_png.c and io_png.h, by Nicolas Limare
 *
 * Conversion operations between color spaces, make use
 * colorspace library, by Pascal Getreuer
 *
 * @author Juan Gabriel Gomila Salas <joan-by@hotmail.com>
 *
 */

/**
 * @mainpage Local Color Correction
 *
 * README.txt:
 * @verbinclude README.txt
 */


/**
 * @file   localcolorcorrection.c
 * @brief  Main executable file
 *
 *
 * @author Juan Gabriel Gomila Salas <joan-by@hotmail.com>
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io_png/io_png.h"
#include "colorspace/colorspace.h"


/*Scale and round a float in [0,1] to an integer in 0,1,...,255.*/
unsigned char ScaleAndRound(float X) 
{
    if (X<=0) {
        return 0;
    } else if (X>=1) {
        return 255;
    } else
        return (unsigned char)(255.0f*X+0.5f);
}



/* 
Converts 'input' array of color data from 
format RRRR...GGGG...BBBB to 3 independent arrays
'RR', 'GG', 'BB' each one of them storing the data
corresponding to Red, Green and Blue channels, respectively.
Data are converted to 'float' type and rescaled in the range [0, 1].
'size' is the length of the 'input' array. 
*/
void input2RGB(unsigned char *input,
               float **RR, float **GG, float **BB,
               int size) 
{
    float *R, *G, *B;
    int n;

    R= (float*) malloc(size*sizeof(float));
    G= (float*) malloc(size*sizeof(float));
    B= (float*) malloc(size*sizeof(float));

    for (n=0; n < size; n++) {
        R[n]=(float)input[n]/255.0f;
        G[n]=(float)input[size+n]/255.0f;
        B[n]=(float)input[2*size+n]/255.0f;
    }

    *RR= R;
    *GG= G;
    *BB= B;
}


/* 
Converts data stored in 3 independent arrays 'R', 'G' and 'B'
to a single array 'output', with format RRR...GGG...BBB.
Data are converted to unsigned char (range [0, 255]).
'size' is the length of the 'output' array. 
*/
void RGB2output(float *R, float *G, float *B,
                unsigned char *output, int size) 
{
    int n;

    for (n=0; n < size; n++) {
        output[n]=ScaleAndRound(R[n]);
        output[size+n]=ScaleAndRound(G[n]);
        output[2*size+n]=ScaleAndRound(B[n]);
    }

}



/*GLOBAL VARIABLES*/

unsigned char *Im, *Imi;
long Dxy;
int Dx, Dy;



/**
 * @brief main function call
 */
int main(int argc, const char **argv) 
{
    const char *namein, *nameout;
    int radius, option, n, nx, ny, t, nc;
    size_t Dxx, Dyy;
    unsigned char *img, *output;
    float *R, *G, *B, *Rout, *Gout, *Bout;
    float *I, *Iout, *Mask, *Mask2, *Pb, *Pr, *H, *S, *Kernel;
    float sum, r;
    num ii, ppb, ppr, hh, ss, ll, rr, gg, bb;


    if (argc < 3) {
        printf("Usage: localcolorcorrection input.png output.png \
[radius=40] [option=1]\n");
        printf(" option=1: process each color channel separately\n");
        printf(" option=2: process intensity channel and add original color \
information\n");
        printf(" option=3: process the Luma channel in the YPbPr color \
space\n");
        printf(" option=4: process the Luminance channel in the HSL color \
space\n");
        return EXIT_FAILURE;
    }

    radius=40;
    option=1;

    /*option=1: process each color channel separately
     option=2: process intensity channel and add original color information
     option=3: process the Luma channel in the YPbPr color space
     option=4: process the Luminance channel in the HSL color space*/

    namein = argv[1];
    nameout = argv[2];
    if (argc > 3) sscanf(argv[3], "%i", &radius);
    if (argc > 4) sscanf(argv[4], "%i", &option);

    if (option!=1&&option!=2&&option!=3&&option!=4){
        return EXIT_FAILURE;
    } 

    /* read the PNG input image */
    img = read_png_u8_rgb(namein, &Dxx, &Dyy);
    Dx=(int) Dxx;
    Dy=(int) Dyy;

    input2RGB(img, &R, &G, &B, Dx*Dy);

    I= (float*) malloc(Dx*Dy*sizeof(float));
    Mask= (float*) malloc(Dx*Dy*sizeof(float));
    Mask2= (float*) malloc(Dx*Dy*sizeof(float));

    Pb=NULL;
    Pr=NULL;
    H=NULL;
    S=NULL;
    Iout=NULL;

    switch (option) {
    case 3:
        Pb= (float*) malloc(Dx*Dy*sizeof(float));
        Pr= (float*) malloc(Dx*Dy*sizeof(float));
        /*I is in this case the Y component*/
        for (n=0; n < Dx*Dy; n++) {
            Rgb2Ypbpr(&ii, &ppb, &ppr, R[n], G[n], B[n]);
            I[n]=ii;
            Pb[n]=ppb;
            Pr[n]=ppr;
        }
        break;
    case 4:
        H= (float*) malloc(Dx*Dy*sizeof(float));
        S= (float*) malloc(Dx*Dy*sizeof(float));
        /*I is in this case the L component*/
        for (n=0; n < Dx*Dy; n++) {
            Rgb2Hsl(&hh, &ss, &ll, R[n], G[n], B[n]);
            H[n]=hh;
            S[n]=ss;
            I[n]=ll;
        }
        break;

    default:
        for (n=0; n < Dx*Dy; n++) {
            I[n]= (R[n]+ G[n]+B[n])/3.0f;
        }
        break;
    }


    for (n=0; n < Dx*Dy; n++) {
        Mask[n] = I[n];
        Mask2[n] =Mask[n];
    }

    /*Create convolution kernel*/
    Kernel =(float*) malloc((2*radius+1)*sizeof(float));

    if (radius > 0) {
        for (n=0; n< 2*radius+1; n++)
            Kernel[n] =
                (float) exp((double)(-2.0f*(-radius+n)*(-radius+n)/(radius*radius)));

        sum=0.0f;
        for (n=0; n< 2*radius+1; n++)
            sum = sum+Kernel[n];
        for (n=0; n< 2*radius+1; n++)
            Kernel[n]=(float) (Kernel[n]/sum);

        /*If the radius of the mask is less than the mid-size of the image,
          we perform the desired algorithm*/
        if (radius<Dx/2&&radius<Dy/2) {
            /*Convolve in the x direction*/
            for (nx=0; nx<Dx; nx++) {
                for (ny=0; ny<Dy; ny++) {
                    sum=0.0f;
                    for (t=-radius; t<=radius; t++) {
                        if (nx+t>=Dx)
                            sum = sum+Mask[2*Dx-2-nx-t+ny*Dx]*Kernel[t+radius];
                        else
                            sum = sum+Mask[abs(nx+t)+ny*Dx]*Kernel[t+radius];
                    }
                    Mask2[nx+ny*Dx]=(float) sum;
                }
            }

            /*Convolve in the y direction*/
            for (ny=0; ny<Dy; ny++)
                for (nx=0; nx<Dx; nx++) {
                    sum=0.0f;
                    for (t=-radius; t<=radius; t++) {
                        if (ny+t>=Dy)
                            sum = sum+Mask2[nx+(2*Dy-2-ny-t)*Dx]*
                                      Kernel[t+radius];
                        else
                            sum = sum+Mask2[nx+abs(ny+t)*Dx]*
                                      Kernel[t+radius];
                    }
                    Mask[nx+ny*Dx]=(float) sum;
                }
        } else {
            /*Otherwise, perform a simple gamma correction*/
            for (ny=0; ny<Dy; ny++) {
                for (nx=0; nx<Dx; nx++) {
                    sum = sum+Mask[nx+ny*Dx];
                }
            }
            sum=sum/(Dx*Dy);
            for (ny=0; ny<Dy; ny++) {
                for (nx=0; nx<Dx; nx++) {
                    Mask[nx+ny*Dx] = sum;
                }
            }
        }
    }


    Rout=(float*) malloc(Dx*Dy*sizeof(float));
    Gout=(float*) malloc(Dx*Dy*sizeof(float));
    Bout=(float*) malloc(Dx*Dy*sizeof(float));

    if (option!=1)
        Iout=(float*) malloc(Dx*Dy*sizeof(float));

    for (n=0; n < Dx*Dy; n++) {
        Mask[n]=(float) (2*Mask[n]-1);
        Mask[n]=(float) pow(2.0,(double)Mask[n]);

        if (option==1) {
            Rout[n]=(float) pow((double)R[n],(double)Mask[n]);
            Gout[n]=(float) pow((double)G[n],(double)Mask[n]);
            Bout[n]=(float) pow((double)B[n],(double)Mask[n]);

        } else {
            Iout[n]=(float) pow((double)I[n],(double)Mask[n]);
            switch (option) {
            case 2:
                Iout[n]=(float) pow((double)I[n],(double)Mask[n]);
                if (I[n]==0)
                    r=0;
                else
                    r=Iout[n]/I[n];

                rr=r*(float) R[n];
                gg=r*(float) G[n];
                bb=r*(float) B[n];
                if ((rr >= 1) || (gg >= 1) || (bb >= 1)) {
                    if ((R[n] >= G[n]) && (R[n] >= B[n])) {
                        r=1/(float) R[n];
                    }
                    if ((G[n] >= R[n]) && (G[n] >= B[n])) {
                        r=1/(float) G[n];
                    }
                    if ((B[n] >= R[n]) && (B[n] >= G[n])) {
                        r=1/(float) B[n];
                    }
                    rr=r*(float) R[n];
                    gg=r*(float) G[n];
                    bb=r*(float) B[n];
                }
                Rout[n]=rr;
                Gout[n]=gg;
                Bout[n]=bb;
                break;

            case 3:
                Ypbpr2Rgb(&rr, &gg, &bb, Iout[n], Pb[n], Pr[n]);
                Rout[n]=rr;
                Gout[n]=gg;
                Bout[n]=bb;
                break;

            case 4:
                Hsl2Rgb(&rr, &gg, &bb, H[n], S[n], Iout[n]);
                Rout[n]=rr;
                Gout[n]=gg;
                Bout[n]=bb;
                break;
            }
        }
    }


    output=(unsigned char *) malloc(3*Dx*Dy);
    RGB2output(Rout, Gout, Bout, output, Dx*Dy);

    /* write the PNG output image */
    nc=3; /*color image*/
    write_png_u8(nameout, output, Dx, Dy, nc);

    free(I) ;
    free(Mask) ;
    free(Mask2) ;
    free(Kernel) ;
    free(R);
    free(G);
    free(B);
    free(Rout);
    free(Gout);
    free(Bout);

    if (H) free(H);
    if (S) free(S);
    if (Pb) free(Pb);
    if (Pr) free(Pr);
    if (Iout) free(Iout);

    free(img);
    free(output);

    return EXIT_SUCCESS;
}

