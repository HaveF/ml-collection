#include "mex.h"
#include "math.h"
#include "image_interpolation.h"
#include "multiple_os_thread.h"


/*
 *  This function is used to transform an 2D image, in a forwards way with an
 *  transformation image.
 *
 *    Iout = image_interpolation_forward(Iin,Tlocalx,Tlocaly,ImageSize)
 *
 *  inputs,
 * 	   Iin : 2D greyscale or color input image
 * 	   Tlocalx,Tlocaly : (Forwards) Transformation images for all image pixels
 *
 * 	(optional)
 * 	   ImageSize:    - Size of output image
 *  outputs,
 *   	   Iout : The transformed image
 *
 *  Function is written by D.Kroon University of Twente (September 2010)
 *
 *  *  Example
 *    I=im2double(imread('d:\matlab\lena.jpg'));
 *    [x,y]=ndgrid(1:size(I,1),1:size(I,2));
 *    ImageSize=[256 256];
 *    Tlocalx=(x-128)*2-(y-128)*2;
 *    Tlocaly=(y-128)*2+(x-128)*2;
 *    J=image_interpolation_forward(I,Tlocalx,Tlocaly,ImageSize);
 *    figure, imshow(J);
 *
 */


/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    double *Iin, *Tlocalx, *Tlocaly,  *Tlocalz, *ImageSized, *Iout;
    double *Tx, *Ty, *Tz;
    int i;
    
    /*  Size of input image */
    const mwSize *Idims;
    int Isize[3]={1,1,1};
    int ImageSize[3]={1,1,1};
    
    /* Check for proper number of arguments. */
    if(nrhs<4) {
        mexErrMsgTxt("Four or Five inputs are required.");
    } else if(nlhs!=1) {
        mexErrMsgTxt("One output is required");
    }
    
    /* Assign pointers to each input. */
    Iin=mxGetPr(prhs[0]);
    Tx=mxGetPr(prhs[1]);
    Ty=mxGetPr(prhs[2]);
	Tz=mxGetPr(prhs[3]);

    Idims = mxGetDimensions(prhs[0]);
    Isize[0]=Idims[0]; Isize[1]=Idims[1]; Isize[2]=Idims[2];
    
        
    Tlocalx=(double*)malloc(Isize[0]*Isize[1]*Isize[2]*sizeof(double));
    Tlocaly=(double*)malloc(Isize[0]*Isize[1]*Isize[2]*sizeof(double));
    Tlocalz=(double*)malloc(Isize[0]*Isize[1]*Isize[2]*sizeof(double));
    
    for(i=0; i<Isize[0]*Isize[1]*Isize[2]; i++)
    {
        Tlocalx[i]=Tx[i]-1; Tlocaly[i]=Ty[i]-1; Tlocalz[i]=Tz[i]-1;
    }
    
    if(nrhs>4)
    {
        ImageSized=mxGetPr(prhs[4]);
        ImageSize[0]=(int)ImageSized[0];
        ImageSize[1]=(int)ImageSized[1];
        ImageSize[2]=(int)ImageSized[2];
    }
    else
    {
        ImageSize[0]=Isize[0];
        ImageSize[1]=Isize[1];
        ImageSize[2]=Isize[2];
    }
    
    plhs[0] = mxCreateNumericArray(3, ImageSize, mxDOUBLE_CLASS, mxREAL);
    
    /* Assign pointer to output. */
    Iout = (double*)mxGetPr(plhs[0]);
    interpolate_forward_3d_double(Iin, Tlocalx, Tlocaly, Tlocalz, Isize, ImageSize, Iout);
    
    free(Tlocalx);
    free(Tlocaly);
    free(Tlocalz);
}






