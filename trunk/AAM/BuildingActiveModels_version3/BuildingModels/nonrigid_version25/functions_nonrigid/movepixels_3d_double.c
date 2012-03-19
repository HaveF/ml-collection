#include "mex.h"
#include "math.h"
#include "image_interpolation.h"
#include "multiple_os_thread.h"

/*  This function movepixels, will translate the pixels of an image
 *  according to x, y and z translation images (bilinear interpolated). 
 * 
 *  Iout = movepixels_3d_double(I,Tx,Ty,Tz);
 *
 *  Function is written by D.Kroon University of Twente (July 2009)
 */

voidthread transformvolume(double **Args) {
    /* I is the input image, Iout the transformed image  */
    /* Tx and Ty images of the translation of every pixel. */
    double *Iin, *Iout, *Tx, *Ty, *Tz;
    double *Nthreadsd;
    int Nthreads;
	/*  if one outside pixels are set to zero. */
	double  *moded;
	int mode=0;
   
   	/* Output image size */
    double *ImageSize_d;
    int ImageSize[3];
	
    /* Cubic and outside black booleans */
    bool black, cubic;
    
    /* 3D index storage*/
    int indexI;
    
    /* Size of input image */
    double *Isize_d;
    int Isize[3]={0,0,0};
        
    /* Location of translated pixel */
    double Tlocalx;
    double Tlocaly;
    double Tlocalz;
    
    /* offset */
    int ThreadOffset=0;
    
    /* The thread ID number/name */
    double *ThreadID;
    
    /* X,Y coordinates of current pixel */
    int x,y,z;
    
    Iin=Args[0];
    Iout=Args[1];
    Tx=Args[2];
    Ty=Args[3];
    Tz=Args[4];
    Isize_d=Args[5];
    ThreadID=Args[6];
	moded=Args[7]; mode=(int) moded[0];
	Nthreadsd=Args[8];  Nthreads=(int)Nthreadsd[0];
	ImageSize_d=Args[9];
	
    if(mode==0||mode==2){ black = false; } else { black = true; }
    if(mode==0||mode==1){ cubic = false; } else { cubic = true; }
    
    Isize[0] = (int)Isize_d[0]; 
    Isize[1] = (int)Isize_d[1]; 
    Isize[2] = (int)Isize_d[2]; 
       
	ImageSize[0] = (int)ImageSize_d[0]; 
    ImageSize[1] = (int)ImageSize_d[1]; 
    ImageSize[2] = (int)ImageSize_d[2]; 
	
    ThreadOffset=(int) ThreadID[0];
	
    /*  Loop through all image pixel coordinates */
    for (z=ThreadOffset; z<ImageSize[2]; z=z+Nthreads)
	{
        for (y=0; y<Isize[1]; y++)
        {
            for (x=0; x<Isize[0]; x++)
            {
				indexI=mindex3(x,y,z,ImageSize[0],ImageSize[1]);
				Tlocalx=((double)x)+Tx[indexI];
                Tlocaly=((double)y)+Ty[indexI];
                Tlocalz=((double)z)+Tz[indexI];
                
                /* Set the current pixel value */
                Iout[indexI]=interpolate_3d_double_gray(Tlocalx, Tlocaly, Tlocalz, Isize, Iin,cubic,black); 
            }
        }
    }

    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
	EndThread;
}


/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* I is the input image, Iout the transformed image  */
    /* Tx and Ty images of the translation of every pixel. */
    double *Iin, *Iout, *Tx, *Ty, *Tz, *ImageSizeT;
	double *moded;
    mxArray *matlabCallOut[1]={0};
    mxArray *matlabCallIn[1]={0};
    double *Nthreadsd;
    int Nthreads;
    double *Tlocalx, *Tlocaly, *Tlocalz;
    int x,y,z;
    int index;
    
	
    /* double pointer array to store all needed function variables) */
    double ***ThreadArgs;
    double **ThreadArgs1;
    
	/* Handles to the worker threads */
		ThreadHANDLE *ThreadList;

    
    /* ID of Threads */
    double **ThreadID;              
    double *ThreadID1;

    /* Size of input image */
    int Isize[3]={1,1,1};
    double Isize_d[3]={0,0,0};
    const mwSize *dims;
 
    /* Size of output image */
    double ImageSize_d[3]={0,0,0};
    int ImageSize[3]={1,1,1};
	
	/* Loop variable  */
	int i;
	
    /* Check for proper number of arguments. */
    if(nrhs<5) {
      mexErrMsgTxt("Five inputs are required.");
    } else if(nlhs!=1) {
      mexErrMsgTxt("One output required");
    }
 
    /* Get the sizes of the input image */
    dims = mxGetDimensions(prhs[0]);   
    Isize[0] = (mwSize)dims[0]; 
    Isize[1] = (mwSize)dims[1];
    Isize[2] = (mwSize)dims[2];

    Isize_d[0]=Isize[0];  Isize_d[1]=Isize[1]; Isize_d[2]=Isize[2];
    

    /* Assign pointers to each input. */
    Iin=(double *)mxGetData(prhs[0]);
    Tx=(double *)mxGetData(prhs[1]);
    Ty=(double *)mxGetData(prhs[2]);
    Tz=(double *)mxGetData(prhs[3]);
	moded=(double *)mxGetData(prhs[4]);
	if(nrhs==6)
    {
      ImageSizeT=mxGetPr(prhs[5]);
    }
	
	if(nrhs==5)
	{
		ImageSize_d[0]=ImageSizeT[0];
		ImageSize_d[1]=ImageSizeT[1];
		ImageSize_d[2]=ImageSizeT[2];
	}
	else
	{
		ImageSize_d[0]= Isize_d[0];
		ImageSize_d[1]= Isize_d[1];
		ImageSize_d[2]= Isize_d[2];
	}
	ImageSize[0]=(int)ImageSize_d[0];
	ImageSize[1]=(int)ImageSize_d[1];
	ImageSize[2]=(int)ImageSize_d[2];
	
	/* Create image matrix for the return arguments with the size of input image  */  
    plhs[0] = mxCreateNumericArray(3, ImageSize, mxDOUBLE_CLASS, mxREAL); 
    
	
    /* Assign pointer to output. */
    Iout = (double *)mxGetData(plhs[0]);

    
    if(moded[0]==4)
    {
        Tlocalx=(double*)malloc(Isize[0]*Isize[1]*Isize[2]*sizeof(double));
        Tlocaly=(double*)malloc(Isize[0]*Isize[1]*Isize[2]*sizeof(double));
        Tlocalz=(double*)malloc(Isize[0]*Isize[1]*Isize[2]*sizeof(double));

        for (z=0; z<Isize[2]; z++)
        {
            for (y=0; y<Isize[1]; y++)
            {
                for (x=0; x<Isize[0]; x++)
                {
                    index=x+y*Isize[0]+z*Isize[0]*Isize[1];
                    Tlocalx[index]=((double)x)+Tx[index];
                    Tlocaly[index]=((double)y)+Ty[index];
                    Tlocalz[index]=((double)z)+Tz[index];
                }
            }
        }
        interpolate_forward_3d_double(Iin, Tlocalx, Tlocaly, Tlocalz, Isize, ImageSize, Iout);
        
        free(Tlocalx);
        free(Tlocaly);
        free(Tlocalz);
    }
    else
    {   
        
        mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
        Nthreadsd=mxGetPr(matlabCallOut[0]);
        Nthreads=(int)Nthreadsd[0];
        /* Reserve room for handles of threads in ThreadList  */
            ThreadList = (ThreadHANDLE*)malloc(Nthreads* sizeof( ThreadHANDLE ));

        ThreadID = (double **)malloc( Nthreads* sizeof(double *) );
        ThreadArgs = (double ***)malloc( Nthreads* sizeof(double **) );


      for (i=0; i<Nthreads; i++)
      {
        /*  Make Thread ID  */
        ThreadID1= (double *)malloc( 1* sizeof(double) );
        ThreadID1[0]=i;
        ThreadID[i]=ThreadID1;  

        /*  Make Thread Structure  */
        ThreadArgs1 = (double **)malloc( 10* sizeof( double * ) );  
        ThreadArgs1[0]=Iin;
        ThreadArgs1[1]=Iout;
        ThreadArgs1[2]=Tx;
        ThreadArgs1[3]=Ty;
        ThreadArgs1[4]=Tz;
        ThreadArgs1[5]=Isize_d;
        ThreadArgs1[6]=ThreadID[i];
        ThreadArgs1[7]=moded;
        ThreadArgs1[8]=Nthreadsd;
		ThreadArgs1[9]=ImageSize_d;
	
        /* Start a Thread  */
        ThreadArgs[i]=ThreadArgs1;
        StartThread(ThreadList[i], &transformvolume, ThreadArgs[i])

      }

      for (i=0; i<Nthreads; i++) { WaitForThreadFinish(ThreadList[i]); }

      for (i=0; i<Nthreads; i++) 
      { 
        free(ThreadArgs[i]);
        free(ThreadID[i]);
      }

      free(ThreadArgs);
      free(ThreadID );
      free(ThreadList);
    }
}
        

