#include "mex.h"
#include "math.h"
#include "image_interpolation.h"
#include "multiple_os_thread.h"
            
voidthread bspline_transform_points(double **Args) 
{
    double *X, *O_trans, *Spacing, *Osized, *Xsized, *Tlocal;
    int indexx, indexy;
    double V[4], U[4];
    double v, u;
    double p, p2, p3;
    double x2, y2;
    int i;
    int ix, iy;
    int ixs, iys;
    int m, l;
    double Cx, Cy;
    double Tx, Ty;
    int ind;
    int inds;
    int Osize[4]={1,1,1,1};
    int Xsize[2]={1,1};
    double *Nthreadsd;
    double *ThreadID;
    int Nthreads;
    int ThreadOffset;
    
    X=Args[0];
    O_trans=Args[1];
    Spacing=Args[2];
    Osized =Args[3];
    Xsized =Args[4];
    Tlocal =Args[5];
    ThreadID=Args[6];
    Nthreadsd=Args[7]; 
    
    Nthreads=(int)Nthreadsd[0];
    ThreadOffset=(int) ThreadID[0];
        
    Osize[0]=(int)Osized[0];
    Osize[1]=(int)Osized[1];
    Osize[2]=(int)Osized[2];
 
    Xsize[0]=(int)Xsized[0];
    Xsize[1]=(int)Xsized[1];
    
    inds=Osize[0]*Osize[1];
    
    /* Make row vectors of input coordinates */
    for (i=ThreadOffset; i<Xsize[0]; i=i+Nthreads)
    {
        
        indexx=i;
        indexy=i+Xsize[0];
        x2=X[indexx]; 
        y2=X[indexy];
        
        ixs=(int)floor(x2/Spacing[0]);
        iys=(int)floor(y2/Spacing[1]);
                            
        /* Calculate the b-spline interpolation constants u,v in the center cell */
        /* range between 0 and 1 */
        v  = (x2-ixs*Spacing[0])/Spacing[0];
        u  = (y2-iys*Spacing[1])/Spacing[1];
       
        /* Get the b-spline coefficients in a matrix V,U,W which contains */
        /* the influence of all knots on the points in (x2,y2) */
        p2=v*v; p3=v*p2;
        V[0] = pow3(1-v)/6; V[1] = ( 3*p3 - 6*p2 + 4)/6;
        V[2] = (-3*p3 + 3*p2 + 3*v + 1)/6; V[3] =  p3/6;
        p2=u*u; p3=u*p2;
        U[0] = pow3(1-u)/6; U[1] = ( 3*p3 - 6*p2 + 4)/6;
        U[2] = (-3*p3 + 3*p2 + 3*u + 1)/6; U[3] =  p3/6;
               
        /* This code calculates for every coordinate in X, the indices of all */
        /* b-spline knots which have influence on the transformation value of */
        /* this point */
        Tx=0;
        Ty=0;
        for(m=0; m<4; m++)
        {
            for(l=0; l<4; l++)
            {
                ix=ixs+m;
                iy=iys+l;
                ix=max(min(ix,Osize[0]-1),0);
                iy=max(min(iy,Osize[1]-1),0);

                /* Look up the b-spline knot values in neighborhood of the points in (x2,y2) */
                ind=ix+iy*Osize[0];
                Cx=O_trans[ind];
                Cy=O_trans[ind+inds];

                /* Calculate the transformation of the points in (x2,y2) by the b-spline grid */
                p=V[m]*U[l];
                Tx+=p*Cx;
                Ty+=p*Cy;
            }
        }
        Tlocal[indexx]=Tx; 
        Tlocal[indexy]=Ty;
    }

     /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
	EndThread;
}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *O_trans, *Spacing, *X, *Tlocal;
    const mwSize *dims;
    double Osized[3]={1,1,1};
    double Xsized[2]={1,1};
    mxArray *matlabCallOut[1]={0};
    mxArray *matlabCallIn[1]={0};
    double *Nthreadsd;
    int Nthreads;
    int i;
    
    /* double pointer array to store all needed function variables) */
    double ***ThreadArgs;
    double **ThreadArgs1;
    
	/* Handles to the worker threads */
	ThreadHANDLE *ThreadList;

    /* ID of Threads */
    double **ThreadID;              
    double *ThreadID1;
            
    O_trans=(double *)mxGetData(prhs[0]);
    Spacing=(double *)mxGetData(prhs[1]);
    X=(double *)mxGetData(prhs[2]);
   
    /* Get the sizes of the grid */
    dims = mxGetDimensions(prhs[0]);   
    Osized[0] = dims[0]; 
    Osized[1] = dims[1];
    Osized[2] = dims[2];

    /* Get the sizes of the points which will be warped*/
    dims = mxGetDimensions(prhs[2]);   
    Xsized[0] = dims[0]; 
    Xsized[1] = dims[1];
   
    /* Create an output aray for the warped points */
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL); 
    Tlocal=(double *)mxGetData(plhs[0]);
    
    /* Get number of allowed threads */
    mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
    Nthreadsd=mxGetPr(matlabCallOut[0]);
    Nthreads=(int)Nthreadsd[0];
    /* Reserve room for handles of threads in ThreadList  */

    ThreadList = (ThreadHANDLE*)malloc(Nthreads* sizeof( ThreadHANDLE ));

    ThreadID = (double **)malloc( Nthreads* sizeof(double *) );
    ThreadArgs = (double ***)malloc( Nthreads* sizeof(double **) );
    
    /* Reserve room for 16 function variables(arrays)   */
    for (i=0; i<Nthreads; i++)
    {
        /*  Make Thread ID  */
        ThreadID1= (double *)malloc( 1* sizeof(double) );
        ThreadID1[0]=i;
        ThreadID[i]=ThreadID1;  

        /*  Make Thread Structure  */
        ThreadArgs1 = (double **)malloc( 8* sizeof( double * ) ); 
        ThreadArgs1[0]=X;
        ThreadArgs1[1]=O_trans;
        ThreadArgs1[2]=Spacing;
        ThreadArgs1[3]=Osized;
        ThreadArgs1[4]=Xsized; 
        ThreadArgs1[5]=Tlocal; 
        ThreadArgs1[6]=ThreadID[i];
        ThreadArgs1[7]=Nthreadsd;
        ThreadArgs[i]=ThreadArgs1;
            
        StartThread(ThreadList[i], &bspline_transform_points, ThreadArgs[i])
    }

    for (i=0; i<Nthreads; i++) { WaitForThreadFinish(ThreadList[i]); }

    for (i=0; i<Nthreads; i++) { free(ThreadArgs[i]); free(ThreadID[i]); }
    free(ThreadArgs);
    free(ThreadID );
    free(ThreadList);
   
}
