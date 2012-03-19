#include "mex.h"
#include "math.h"
#include "image_interpolation.h"
#include "multiple_os_thread.h"

            
voidthread bspline_transform_points(double **Args) 
{
    double *X, *O_trans, *Spacing, *Osized, *Xsized, *Tlocal;
    int indexx, indexy, indexz;
    double V[4], U[4], W[4];
    double v, u, w;
    double p, p2, p3;
    double x2, y2, z2;
    int i;
    int ix, iy, iz;
    int ixs, iys, izs;
    int m, l, k;
    double Cx, Cy, Cz;
    double Tx, Ty, Tz;
    int ind, inds;
    int Osize[4]={1,1,1,1};
    int Xsize[2]={1,1};
    double *Nthreadsd;
    double *ThreadID;
    int Nthreads;
    int ThreadOffset;
    double tx2,ty2,tz2;
	double d6=1.0/6.0;
    int offset1, offset2;

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
    Osize[3]=(int)Osized[3];

    Xsize[0]=(int)Xsized[0];
    Xsize[1]=(int)Xsized[1];

	offset2=Osize[0]*Osize[1];

    inds=Osize[0]*Osize[1]*Osize[2];
    /* Make row vectors of input coordinates */
    for (i=ThreadOffset; i<Xsize[0]; i=i+Nthreads)
    {
        indexx=i;
        indexy=i+Xsize[0];
        indexz=i+Xsize[0]*2;
        x2=X[indexx]; 
        y2=X[indexy];
        z2=X[indexz];
        
		tx2=x2/Spacing[0];
		ty2=y2/Spacing[1];
		tz2=z2/Spacing[2];

        ixs=(int)floor(tx2);
        iys=(int)floor(ty2);
        izs=(int)floor(tz2);
                            
        /* Calculate the b-spline interpolation constants u,v in the center cell */
        /* range between 0 and 1 */
        v  = tx2-ixs;
        u  = ty2-iys;
        w  = tz2-izs;

        /* Get the b-spline coefficients in a matrix V,U,W which contains */
        /* the influence of all knots on the points in (x2,y2) */
        p2=v*v; p3=v*p2;
        V[0] = pow3(1-v)*d6; V[1] = ( 3*p3 - 6*p2 + 4)*d6;
        V[2] = (-3*p3 + 3*p2 + 3*v + 1)*d6; V[3] =  p3*d6;
        p2=u*u; p3=u*p2;
        U[0] = pow3(1-u)*d6; U[1] = ( 3*p3 - 6*p2 + 4)*d6;
        U[2] = (-3*p3 + 3*p2 + 3*u + 1)*d6; U[3] =  p3*d6;
        p2=w*w; p3=w*p2;
        W[0] = pow3(1-w)*d6; W[1] = ( 3*p3 - 6*p2 + 4)*d6;
        W[2] = (-3*p3 + 3*p2 + 3*w + 1)*d6; W[3] =  p3*d6;

        
        /* This code calculates for every coordinate in X, the indices of all */
        /* b-spline knots which have influence on the transformation value of */
        /* this point */
        Tx=0;
        Ty=0;
        Tz=0;
        for(m=0; m<4; m++)
        {
		    ix=ixs+m;
			ix=max(min(ix,Osize[0]-1),0);
    		for(l=0; l<4; l++)
            {
                iy=iys+l;
    		    iy=max(min(iy,Osize[1]-1),0);
				offset1=ix+iy*Osize[0];
				for(k=0; k<4; k++)
                {
                    iz=izs+k;
                    iz=max(min(iz,Osize[2]-1),0);

                    /* Look up the b-spline knot values in neighborhood of the points in (x2,y2) */
                    ind=offset1+iz*offset2;
                    Cx=O_trans[ind];
                    Cy=O_trans[ind + inds];
                    Cz=O_trans[ind + 2*inds];

                    /* Calculate the transformation of the points in (x2,y2) by the b-spline grid */
                    p=V[m]*U[l]*W[k];
                    Tx+=p*Cx;
                    Ty+=p*Cy;
                    Tz+=p*Cz;
                }
            }
        }
        Tlocal[indexx]=Tx; 
        Tlocal[indexy]=Ty;
        Tlocal[indexz]=Tz;
    }
   
    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
	EndThread;
}
/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

    double *O_trans, *Spacing, *X, *Tlocal;
    const mwSize *dims;
    double Osized[4]={1,1,1,1};
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
            
    /*function Tlocal=bspline_transform_points_3d(O_trans,Spacing,X) */

    O_trans=(double *)mxGetData(prhs[0]);
    Spacing=(double *)mxGetData(prhs[1]);
    X=(double *)mxGetData(prhs[2]);
      
    /* Get the sizes of the grid */
    dims = mxGetDimensions(prhs[0]);   
    Osized[0] = dims[0]; 
    Osized[1] = dims[1];
    Osized[2] = dims[2];
    Osized[3] = dims[3];

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
