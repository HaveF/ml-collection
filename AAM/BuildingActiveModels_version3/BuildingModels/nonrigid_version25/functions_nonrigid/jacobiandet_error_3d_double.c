#include "mex.h"
#include "math.h"
#include "image_interpolation.h"
#include "multiple_os_thread.h"
#include "bspline_functions.h"

double * doublemalloc(int a)
{
   void *ptr;
   ptr = malloc(a*sizeof(double));
   if (ptr == NULL) {
        mexErrMsgTxt("Not enough memory");
        abort();
        return 0;  
   }
   return (double*)ptr;
}

double deterror(double Idet)
{
    if(Idet<0.01) return 5.6002032713-99.5033085317*Idet; 
    if(Idet==0) return 0; 
    return fabs(log(Idet)); 
}


voidthread jacobian_error(double **Args) {
    double *Bu, *Bv, *Bw, *Bdu, *Bdv, *Bdw, *dxa, *dya, *dza, *ThreadID, *Ox, *Oy, *Oz, *ThreadErrorOut;
    double *Isize_d;
    double *Osize_d;
    int Isize[3]={0, 0, 0};
    int Osize[3]={0, 0, 0};
    int Onumel;
    double *Nthreadsd;
    int Nthreads;
    /* Multiple threads, one does the odd the other even indexes */
    int ThreadOffset;
    /* Location of pixel which will be come the current pixel */
    double Tlocald[9]={0,0,0,0,0,0,0,0,0};
    /* Grid distance */
    int dx, dy, dz;
    /* X,Y coordinates of current pixel */
    int x, y, z;
    /* B-spline variablesl */
    int u_index=0, v_index=0, w_index=0;
    int i, j, k;

    /* Look up tables index */
    int *u_index_array, *i_array;
    int *v_index_array, *j_array;
    int *w_index_array, *k_array;
       /* current accumlated image error */
    double err=0;
    /* Current determinant /error  pixel */
    double Idet;
    
    /* Split input into variables */
    Bu=Args[0];
    Bv=Args[1];
    Bw=Args[2];
    Isize_d=Args[3];
    Osize_d=Args[4];
    ThreadErrorOut=Args[5];
    dxa=Args[6];
    dya=Args[7];
    dza=Args[8];
    ThreadID=Args[9];
    Ox=Args[10];
    Oy=Args[11];
    Oz=Args[12];
    Nthreadsd=Args[13];  Nthreads=(int)Nthreadsd[0];
    Bdu=Args[14];
    Bdv=Args[15];
    Bdw=Args[16];
        
    Isize[0] = (int)Isize_d[0];
    Isize[1] = (int)Isize_d[1];
    Isize[2] = (int)Isize_d[2];
    Osize[0] = (int)Osize_d[0];
    Osize[1] = (int)Osize_d[1];
    Osize[2] = (int)Osize_d[2];
    
    Onumel=Osize[0]*Osize[1]*Osize[2];
    
    /* Get the spacing of the uniform b-spline grid */
    dx=(int)dxa[0]; dy=(int)dya[0]; dz=(int)dza[0];
    
    ThreadOffset=(int) ThreadID[0];
    /*  Calculate the indexs need to look up the B-spline values. */
    u_index_array= (int*)malloc(Isize[0]* sizeof(int));
    i_array= (int*)malloc(Isize[0]* sizeof(int));
    v_index_array= (int*)malloc(Isize[1]* sizeof(int));
    j_array= (int*)malloc(Isize[1]* sizeof(int));
    w_index_array= (int*)malloc(Isize[2]* sizeof(int));
    k_array= (int*)malloc(Isize[2]* sizeof(int));
    for (x=0; x<Isize[0]; x++) {
        u_index_array[x]=(x%dx)*4; /* Already multiplied by 4, because it specifies the y dimension */
        i_array[x]=(int)floor((double)x/dx); /*  (first row outside image against boundary artefacts)  */
    }
    for (y=0; y<Isize[1]; y++) {
        v_index_array[y]=(y%dy)*4; j_array[y]=(int)floor((double)y/dy);
    }
    for (z=ThreadOffset; z<Isize[2]; z++) {
        w_index_array[z]=(z%dz)*4; k_array[z]=(int)floor((double)z/dz);
    }
    
    /*  Loop through all image pixel coordinates */
    for (z=ThreadOffset; z<Isize[2]; z=z+Nthreads) {
        w_index=w_index_array[z]; k=k_array[z];
        for (y=ThreadOffset; y<Isize[1]; y=y+Nthreads) {
            v_index=v_index_array[y]; j=j_array[y];
            for (x=0; x<Isize[0]; x++) {
                u_index=u_index_array[x]; i=i_array[x];
                /*  This part calculates the coordinates of the pixel */
                /*  which will be transformed to the current x,y pixel. */
                LocalTransformation3D(i, j, k, Bu, Bv, Bw, Bdu, Bdv, Bdw, u_index, v_index, w_index, Ox, Oy, Oz, Osize, Tlocald);
                /* Set the current pixel value */
				Idet=Tlocald[0]*Tlocald[4]*Tlocald[8] + Tlocald[1]*Tlocald[5]*Tlocald[6] +  Tlocald[2]*Tlocald[3]*Tlocald[7] - Tlocald[2]*Tlocald[4]*Tlocald[6] -Tlocald[1]*Tlocald[3]*Tlocald[8] - Tlocald[0]*Tlocald[5]*Tlocald[7];
                err+=deterror(Idet);
            }
        }
    }
    ThreadErrorOut[ThreadOffset]=err;
    
    /* Free memory index look up tables */
    free(u_index_array);
    free(i_array);
    free(v_index_array);
    free(j_array);
    free(w_index_array);
    free(k_array);
    
    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    EndThread;
}

voidthread jacobian_errorgradient(double **Args) {
    double *Bu, *Bv, *Bw, *Bdu, *Bdv, *Bdw, *dxa, *dya, *dza, *ThreadID, *Ox, *Oy, *Oz, *ThreadErrorOut;
    double *ThreadGradientOutX, *ThreadGradientOutY, *ThreadGradientOutZ;
	double *Isize_d;
    double *Osize_d;
    int Isize[3]={0, 0, 0};
    int Osize[3]={0, 0, 0};
    int Onumel;
    double *Nthreadsd;
    int Nthreads;
    /* Multiple threads, one does the odd the other even indexes */
    int ThreadOffset;
    /* index offset in grid gradient */
    int IndexOffGrad[64];
    /* The control points which influence a pixel */
    int nCPinvolved[1];
    int CPinvolved[64];
    int *membergrid;
    /* Step finite difference */
    double step=0.001;
    /* Location of pixel which will be come the current pixel */
    double TlocalGradxx[64];
    double TlocalGradxy[64];
    double TlocalGradxz[64];
    double TlocalGradyx[64];
    double TlocalGradyy[64];
    double TlocalGradyz[64];
    double TlocalGradzx[64];
    double TlocalGradzy[64];
    double TlocalGradzz[64];
    double Tlocald[9]={0,0,0,0,0,0,0,0,0};
    
    /* Loop variables */
    int ix, iy, iz, jx, jy, jz, k1, k2, k3;
    /* Variables to store 1D index */
    int indexO;
    /* Grid distance */
    int dx, dy, dz;
    /* X,Y coordinates of current pixel */
    int x, y, z;
    /* B-spline variablesl */
    int u_index=0, v_index=0, w_index=0;
    int i, j, k;
    /* index offsets */
    int offset1;
      /* Look up tables index */
    int *u_index_array, *i_array;
    int *v_index_array, *j_array;
    int *w_index_array, *k_array;
     /* current accumlated image error */
    double err=0;
    double *err_gradientx, *err_gradienty, *err_gradientz;
    /* Current determinant /error  pixel */
    double Idet, Ierr;
    
    /* Split input into variables */
    Bu=Args[0];
    Bv=Args[1];
    Bw=Args[2];
    Isize_d=Args[3];
    Osize_d=Args[4];
    ThreadErrorOut=Args[5];
    dxa=Args[6];
    dya=Args[7];
    dza=Args[8];
    ThreadID=Args[9];
    Ox=Args[10];
    Oy=Args[11];
    Oz=Args[12];
    Nthreadsd=Args[13];  Nthreads=(int)Nthreadsd[0];
    Bdu=Args[14];
    Bdv=Args[15];
    Bdw=Args[16];
    ThreadGradientOutX=Args[17];
    ThreadGradientOutY=Args[18];
	ThreadGradientOutZ=Args[19];
	
	
    Isize[0] = (int)Isize_d[0];
    Isize[1] = (int)Isize_d[1];
    Isize[2] = (int)Isize_d[2];
    Osize[0] = (int)Osize_d[0];
    Osize[1] = (int)Osize_d[1];
    Osize[2] = (int)Osize_d[2];
    Onumel=Osize[0]*Osize[1]*Osize[2];
    
    /* Get the spacing of the uniform b-spline grid */
    dx=(int)dxa[0]; dy=(int)dya[0]; dz=(int)dza[0];
    
    ThreadOffset=(int) ThreadID[0];
    /*  Calculate the indexs need to look up the B-spline values. */
    u_index_array= (int*)malloc(Isize[0]* sizeof(int));
    i_array= (int*)malloc(Isize[0]* sizeof(int));
    v_index_array= (int*)malloc(Isize[1]* sizeof(int));
    j_array= (int*)malloc(Isize[1]* sizeof(int));
    w_index_array= (int*)malloc(Isize[2]* sizeof(int));
    k_array= (int*)malloc(Isize[2]* sizeof(int));
    for (x=0; x<Isize[0]; x++) {
        u_index_array[x]=(x%dx)*4; /* Already multiplied by 4, because it specifies the y dimension */
        i_array[x]=(int)floor((double)x/dx); /*  (first row outside image against boundary artefacts)  */
    }
    for (y=0; y<Isize[1]; y++) {
        v_index_array[y]=(y%dy)*4; j_array[y]=(int)floor((double)y/dy);
    }
    for (z=ThreadOffset; z<Isize[2]; z++) {
        w_index_array[z]=(z%dz)*4; k_array[z]=(int)floor((double)z/dz);
    }
    
    /* Initialize gradient error storage */
    err_gradientx=doublemalloc(Onumel);
    err_gradienty=doublemalloc(Onumel);
    err_gradientz=doublemalloc(Onumel);
    
    for(k1=0; k1<Onumel; k1++)	{ err_gradientx[k1]=0; err_gradienty[k1]=0; err_gradientz[k1]=0; }
    
    /* Make the grids for the finite difference gradient*/
    membergrid=(int*)malloc(Onumel*sizeof(int));
    
    /* Copy the current grid to all gradient grid arrays */
    for(k1=0; k1<64; k1++) {
        IndexOffGrad[k1]=k1*Onumel;
    }
    
    /* Move every 4th node in the grid arrays */
    for (iz=0; iz<4; iz++) {
        for (iy=0; iy<4; iy++) {
            for (ix=0; ix<4; ix++) {
                k2=ix+iy*4+iz*16;
                k3=IndexOffGrad[k2];
                for (jz=iz; jz<Osize[2]; jz+=4) {
                    for (jy=iy; jy<Osize[1]; jy+=4) {
                        for (jx=ix; jx<Osize[0]; jx+=4) {
                            k1=jx+jy*Osize[0]+jz*Osize[0]*Osize[1];
                            membergrid[k1]=k2;
                        }
                    }
                }
            }
        }
    }
    
    
    /*  Loop through all image pixel coordinates */
    for (z=ThreadOffset; z<Isize[2]; z=z+Nthreads) {
        w_index=w_index_array[z]; k=k_array[z];
        for (y=ThreadOffset; y<Isize[1]; y=y+Nthreads) {
            v_index=v_index_array[y]; j=j_array[y];
            for (x=0; x<Isize[0]; x++) {
                /* Calculate the index needed to loop up the B-spline values. */
                u_index=u_index_array[x]; i=i_array[x];
                /*  This part calculates the coordinates of the pixel */
                /*  which will be transformed to the current x,y pixel. */
                LocalTransformation3DGradient( i, j,  k,  Bu,  Bv,  Bw,  Bdu, Bdv,  Bdw,  u_index, v_index,  w_index,  Ox,  Oy,  Oz, Osize,  nCPinvolved,  Tlocald,  TlocalGradxx,  TlocalGradxy, TlocalGradxz,  TlocalGradyx,  TlocalGradyy, TlocalGradyz,  TlocalGradzx,  TlocalGradzy, TlocalGradzz, CPinvolved);
				
				for(k1=0; k1<64; k1++)
				{
					TlocalGradxx[k1]+=Tlocald[0]; 
					TlocalGradxy[k1]+=Tlocald[1];
					TlocalGradxz[k1]+=Tlocald[2];
					TlocalGradyx[k1]+=Tlocald[3]; 
					TlocalGradyy[k1]+=Tlocald[4];
					TlocalGradyz[k1]+=Tlocald[5];
					TlocalGradzx[k1]+=Tlocald[6]; 
					TlocalGradzy[k1]+=Tlocald[7];
					TlocalGradzz[k1]+=Tlocald[8];
				}
			
                /* Set the current pixel value */
                Idet=Tlocald[0]*Tlocald[4]*Tlocald[8] + Tlocald[1]*Tlocald[5]*Tlocald[6] +  Tlocald[2]*Tlocald[3]*Tlocald[7] - Tlocald[2]*Tlocald[4]*Tlocald[6] -Tlocald[1]*Tlocald[3]*Tlocald[8] - Tlocald[0]*Tlocald[5]*Tlocald[7];
                Ierr=deterror(Idet);
                err+=Ierr;
                for(k1=0; k1<nCPinvolved[0]; k1++) {
                    indexO=CPinvolved[k1];
                    k3=membergrid[indexO];
                    
                    Idet=TlocalGradxx[k3]*Tlocald[4]*Tlocald[8] + TlocalGradxy[k3]*Tlocald[5]*Tlocald[6] +  TlocalGradxz[k3]*Tlocald[3]*Tlocald[7] - TlocalGradxz[k3]*Tlocald[4]*Tlocald[6] -TlocalGradxy[k3]*Tlocald[3]*Tlocald[8] - TlocalGradxx[k3]*Tlocald[5]*Tlocald[7];
                    err_gradientx[indexO]+=deterror(Idet)-Ierr;
                    
                    Idet=Tlocald[0]*TlocalGradyy[k3]*Tlocald[8] + Tlocald[1]*TlocalGradyz[k3]*Tlocald[6] +  Tlocald[2]*TlocalGradyx[k3]*Tlocald[7] - Tlocald[2]*TlocalGradyy[k3]*Tlocald[6] -Tlocald[1]*TlocalGradyx[k3]*Tlocald[8] - Tlocald[0]*TlocalGradyz[k3]*Tlocald[7];
                    err_gradienty[indexO]+=deterror(Idet)-Ierr;
                    
                    Idet=Tlocald[0]*Tlocald[4]*TlocalGradzz[k3] + Tlocald[1]*Tlocald[5]*TlocalGradzx[k3] +  Tlocald[2]*Tlocald[3]*TlocalGradzy[k3] - Tlocald[2]*Tlocald[4]*TlocalGradzx[k3] -Tlocald[1]*Tlocald[3]*TlocalGradzz[k3] - Tlocald[0]*Tlocald[5]*TlocalGradzy[k3];
                    err_gradientz[indexO]+=deterror(Idet)-Ierr;
                }
            }
        }
    }
    /* Return error outputs */
  /* Return error outputs */
	ThreadErrorOut[ThreadOffset]=err;
              
	offset1=ThreadOffset*Onumel;
    for(j=0; j<Onumel; j++) {
		ThreadGradientOutX[j+offset1]=err_gradientx[j];
		ThreadGradientOutY[j+offset1]=err_gradienty[j];
		ThreadGradientOutZ[j+offset1]=err_gradientz[j];
	}
	    
    /* Empty arrays made with Malloc */
    free(err_gradientx);
    free(err_gradienty);
    free(err_gradientz);
        
    /* Free memory index look up tables */
    free(u_index_array);
    free(v_index_array);
    free(w_index_array);
    
    free(i_array);
    free(j_array);
    free(k_array);
    
    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    EndThread;
}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] ) {
    /* Ox and Oy are the grid points */
    /* Zo is the input image */
    /* Zi is the transformed image */
    
    /* nx and ny are the number of grid points (inside the image) */
    double *Ox, *Oy, *Oz, *dxa, *dya, *dza, *E, *Egradient, *ThreadErrorOut;
	double *ThreadGradientOutX, *ThreadGradientOutY, *ThreadGradientOutZ;
    mxArray *matlabCallOut[1]={0};
    mxArray *matlabCallIn[1]={0};
    double *Nthreadsd;
    int Nthreads;
    /* Finite difference step size */
    double step=0.001;
    /* index offsets */
    int offset1;
    /* Dims outputs */
    const int dims_error[2]={1, 1};
    int dims_error_gradient[4]={1, 1, 1, 3};
    /* double pointer array to store all needed function variables) */
    double ***ThreadArgs;
    double **ThreadArgs1;
    /* Handles to the worker threads */
	ThreadHANDLE *ThreadList;
            
    /* ID of Threads */
    double **ThreadID;
    double *ThreadID1;
    
    /* Size of input image */
    double *Isize_d;
    const mwSize *dims;
    
    /* Size of grid */
    mwSize  Osizex, Osizey, Osizez;
    int Onumel;
	double Inumel;
    double Osize_d[3]={0, 0, 0};
    
    /* B-spline variablesl */
    double u, v, w;
    int u_index=0;
    int v_index=0;
    int w_index=0;
    
    double *Bu, *Bv, *Bw;
    double *Bdu, *Bdv, *Bdw;
    
    /* Loop variables  */
    int i, j;
    /* Grid distance */
    int dx, dy, dz;
    /* X,Y,Z coordinates of current pixel */
    int x, y, z;
    
    /* Check for proper number of arguments. */
    if(nrhs!=7) {
        mexErrMsgTxt("Seven inputs are required.");
    }
    
    /* Get the sizes of the grid */
    dims = mxGetDimensions(prhs[0]);
    Osizex = dims[0];
    Osizey = dims[1];
    Osizez = dims[2];
    Onumel = Osizex*Osizey*Osizez;
    
    /* Create image matrix for the Error return argument  */
    plhs[0] = mxCreateNumericArray(2, dims_error, mxDOUBLE_CLASS, mxREAL);
    if(nlhs>1) {
        dims_error_gradient[0]=Osizex; dims_error_gradient[1]=Osizey; dims_error_gradient[2]=Osizez;
        /* Error Gradient needed */
        plhs[1] = mxCreateNumericArray(4, dims_error_gradient, mxDOUBLE_CLASS, mxREAL);
    }
    
    /* Assign pointers to each input. */
    Ox=(double *)mxGetData(prhs[0]);
    Oy=(double *)mxGetData(prhs[1]);
    Oz=(double *)mxGetData(prhs[2]);
    Isize_d=(double *)mxGetData(prhs[3]);
    dxa=(double *)mxGetData(prhs[4]);
    dya=(double *)mxGetData(prhs[5]);
    dza=(double *)mxGetData(prhs[6]);
    Inumel=Isize_d[0]*Isize_d[1]*Isize_d[2];
    
    /* Get the spacing of the uniform b-spline grid */
    dx=(int)dxa[0]; dy=(int)dya[0]; dz=(int)dza[0];
    
    /* Get number of allowed threads */
    mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
    Nthreadsd=mxGetPr(matlabCallOut[0]);
    Nthreads=(int)Nthreadsd[0];
    
    /* Reserve room for handles of threads in ThreadList  */
    #ifdef _WIN32
            ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
    #else
            ThreadList = (ThreadHANDLE*)malloc(Nthreads* sizeof( ThreadHANDLE ));
    #endif
            
            ThreadID = (double **)malloc( Nthreads* sizeof(double *) );
    ThreadArgs = (double ***)malloc( Nthreads* sizeof(double **) );
    ThreadErrorOut= (double *)malloc(Nthreads* sizeof(double) ); 
    if(nlhs==1)
    {	
        ThreadGradientOutX=NULL;
        ThreadGradientOutY=NULL;
		ThreadGradientOutZ=NULL;
    }
    else 
    {   
        ThreadGradientOutX= (double *)malloc(Nthreads*Onumel*sizeof(double));
        ThreadGradientOutY= (double *)malloc(Nthreads*Onumel*sizeof(double));
        ThreadGradientOutZ= (double *)malloc(Nthreads*Onumel*sizeof(double));
		}
    /* Assign pointer to output. */
    E = mxGetPr(plhs[0]);
    if(nlhs>1) { Egradient = mxGetPr(plhs[1]); }
    
    /*  Make polynomial look up tables   */
    Bu=malloc(dx*4*sizeof(double));
    Bv=malloc(dy*4*sizeof(double));
    Bw=malloc(dz*4*sizeof(double));
    Bdu=malloc(dx*4*sizeof(double));
    Bdv=malloc(dy*4*sizeof(double));
    Bdw=malloc(dz*4*sizeof(double));
    
    for (x=0; x<dx; x++) {
        u=((double)x/(double)dx)-floor((double)x/(double)dx);
        Bu[mindex2(0, x, 4)] = BsplineCoefficient(u, 0);
        Bu[mindex2(1, x, 4)] = BsplineCoefficient(u, 1);
        Bu[mindex2(2, x, 4)] = BsplineCoefficient(u, 2);
        Bu[mindex2(3, x, 4)] = BsplineCoefficient(u, 3);
        Bdu[mindex2(0, x, 4)] = BsplineCoefficientDerivative(u, 0)/dxa[0];
        Bdu[mindex2(1, x, 4)] = BsplineCoefficientDerivative(u, 1)/dxa[0];
        Bdu[mindex2(2, x, 4)] = BsplineCoefficientDerivative(u, 2)/dxa[0];
        Bdu[mindex2(3, x, 4)] = BsplineCoefficientDerivative(u, 3)/dxa[0];
    }
    
    for (y=0; y<dy; y++) {
        v=((double)y/(double)dy)-floor((double)y/(double)dy);
        Bv[mindex2(0, y, 4)] = BsplineCoefficient(v, 0);
        Bv[mindex2(1, y, 4)] = BsplineCoefficient(v, 1);
        Bv[mindex2(2, y, 4)] = BsplineCoefficient(v, 2);
        Bv[mindex2(3, y, 4)] = BsplineCoefficient(v, 3);
        Bdv[mindex2(0, y, 4)] = BsplineCoefficientDerivative(v, 0)/dya[0];
        Bdv[mindex2(1, y, 4)] = BsplineCoefficientDerivative(v, 1)/dya[0];
        Bdv[mindex2(2, y, 4)] = BsplineCoefficientDerivative(v, 2)/dya[0];
        Bdv[mindex2(3, y, 4)] = BsplineCoefficientDerivative(v, 3)/dya[0];
    }
    
    for (z=0; z<dz; z++) {
        w=((double)z/(double)dz)-floor((double)z/(double)dz);
        Bw[mindex2(0, z, 4)] = BsplineCoefficient(w, 0);
        Bw[mindex2(1, z, 4)] = BsplineCoefficient(w, 1);
        Bw[mindex2(2, z, 4)] = BsplineCoefficient(w, 2);
        Bw[mindex2(3, z, 4)] = BsplineCoefficient(w, 3);
        Bdw[mindex2(0, z, 4)] = BsplineCoefficientDerivative(w, 0)/dza[0];
        Bdw[mindex2(1, z, 4)] = BsplineCoefficientDerivative(w, 1)/dza[0];
        Bdw[mindex2(2, z, 4)] = BsplineCoefficientDerivative(w, 2)/dza[0];
        Bdw[mindex2(3, z, 4)] = BsplineCoefficientDerivative(w, 3)/dza[0];
        
    }
    
    
    Osize_d[0]=(double)Osizex;  Osize_d[1]=(double)Osizey; Osize_d[2]=(double)Osizez;
    
    /* Reserve room for 16 function variables(arrays)   */
    for (i=0; i<Nthreads; i++) {
        /*  Make Thread ID  */
        ThreadID1= (double *)malloc( 1* sizeof(double) );
        ThreadID1[0]=(double)i;
        ThreadID[i]=ThreadID1;
        
        /*  Make Thread Structure  */
        ThreadArgs1 = (double **)malloc( 20 * sizeof( double * ) );
        ThreadArgs1[0]=Bu;
        ThreadArgs1[1]=Bv;
        ThreadArgs1[2]=Bw;
        ThreadArgs1[3]=Isize_d;
        ThreadArgs1[4]=Osize_d;
        ThreadArgs1[5]=ThreadErrorOut;
        ThreadArgs1[6]=dxa;
        ThreadArgs1[7]=dya;
        ThreadArgs1[8]=dza;
        ThreadArgs1[9]=ThreadID[i];
        ThreadArgs1[10]=Ox;
        ThreadArgs1[11]=Oy;
        ThreadArgs1[12]=Oz;
        ThreadArgs1[13]=Nthreadsd;
        ThreadArgs1[14]=Bdu;
        ThreadArgs1[15]=Bdv;
        ThreadArgs1[16]=Bdw;
        ThreadArgs1[17]=ThreadGradientOutX;
        ThreadArgs1[18]=ThreadGradientOutY;
        ThreadArgs1[19]=ThreadGradientOutZ;
	
        ThreadArgs[i]=ThreadArgs1;
        
        if(nlhs>1) {
            StartThread(ThreadList[i], &jacobian_errorgradient, ThreadArgs[i])
        }
        else {
            StartThread(ThreadList[i], &jacobian_error, ThreadArgs[i])
        }
    }
    
    for (i=0; i<Nthreads; i++) { WaitForThreadFinish(ThreadList[i]); }
    
    /* Add accumlated error of all threads */
	for (i=0; i<Nthreads; i++) 
    {
        E[0]+=ThreadErrorOut[i]; 
    } 
    E[0]/=Inumel;

    if(nlhs>1) {
        for (i=0; i<Nthreads; i++) {
            offset1=i*Onumel;
            for(j=0; j<Onumel; j++) {
                Egradient[j]+=ThreadGradientOutX[j+offset1];
                Egradient[j+Onumel]+=ThreadGradientOutY[j+offset1];
				Egradient[j+Onumel+Onumel]+=ThreadGradientOutZ[j+offset1];
            }
        }
        for(j=0; j<Onumel; j++) {
            Egradient[j]/=Inumel*step;
            Egradient[j+Onumel]/=Inumel*step;
			Egradient[j+Onumel+Onumel]/=Inumel*step;
        }
    }
    	    
    free(ThreadErrorOut);
    free(ThreadGradientOutX);
    free(ThreadGradientOutY);
	free(ThreadGradientOutZ);
    
    free(ThreadArgs);
    free(ThreadID );
    free(ThreadList);
    
    free(Bu);
    free(Bv);
    free(Bw);
    free(Bdu);
    free(Bdv);
    free(Bdw);
}


