#include "mex.h"
#include "math.h"

/*
 * CT slice to radial slice (dental scan like)
 *
 *
 * slicerad=slice2_rad_double(slice)
 *
 *
 */

// Convert 2D/3D matrix index to 1D index
int mindex2(int x, int y, int sizx, int sizy) {
    if(x<0) { x=0; }
    if(y<0) { y=0; }
    if(x>(sizx-1)) { x=sizx-1; }
    if(y>(sizy-1)) { y=sizy-1; }
    
    return y*sizx+x;
}

double pow2(double val){ return val*val; }
double pow3(double val){ return val*val*val; }
double pow4(double val){ return pow2(val)*pow2(val); }


// The matlab mex function
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    // Ox and Oy are the grid points
    // Zo is the input image
    // Zi is the transformed image
    // nx and ny are the number of grid points (inside the image)
    double *Iin, *x1, *y1, *x2, *y2, *steps, *Iout;
    
    double Pi=3.14159265358979323846;
    
    // Dimensions of in and outputs
    const mwSize *Iin_dims;
    const mwSize *x1_dims;
    mwSize Iout_dims[3]={1, 1, 1};
    
    // Center of the image
    double center[2]={0, 0};
    
    // Position of current pixel
    double Tlocalx, Tlocaly;
    double current_angle=0;
    
    // Floor of coordinate
    double fTlocalx, fTlocaly;
    // Zero neighbor
    int xBas0, yBas0;
    // The location in between the pixels 0..1
    double tx, ty;
    // Neighbor loccations
    int xn[4], yn[4];
    
    // The vectors
    double vector_tx[4], vector_ty[4];
    double vector_qx[4]; //={0.25,0.25,0.25,0.25};
    double vector_qy[4]; //={0.25,0.25,0.25,0.25};
    double vector_b[4];
    // Interpolated Intensity;
    double Ipixel=0;
    // Temporary value boundary
    int b;
    
    // Percentage pixel position between center and outer position
    double perc=0;
    
    // Loop variables
    int i, j, k,l;
    
    // 2D to 1D index
    int index, indexoff1, indexoff2;
    
    /* Check for proper number of arguments. */
    if(nrhs!=6) {
        mexErrMsgTxt("6 inputs are required.");
    } else if(nlhs!=1) {
        mexErrMsgTxt("One output required");
    }
    
    
    // Get the sizes of the image
    Iin_dims = mxGetDimensions(prhs[0]);
    
    // Get the sizes of the position vectors
    x1_dims = mxGetDimensions(prhs[1]);
    
    /* Assign pointers to input. */
    Iin=mxGetPr(prhs[0]);
    x1=mxGetPr(prhs[1]);
    y1=mxGetPr(prhs[2]);
    x2=mxGetPr(prhs[3]);
    y2=mxGetPr(prhs[4]);
    steps=mxGetPr(prhs[5]);
    
    
    
    Iout_dims[0]=(int)steps[0];
    Iout_dims[1]=(int)x1_dims[0];
    
    // Make output array
    if(mxGetNumberOfDimensions(prhs[0])>2) {
        Iout_dims[2]=Iin_dims[2];
    }
    plhs[0] = mxCreateNumericArray(3, Iout_dims, mxDOUBLE_CLASS, mxREAL);
    
    /* Assign pointer to output. */
    Iout = mxGetPr(plhs[0]);
    
    
    // Loop through all positions
    for (i=0; i<(int)x1_dims[0]; i++) {
        for (j=0; j<(int)steps[0]; j++) {
            perc=((double)j)/steps[0];
            Tlocalx=x1[i]*(1-perc)+x2[i]*perc;
            Tlocaly=y1[i]*(1-perc)+y2[i]*perc;
            
            // Determine of the zero neighbor
            fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
            xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
            
            // Determine the location in between the pixels 0..1
            tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly;
            
            // Determine the t vectors
            vector_tx[0]= 0.5; vector_tx[1]= 0.5*tx; vector_tx[2]= 0.5*pow2(tx); vector_tx[3]= 0.5*pow3(tx);
            vector_ty[0]= 0.5; vector_ty[1]= 0.5*ty; vector_ty[2]= 0.5*pow2(ty); vector_ty[3]= 0.5*pow3(ty);
            
            // t vector multiplied with 4x4 bicubic kernel gives the to q vectors
            vector_qx[0]= -1.0*vector_tx[1]+2.0*vector_tx[2]-1.0*vector_tx[3];
            vector_qx[1]= 2.0*vector_tx[0]-5.0*vector_tx[2]+3.0*vector_tx[3];
            vector_qx[2]= 1.0*vector_tx[1]+4.0*vector_tx[2]-3.0*vector_tx[3];
            vector_qx[3]= -1.0*vector_tx[2]+1.0*vector_tx[3];
            vector_qy[0]= -1.0*vector_ty[1]+2.0*vector_ty[2]-1.0*vector_ty[3];
            vector_qy[1]= 2.0*vector_ty[0]-5.0*vector_ty[2]+3.0*vector_ty[3];
            vector_qy[2]= 1.0*vector_ty[1]+4.0*vector_ty[2]-3.0*vector_ty[3];
            vector_qy[3]= -1.0*vector_ty[2]+1.0*vector_ty[3];
            
            // Determine 1D neighbour coordinates
            xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
            yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
            
            // Clamp to image boundary if outside image
            if(xn[0]<0) { xn[0]=0;if(xn[1]<0) { xn[1]=0;if(xn[2]<0) { xn[2]=0; if(xn[3]<0) { xn[3]=0; }}}}
            if(yn[0]<0) { yn[0]=0;if(yn[1]<0) { yn[1]=0;if(yn[2]<0) { yn[2]=0; if(yn[3]<0) { yn[3]=0; }}}}
            b=Iin_dims[0]-1;
            if(xn[3]>b) { xn[3]=b;if(xn[2]>b) { xn[2]=b;if(xn[1]>b) { xn[1]=b; if(xn[0]>b) { xn[0]=b; }}}}
            b=Iin_dims[1]-1;
            if(yn[3]>b) { yn[3]=b;if(yn[2]>b) { yn[2]=b;if(yn[1]>b) { yn[1]=b; if(yn[0]>b) { yn[0]=b; }}}}
            
            // First do interpolation in the x direction followed by interpolation in the y direction
            index=mindex2(j, i, Iout_dims[0], Iout_dims[1]);
            for (k=0; k<Iout_dims[2]; k++) {
                indexoff1=k*Iout_dims[0]*Iout_dims[1];
                indexoff2=k*Iin_dims[0]*Iin_dims[1];
                Iout[index+indexoff1]=0;
                for(l=0; l<4; l++) {
                    vector_b[l] =vector_qx[0]*Iin[xn[0]+yn[l]*Iin_dims[0]+indexoff2];
                    vector_b[l]+=vector_qx[1]*Iin[xn[1]+yn[l]*Iin_dims[0]+indexoff2];
                    vector_b[l]+=vector_qx[2]*Iin[xn[2]+yn[l]*Iin_dims[0]+indexoff2];
                    vector_b[l]+=vector_qx[3]*Iin[xn[3]+yn[l]*Iin_dims[0]+indexoff2];
                    Iout[index+indexoff1]+= vector_qy[l]*vector_b[l];
                }
            }
            
        }
    }
}















