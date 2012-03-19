

void LocalTransformation3D(int i,int j,int k,double *Bu, double *Bv,  double *Bw, double *Bdu,double *Bdv,  double *Bdw, int u_index,int v_index, int w_index,double *Ox, double *Oy, double *Oz, int *Osize, double *Tlocald)
{
	/*  B-Spline loop variabels */
	int l,m,n;
	/* index in grid */
	int indexO;
	/* temporary value */
	double valx, valy, valz;
	Tlocald[0]=0; Tlocald[1]=0; Tlocald[2]=0;
	Tlocald[3]=0; Tlocald[4]=0; Tlocald[5]=0;
	Tlocald[6]=0; Tlocald[7]=0; Tlocald[8]=0;
	
	for(l=0; l<4; l++) {
		if(((i+l)>=0)&&((i+l)<Osize[0])) {
			for(m=0; m<4; m++) {
				if(((j+m)>=0)&&((j+m)<Osize[1])) {
					for(n=0; n<4; n++) {
						if(((k+n)>=0)&&((k+n)<Osize[2])) {
							indexO=(i+l)+(j+m)*Osize[0]+(k+n)*Osize[0]*Osize[1];
							valx=Bdu[l+u_index]*Bv[m+v_index]*Bw[n+w_index];
							valy=Bu[l+u_index]*Bdv[m+v_index]*Bw[n+w_index];
							valz=Bu[l+u_index]*Bv[m+v_index]*Bdw[n+w_index];
							
							Tlocald[0]+=valx*Ox[indexO];
							Tlocald[1]+=valy*Ox[indexO];
							Tlocald[2]+=valz*Ox[indexO];
							Tlocald[3]+=valx*Oy[indexO];
							Tlocald[4]+=valy*Oy[indexO];
							Tlocald[5]+=valz*Oy[indexO];
							Tlocald[6]+=valx*Oz[indexO];
							Tlocald[7]+=valy*Oz[indexO];
							Tlocald[8]+=valz*Oz[indexO];
						}
					}
				}
			}
		}
	}
}


void LocalTransformation2D(int i,int j,double *Bu, double *Bv, double *Bdu,double *Bdv, int u_index,int v_index, double *Ox, double *Oy, int *Osize, double *Tlocald)
{
    /*  B-Spline loop variabels */
	int l,m;
	 /* index in grid */
	int indexO;
	 /* temporary value */
	double valx, valy;
	Tlocald[0]=0; Tlocald[1]=0; Tlocald[3]=0; Tlocald[2]=0;
	for(l=0; l<4; l++) {
		if(((i+l)>=0)&&((i+l)<Osize[0])) {
			for(m=0; m<4; m++) {
				if(((j+m)>=0)&&((j+m)<Osize[1])) {
					indexO=(i+l)+(j+m)*Osize[0];
					valx=Bdu[l+u_index]*Bv[m+v_index];
					valy=Bu[l+u_index]*Bdv[m+v_index];
					
					Tlocald[0]+=valx*Ox[indexO]; Tlocald[1]+=valy*Ox[indexO];
					Tlocald[2]+=valx*Oy[indexO]; Tlocald[3]+=valy*Oy[indexO];
				}
			}
		}
	}
}

void LocalTransformation2DGradient(int i,int j,double *Bu, double *Bv, double *Bdu,double *Bdv, int u_index,int v_index, double *Ox, double *Oy, int *Osize, int *nCPinvolved, double *Tlocald, double *TlocalGradxx, double *TlocalGradxy,double *TlocalGradyx, double *TlocalGradyy, int *CPinvolved)
{
    /*  B-Spline loop variabels */
	int l,m;
	 /* index in grid */
	int indexO;
	 /* temporary value */
	double valx, valy;
	int ix,iy;
	int k1;
	double step=0.001;
    nCPinvolved[0]=0;
    Tlocald[0]=0; Tlocald[1]=0; Tlocald[3]=0; Tlocald[2]=0;
	for(l=0; l<4; l++) {
		if(((i+l)>=0)&&((i+l)<Osize[0])) {
			for(m=0; m<4; m++) {
				if(((j+m)>=0)&&((j+m)<Osize[1])) {
                    indexO=(i+l)+(j+m)*Osize[0];
					CPinvolved[nCPinvolved[0]]=indexO; nCPinvolved[0]++;
					valx=Bdu[l+u_index]*Bv[m+v_index];
					valy=Bu[l+u_index]*Bdv[m+v_index];
					
					Tlocald[0]+=valx*Ox[indexO]; Tlocald[1]+=valy*Ox[indexO];
					Tlocald[2]+=valx*Oy[indexO]; Tlocald[3]+=valy*Oy[indexO];
					
					ix=(i+l)%4; iy=(j+m)%4;
					k1=ix+iy*4;
					TlocalGradxx[k1]=valx*step; TlocalGradxy[k1]=valy*step; 
					TlocalGradyx[k1]=valx*step; TlocalGradyy[k1]=valy*step;
				}
			}
		}
	}
}


void LocalTransformation3DGradient(int i, int j, int k, double *Bu, double *Bv, double *Bw, double *Bdu,double *Bdv, double *Bdw, int u_index,int v_index, int w_index, double *Ox, double *Oy, double *Oz,int *Osize, int *nCPinvolved, double *Tlocald, double *TlocalGradxx, double *TlocalGradxy,double *TlocalGradxz, double *TlocalGradyx, double *TlocalGradyy,double *TlocalGradyz, double *TlocalGradzx, double *TlocalGradzy,double *TlocalGradzz, int *CPinvolved)
{
    /*  B-Spline loop variabels */
	int l,m,n;
	 /* index in grid */
	int indexO;
	 /* temporary value */
	double valx, valy, valz;
	int ix,iy, iz;
	int k1;
	double step=0.001;
	Tlocald[0]=0; Tlocald[1]=0; Tlocald[2]=0;
    Tlocald[3]=0; Tlocald[4]=0; Tlocald[5]=0;
    Tlocald[6]=0; Tlocald[7]=0; Tlocald[8]=0;
    nCPinvolved[0]=0;
	for(l=0; l<4; l++) {
		if(((i+l)>=0)&&((i+l)<Osize[0])) {
			for(m=0; m<4; m++) {
				if(((j+m)>=0)&&((j+m)<Osize[1])) {
					for(n=0; n<4; n++) {
						if(((k+n)>=0)&&((k+n)<Osize[2])) {
							indexO=(i+l)+(j+m)*Osize[0]+(k+n)*Osize[0]*Osize[1];
							CPinvolved[nCPinvolved[0]]=indexO; nCPinvolved[0]++;
							valx=Bdu[l+u_index]*Bv[m+v_index]*Bw[n+w_index];
							valy=Bu[l+u_index]*Bdv[m+v_index]*Bw[n+w_index];
							valz=Bu[l+u_index]*Bv[m+v_index]*Bdw[n+w_index];
							
							Tlocald[0]+=valx*Ox[indexO];
							Tlocald[1]+=valy*Ox[indexO];
							Tlocald[2]+=valz*Ox[indexO];
							
							Tlocald[3]+=valx*Oy[indexO];
							Tlocald[4]+=valy*Oy[indexO];
							Tlocald[5]+=valz*Oy[indexO];
							
							Tlocald[6]+=valx*Oz[indexO];
							Tlocald[7]+=valy*Oz[indexO];
							Tlocald[8]+=valz*Oz[indexO];

							ix=(i+l)%4;
							iy=(j+m)%4;
							iz=(k+n)%4;
							k1=ix+iy*4+iz*16;
				
							TlocalGradxx[k1]=valx*step;
							TlocalGradxy[k1]=valy*step;
							TlocalGradxz[k1]=valz*step;
								
							TlocalGradyx[k1]=valx*step;
							TlocalGradyy[k1]=valy*step;
							TlocalGradyz[k1]=valz*step;
								
							TlocalGradzx[k1]=valx*step;
							TlocalGradzy[k1]=valy*step;
							TlocalGradzz[k1]=valz*step;

						}
					}
				}
			}
		}
	}
}
