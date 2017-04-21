// http://www.caam.rice.edu/~timwar/CAAM210/Flows.html

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern "C"
{
#include "png_util.h"
}

#define dfloat double

#define FLUID 0
#define WALL 1

#define NSPECIES 9

#include "cuda.h"

// loop up 1D array index from 2D node coordinates
__host__ __device__ int idx(int N, int n, int m){
  return n + m*(N+2);
}

void lbmInput(const char *imageFileName,
	      dfloat threshold,
	      int *outN,
	      int *outM,
	      unsigned char **rgb,
	      unsigned char **alpha,
	      int **nodeType){

  int n,m, N,M;

  // read png file
  read_png(imageFileName, &N, &M, rgb, alpha);

  // pad to guarantee space around obstacle and extend the wake
  int Npad = 3*N;
  int Mpad = 2*M;
  //int Mpad = M;

  if(Npad>8192) Npad = 8192;
  if(Mpad>8192) Mpad = 8192;

  // threshold walls based on gray scale
  *nodeType = (int*) calloc((Npad+2)*(Mpad+2), sizeof(int));

  // mark pixels by gray scale intensity
  unsigned char *rgbPad   = (unsigned char*) calloc(3*(Npad+2)*(Mpad+2), sizeof(unsigned char));
  unsigned char *alphaPad = (unsigned char*) calloc((Npad+2)*(Mpad+2),   sizeof(unsigned char));
  int wallCount = 0;
  for(m=1;m<=M;++m){
    for(n=1;n<=N;++n){
      int offset = ((n-1)+(m-1)*N);
      dfloat r = (*rgb)[3*offset+0];
      dfloat g = (*rgb)[3*offset+1];
      dfloat b = (*rgb)[3*offset+2];
      dfloat a = (*alpha) ? (*alpha)[offset]:255;
      // center image in padded region (including halo zone)
      //   int hoffset = N/4, yoffset = M/2;
      int hoffset = N/4, yoffset = M/2;
      int id = idx(Npad,n+hoffset,m+yoffset);

      if(a==0)
	(*nodeType)[id] = FLUID;
      else
	(*nodeType)[id] = WALL*(sqrt(r*r+g*g+b*b)<threshold);

      wallCount += (*nodeType)[id];
      rgbPad[3*id+0] = r;
      rgbPad[3*id+1] = g;
      rgbPad[3*id+2] = b;
      alphaPad[id] = 255;
    }
  }

  for(n=1;n<=Npad;++n){
    (*nodeType)[idx(Npad,n,1)] = WALL;
    (*nodeType)[idx(Npad,n,Mpad)] = WALL;
  }
  
  free(*rgb); free(*alpha);
  *rgb = rgbPad;
  *alpha = alphaPad;

  
  
  printf("wallCount = %d (%g percent of %d x %d nodes)\n", wallCount, 100.*((dfloat)wallCount/((Npad+2)*(Mpad+2))), Npad, Mpad);
  
  *outN = Npad;
  *outM = Mpad;
}

  
void lbmOutput(const char *fname,
	       const int *nodeType,
	       unsigned char *rgb,
	       unsigned char *alpha,
	       const dfloat c,
	       const dfloat dx,
	       int N,
	       int M,
	       const dfloat *f){
  int n,m,s;
  FILE *bah = fopen(fname, "w");

  // compute vorticity
  dfloat *Ux = (dfloat*) calloc((N+2)*(M+2), sizeof(dfloat));
  dfloat *Uy = (dfloat*) calloc((N+2)*(M+2), sizeof(dfloat));

  dfloat fnm[NSPECIES];
  for(m=1;m<=M;++m){
    for(n=1;n<=N;++n){
      int base = idx(N, n, m);
      for(s=0;s<NSPECIES;++s)
	fnm[s] = f[base+s*(N+2)*(M+2)];

      const dfloat rho = fnm[0]+fnm[1]+fnm[2]+fnm[3]+fnm[4]+fnm[5]+fnm[6]+fnm[7]+fnm[8];
      // macroscopic momentum
      Ux[base] = (fnm[1] - fnm[3] + fnm[5] - fnm[6] - fnm[7] + fnm[8])*c/rho;
      Uy[base] = (fnm[2] - fnm[4] + fnm[5] + fnm[6] - fnm[7] - fnm[8])*c/rho;
    }
  }

  
  dfloat plotMin = -4, plotMax = 4;
  for(m=1;m<=M;++m){
    for(n=1;n<=N;++n){
      int id = idx(N,n,m);

      // over write pixels in fluid region
      if(nodeType[id]==FLUID){
	unsigned char r,g,b,a;

	// reconstruct macroscopic density
	dfloat rho = 0;
	for(s=0;s<NSPECIES;++s)
	  rho += f[id+s*(N+2)*(M+2)];
	rho = ((rho-plotMin)/(plotMax-plotMin)); // rescale

	dfloat dUxdy = (Ux[idx(N,n,m+1)]-Ux[idx(N,n,m-1)])/(2.*dx);
	dfloat dUydx = (Uy[idx(N,n+1,m)]-Uy[idx(N,n-1,m)])/(2.*dx);
	
	dfloat curlU = dUydx-dUxdy;
	curlU = ((curlU-plotMin)/(plotMax-plotMin));

#if 0
	r = 255*curlU;
	g = 255*curlU;
	b = 255*curlU;
	a = 255;
#else
	a = 255;
	if(curlU>.55){
	  r = 255*curlU;
	  g = 0;
	  b = 0;
	}
	else if(curlU<.45){
	  r = 0;
	  g = 0;
	  b = 255*curlU;
	}
	else{
	  r = 255;
	  g = 255;
	  b = 255;
	}
#endif
	rgb[idx(N,n,m)*3+0] = r;
	rgb[idx(N,n,m)*3+1] = g;
	rgb[idx(N,n,m)*3+2] = b;
	alpha[idx(N,n,m)] = a;
      }
    }
  }
  
  write_png(bah, N+2, M+2, rgb, alpha);

  fclose(bah);
  free(Ux);
  free(Uy);
}

#if 0
void lbmOutputMicro(const char *fname,
		    const int *nodeType,
		    unsigned char *rgb,
		    unsigned char *alpha,
		    const dfloat c,
		    const dfloat dx,
		    int N,
		    int M,
		    const dfloat *f){
  int n,m,s;
  FILE *bah = fopen(fname, "w");

  unsigned char *rgb3 = (unsigned char*) calloc(3*NSPECIES*(N+2)*(M+2), sizeof(unsigned char));
  unsigned char *alpha3 = (unsigned char*) calloc(3*NSPECIES*(N+2)*(M+2), sizeof(unsigned char));

  int soffsets[NSPECIES];
  soffsets[0] = 0;
  soffsets[1] = 1;
  soffsets[2] = 3*(N+2);
  soffsets[3] = -1;
  soffsets[4] = -3*(N+2);
  soffsets[5] = 3*(N+2)+1;
  soffsets[6] = 3*(N+2)-1;
  soffsets[7] = -3*(N+2)-1;
  soffsets[8] = -3*(N+2)+1;
  
  // compute vorticity
  dfloat plotMin = -1, plotMax = 1;
  for(s=0;s<NSPECIES;++s){
    for(m=1;m<=M;++m){
      for(n=1;n<=N;++n){
	int id = idx(N,n,m);
	
	int soffset = soffsets[s];

	unsigned char r,g,b,a;
	r = rgb[idx(N,n,m)*3+0];
	g = rgb[idx(N,n,m)*3+1];
	b = rgb[idx(N,n,m)*3+2];

	// over write pixels in fluid region
	if(nodeType[id]==FLUID){

	  dfloat fs = f[idx(N,n,m)+s*(N+2)*(M+2)];
	  fs = ((fs-plotMin)/(plotMax-plotMin));
	  
	  a = 255;
	  if(fs>.6){
	    r = 255*fs;
	    g = 0;
	    b = 0;
	  }
	  else if(fs<.4){
	    r = 0;
	    g = 0;
	    b = 255*fs;
	  }
	  else{
	    r = 255;
	    g = 255;
	    b = 255;
	  }
	  int base = 3*n + 9*(N+2)*m+soffset;
	  rgb3[base*3+0] = r;
	  rgb3[base*3+1] = g;
	  rgb3[base*3+2] = b;
	  alpha3[base] = a;
	}
      }
    }
  }
  
  write_png(bah, 3*(N+2), 3*(M+2), rgb3, alpha3);

  fclose(bah);
  free(rgb3);
  free(alpha3);
}
#endif

// weights used to compute equilibrium distribution (post collision)
const dfloat w0 = 4.f/9.f, w1 = 1.f/9.f, w2 = 1.f/9.f, w3 =  1.f/9.f;
const dfloat w4 = 1.f/9.f, w5 = 1.f/36.f, w6 = 1.f/36.f, w7 = 1.f/36.f, w8 = 1.f/36.f;

__host__ __device__ void lbmEquilibrium(const dfloat c,
					const dfloat rho,
					const dfloat Ux, 
					const dfloat Uy,
					dfloat *  feq){

  // resolve macroscopic velocity into lattice particle velocity directions
  const dfloat U2 = Ux*Ux+Uy*Uy;
  
  const dfloat v0 = 0;
  const dfloat v1 = +Ux/c;
  const dfloat v2 = +Uy/c;
  const dfloat v3 = -Ux/c;
  const dfloat v4 = -Uy/c;
  const dfloat v5 =  (+Ux+Uy)/c;
  const dfloat v6 =  (-Ux+Uy)/c;
  const dfloat v7 =  (-Ux-Uy)/c;
  const dfloat v8 =  (+Ux-Uy)/c;
  
  // compute LBM post-collisional 
  feq[0] = rho*w0*(1.f + 3.f*v0 + 4.5f*v0*v0 - 1.5f*U2/(c*c));
  feq[1] = rho*w1*(1.f + 3.f*v1 + 4.5f*v1*v1 - 1.5f*U2/(c*c));
  feq[2] = rho*w2*(1.f + 3.f*v2 + 4.5f*v2*v2 - 1.5f*U2/(c*c));
  feq[3] = rho*w3*(1.f + 3.f*v3 + 4.5f*v3*v3 - 1.5f*U2/(c*c));
  feq[4] = rho*w4*(1.f + 3.f*v4 + 4.5f*v4*v4 - 1.5f*U2/(c*c));
  feq[5] = rho*w5*(1.f + 3.f*v5 + 4.5f*v5*v5 - 1.5f*U2/(c*c));
  feq[6] = rho*w6*(1.f + 3.f*v6 + 4.5f*v6*v6 - 1.5f*U2/(c*c));
  feq[7] = rho*w7*(1.f + 3.f*v7 + 4.5f*v7*v7 - 1.5f*U2/(c*c));
  feq[8] = rho*w8*(1.f + 3.f*v8 + 4.5f*v8*v8 - 1.5f*U2/(c*c));
}

#define TX 32
#define TY 8

// perform lattice streaming and collision steps
__global__ void lbmUpdate(const int N,                  // number of nodes in x
			  const int M,                  // number of nodes in y
			  const dfloat c,                // speed of sound
			  const dfloat * __restrict__ tau,           // relaxation rate
			  const int    * __restrict__ nodeType,      // (N+2) x (M+2) node types 
			  const dfloat * __restrict__ f,             // (N+2) x (M+2) x 9 fields before streaming and collisions
			  dfloat * __restrict__ fnew){               // (N+2) x (M+2) x 9 fields after streaming and collisions
  
  // number of nodes in whole array including halo
  int Nall = (N+2)*(M+2);
  
  // loop over all non-halo nodes in lattice
  int n = 1 + threadIdx.x + blockIdx.x*TX;
  int m = 1 + threadIdx.y + blockIdx.y*TY;

  if(m<M+1 && n<=N+1){

    // physics paramaters
    dfloat tauinv = 1.f/tau[idx(N,n,m)];
    
    // discover type of node (WALL or FLUID)
    const int nt = nodeType[idx(N,n,m)];
    dfloat fnm[NSPECIES];
    
    // OUTFLOW
    if(n==N+1){
      fnm[0] = f[idx(N,n,  m)   + 0*Nall]; // stationary 
      fnm[1] = f[idx(N,n-1,m)   + 1*Nall]; // E bound from W
      fnm[2] = f[idx(N,n,m-1)   + 2*Nall]; // N bound from S
      fnm[3] = f[idx(N,n,m)     + 3*Nall]; // W bound from E
      fnm[4] = f[idx(N,n,m+1)   + 4*Nall]; // S bound from N
      fnm[5] = f[idx(N,n-1,m-1) + 5*Nall]; // NE bound from SW
      fnm[6] = f[idx(N,n,m-1)   + 6*Nall]; // NW bound from SE
      fnm[7] = f[idx(N,n,m+1)   + 7*Nall]; // SW bound from NE
      fnm[8] = f[idx(N,n-1,m+1) + 8*Nall]; // SE bound from NW      
    }
    else if(nt == FLUID){
      fnm[0] = f[idx(N,n,  m)   + 0*Nall]; // stationary 
      fnm[1] = f[idx(N,n-1,m)   + 1*Nall]; // E bound from W
      fnm[2] = f[idx(N,n,m-1)   + 2*Nall]; // N bound from S
      fnm[3] = f[idx(N,n+1,m)   + 3*Nall]; // W bound from E
      fnm[4] = f[idx(N,n,m+1)   + 4*Nall]; // S bound from N
      fnm[5] = f[idx(N,n-1,m-1) + 5*Nall]; // NE bound from SW
      fnm[6] = f[idx(N,n+1,m-1) + 6*Nall]; // NW bound from SE
      fnm[7] = f[idx(N,n+1,m+1) + 7*Nall]; // SW bound from NE
      fnm[8] = f[idx(N,n-1,m+1) + 8*Nall]; // SE bound from NW
    }
    else{
      // WALL reflects particles
      fnm[0] = f[idx(N,n,m) + 0*Nall]; // stationary 
      fnm[1] = f[idx(N,n,m) + 3*Nall]; // E bound from W
      fnm[2] = f[idx(N,n,m) + 4*Nall]; // N bound from S
      fnm[3] = f[idx(N,n,m) + 1*Nall]; // W bound from E
      fnm[4] = f[idx(N,n,m) + 2*Nall]; // S bound from N
      fnm[5] = f[idx(N,n,m) + 7*Nall]; // NE bound from SW
      fnm[6] = f[idx(N,n,m) + 8*Nall]; // NW bound from SE
      fnm[7] = f[idx(N,n,m) + 5*Nall]; // SW bound from NE
      fnm[8] = f[idx(N,n,m) + 6*Nall]; // SE bound from NW
    }
    
    // macroscopic density
    const dfloat rho = fnm[0]+fnm[1]+fnm[2]+fnm[3]+fnm[4]+fnm[5]+fnm[6]+fnm[7]+fnm[8];
    
    //    if(rho<1e-4){ printf("rho(%d,%d)=%g\n", n,m,rho); exit(-1); }
    
    // macroscopic momentum
    const dfloat delta2 = 1e-8;
    const dfloat Ux = (fnm[1] - fnm[3] + fnm[5] - fnm[6] - fnm[7] + fnm[8])*c/sqrt(rho*rho+delta2);
    const dfloat Uy = (fnm[2] - fnm[4] + fnm[5] + fnm[6] - fnm[7] - fnm[8])*c/sqrt(rho*rho+delta2);

    // compute equilibrium distribution
    dfloat feq[NSPECIES];
    lbmEquilibrium(c, rho, Ux, Uy, feq);

    // MRT stabilization
    const dfloat g0 = 1.f, g1 = -2.f, g2 = -2.f, g3 = -2.f, g4 = -2.f;
    const dfloat g5 = 4.f, g6 = 4.f, g7 = 4.f, g8 = 4.f;

    const dfloat R = g0*fnm[0] + g1*fnm[1] + g2*fnm[2]+ g3*fnm[3] + g4*fnm[4] + g5*fnm[5] + g6*fnm[6] + g7*fnm[7] + g8*fnm[8];
        
    // post collision densities
    fnm[0] -= tauinv*(fnm[0]-feq[0]) + (1.f-tauinv)*w0*g0*R*0.25f;
    fnm[1] -= tauinv*(fnm[1]-feq[1]) + (1.f-tauinv)*w1*g1*R*0.25f;
    fnm[2] -= tauinv*(fnm[2]-feq[2]) + (1.f-tauinv)*w2*g2*R*0.25f;
    fnm[3] -= tauinv*(fnm[3]-feq[3]) + (1.f-tauinv)*w3*g3*R*0.25f;
    fnm[4] -= tauinv*(fnm[4]-feq[4]) + (1.f-tauinv)*w4*g4*R*0.25f;
    fnm[5] -= tauinv*(fnm[5]-feq[5]) + (1.f-tauinv)*w5*g5*R*0.25f;
    fnm[6] -= tauinv*(fnm[6]-feq[6]) + (1.f-tauinv)*w6*g6*R*0.25f;
    fnm[7] -= tauinv*(fnm[7]-feq[7]) + (1.f-tauinv)*w7*g7*R*0.25f;
    fnm[8] -= tauinv*(fnm[8]-feq[8]) + (1.f-tauinv)*w8*g8*R*0.25f;
    
    // store new densities
    const int base = idx(N,n,m);
    fnew[base+0*Nall] = fnm[0];
    fnew[base+1*Nall] = fnm[1];
    fnew[base+2*Nall] = fnm[2];
    fnew[base+3*Nall] = fnm[3];
    fnew[base+4*Nall] = fnm[4];
    fnew[base+5*Nall] = fnm[5];
    fnew[base+6*Nall] = fnm[6];
    fnew[base+7*Nall] = fnm[7];
    fnew[base+8*Nall] = fnm[8];
  }
  
}

void lbmCheck(int N, int M, dfloat *f){

  int n,m,s;
  int nanCount = 0;
  for(s=0;s<NSPECIES;++s){
    for(m=0;m<=M+1;++m){
      for(n=0;n<=N+1;++n){
	nanCount += isnan(f[idx(N,n,m)+s*(N+2)*(M+2)]);
      }
    }
  }
  
  if(nanCount){   printf("found %d nans\n", nanCount); exit(-1); }
}



// set initial conditions (use uniform flow f everywhere)
void lbmInitialConditions(dfloat c, int N, int M, int *nodeType, dfloat *f){
  int n,m,s;
  dfloat feqIC[NSPECIES];
  dfloat feqWALL[NSPECIES];
  dfloat rhoIC = 1.;
  dfloat UxIC = 1.;
  dfloat UyIC = 0.;

  lbmEquilibrium(c, rhoIC, UxIC, UyIC, feqIC);
  lbmEquilibrium(c, rhoIC,    0.,  0., feqWALL);

  for(m=0;m<=M+1;++m){
    for(n=0;n<=N+1;++n){
      int base = idx(N, n, m);
      int s;

      if(n==0){
	for(s=0;s<NSPECIES;++s){
	  f[idx(N,n,m)+s*(N+2)*(M+2)] = feqIC[s];
	}
      }
      else{
	for(s=0;s<NSPECIES;++s){
	  f[idx(N,n,m)+s*(N+2)*(M+2)] = feqWALL[s];
	}
      }
    }
  }
}

int main(int argc, char **argv){

  if(argc!=3){
    printf("usage: ./lbm foo.png threshold\n");
    exit(-1);
  }

  // read threshold 
  dfloat threshold = atof(argv[2]);
  char *imageFileName = strdup(argv[1]);

  int N, M; // size of lattice
  unsigned char *rgb, *alpha;
  int *nodeType;
  lbmInput(imageFileName, threshold, &N, &M, &rgb, &alpha, &nodeType); 
  
  // physical parameters
  dfloat dx = .01;    // lattice node spacings 
  dfloat dt = dx*.1; // time step (also determines Mach number)
  dfloat c  = dx/dt; // speed of sound
  dfloat tau = .63; // relaxation rate
  dfloat Reynolds = 2./((tau-.5)*c*c*dt/3.);

  printf("Reynolds number %g\n", Reynolds);

  // create lattice storage
  dfloat *h_f    = (dfloat*) calloc((N+2)*(M+2)*NSPECIES, sizeof(dfloat));
  dfloat *h_fnew = (dfloat*) calloc((N+2)*(M+2)*NSPECIES, sizeof(dfloat));
  dfloat *h_tau  = (dfloat*) calloc((N+2)*(M+2), sizeof(dfloat));
  
  // set initial flow densities
  lbmInitialConditions(c, N, M, nodeType, h_f);
  lbmInitialConditions(c, N, M, nodeType, h_fnew);

  // set tau based on n index
  dfloat xo = .95;
  int n,m;
  for(m=0;m<=M+1;++m){
    for(n=0;n<=N+1;++n){
      dfloat x = ((double)n)/N;
      dfloat taunm = tau*(1 + 4*(1+tanh(20*(x-xo))));
      h_tau[idx(N,n,m)] = taunm;
    }
  }

  // DEVICE storage
  dfloat *c_f, *c_fnew, *c_tau;
  int *c_nodeType;
  
  cudaMalloc(&c_f,    (N+2)*(M+2)*NSPECIES*sizeof(dfloat));
  cudaMalloc(&c_fnew, (N+2)*(M+2)*NSPECIES*sizeof(dfloat));
  cudaMalloc(&c_nodeType, (N+2)*(M+2)*sizeof(int));
  cudaMalloc(&c_tau,      (N+2)*(M+2)*sizeof(dfloat));

  cudaMemcpy(c_f,    h_f,    (N+2)*(M+2)*NSPECIES*sizeof(dfloat), cudaMemcpyHostToDevice);
  cudaMemcpy(c_fnew, h_fnew, (N+2)*(M+2)*NSPECIES*sizeof(dfloat), cudaMemcpyHostToDevice);
  cudaMemcpy(c_nodeType, nodeType, (N+2)*(M+2)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(c_tau,  h_tau, (N+2)*(M+2)*sizeof(dfloat), cudaMemcpyHostToDevice);
  
  int Nsteps = 480000/2, tstep = 0, iostep = 100;

  // time step
  for(tstep=0;tstep<Nsteps;++tstep){

    // perform two updates
    dim3 T(TX,TY,1);
    dim3 B( (N+1+TX-1)/TX, (M+1+TY-1)/TY, 1);
    
    lbmUpdate <<< B, T >>> (N, M, c, c_tau, c_nodeType, c_f, c_fnew);
    lbmUpdate <<< B, T >>> (N, M, c, c_tau, c_nodeType, c_fnew, c_f);

    if(!(tstep%iostep)){ // output an image every iostep
      printf("tstep = %d\n", tstep);
      char fname[BUFSIZ];
      sprintf(fname, "bah%06d.png", tstep/iostep);

      cudaMemcpy(h_f, c_f, (N+2)*(M+2)*NSPECIES*sizeof(dfloat), cudaMemcpyDeviceToHost);
      lbmOutput(fname, nodeType, rgb, alpha, c, dx, N, M, h_f);

      //      sprintf(fname, "vel%06d.png", tstep/iostep);
      //      lbmOutputMicro(fname, nodeType, rgb, alpha, c, dx, N, M, h_f);

      lbmCheck(N,M,h_f);
    }
  }

  // output final result as image
  cudaMemcpy(h_f, c_f, (N+2)*(M+2)*NSPECIES*sizeof(dfloat), cudaMemcpyDeviceToHost);
  lbmOutput("bahFinal.png", nodeType, rgb, alpha, c, dx, N, M, h_f);

  exit(0);
  return 0;
}
  
