// http://www.caam.rice.edu/~timwar/CAAM210/Flows.html

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "occa.hpp"

extern "C"
{
#include "png_util.h"
}

#define dfloat double
#define FLUID 0
#define WALL 1
#define NSPECIES 9

// loop up 1D array index from 2D node coordinates
int idx(int N, int n, int m){
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

// weights used to compute equilibrium distribution (post collision)
const dfloat w0 = 4.f/9.f, w1 = 1.f/9.f, w2 = 1.f/9.f, w3 =  1.f/9.f;
const dfloat w4 = 1.f/9.f, w5 = 1.f/36.f, w6 = 1.f/36.f, w7 = 1.f/36.f, w8 = 1.f/36.f;

void lbmEquilibrium(const dfloat c,
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

  occa::device device;
  device.setup("mode=OpenMP, deviceID=0");
  
  // read threshold 
  dfloat threshold = atof(argv[2]);
  char *imageFileName = strdup(argv[1]);

  int N, M; // size of lattice
  unsigned char *rgb, *alpha;
  int *h_nodeType;
  lbmInput(imageFileName, threshold, &N, &M, &rgb, &alpha, &h_nodeType); 
  
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
  lbmInitialConditions(c, N, M, h_nodeType, h_f);
  lbmInitialConditions(c, N, M, h_nodeType, h_fnew);

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

  // OCCA DEVICE storage
  occa::memory o_f        = device.malloc((N+2)*(M+2)*NSPECIES*sizeof(dfloat),    h_f);
  occa::memory o_fnew     = device.malloc((N+2)*(M+2)*NSPECIES*sizeof(dfloat), h_fnew);
  occa::memory o_nodeType = device.malloc((N+2)*(M+2)*sizeof(int),         h_nodeType);
  occa::memory o_tau      = device.malloc((N+2)*(M+2)*sizeof(dfloat),           h_tau);

  // OCCA DEVICE kernel
  occa::kernel lbmUpdate = device.buildKernelFromSource("occaLBM.okl", "lbmUpdate");
  
  int Nsteps = 480000/2, tstep = 0, iostep = 100;

  // time step
  for(tstep=0;tstep<Nsteps;++tstep){

    // perform two updates
    lbmUpdate(N, M, c, o_tau, o_nodeType, o_f, o_fnew);
    lbmUpdate(N, M, c, o_tau, o_nodeType, o_fnew, o_f);

    if(!(tstep%iostep)){ // output an image every iostep
      printf("tstep = %d\n", tstep);
      char fname[BUFSIZ];
      sprintf(fname, "bah%06d.png", tstep/iostep);

      o_f.copyTo(h_f);
      lbmOutput(fname, h_nodeType, rgb, alpha, c, dx, N, M, h_f);

      lbmCheck(N,M,h_f);
    }
  }

  // output final result as image
  o_f.copyTo(h_f);
  lbmOutput("bahFinal.png", h_nodeType, rgb, alpha, c, dx, N, M, h_f);

  exit(0);
  return 0;
}
  
