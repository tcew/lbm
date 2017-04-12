// http://www.caam.rice.edu/~timwar/CAAM210/Flows.html

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "png_util.h"
#include <omp.h>

#define FLUID 0
#define WALL 1

#define NSPECIES 9

// loop up 1D array index from 2D node coordinates
int idx(int N, int n, int m){
  return n + m*(N+2);
}

void lbmInput(const char *imageFileName,
	      double threshold,
	      int *outN,
	      int *outM,
	      unsigned char **rgb,
	      unsigned char **alpha,
	      int **nodeType){

  int n,m, N,M;

  printf("*rgb = %p\n", rgb[0]);
  // read png file
  read_png(imageFileName, &N, &M, rgb, alpha);

  printf("after *rgb = %p\n", rgb[0]);

  // pad to guarantee space around obstacle
  int Npad = 2*N;
  int Mpad = 2*M;
  
  // threshold walls based on gray scale
  *nodeType = (int*) calloc((Npad+2)*(Mpad+2), sizeof(int));

  // mark pixels by gray scale intensity
  unsigned char *rgbPad   = (unsigned char*) calloc(3*(Npad+2)*(Mpad+2), sizeof(unsigned char));
  unsigned char *alphaPad = (unsigned char*) calloc((Npad+2)*(Mpad+2),   sizeof(unsigned char));
  int wallCount = 0;
  for(m=1;m<=M;++m){
    for(n=1;n<=N;++n){
      int offset = ((n-1)+(m-1)*N);
      double r = (*rgb)[3*offset+0];
      double g = (*rgb)[3*offset+1];
      double b = (*rgb)[3*offset+2];

      // center image in padded region (including halo zone)
      int id = idx(Npad,n+(N/4),m+(M/2));
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

  
  
  printf("wallCount = %d (%g percent of %d x %d nodes)\n", wallCount, 100.*((double)wallCount/((Npad+2)*(Mpad+2))), Npad, Mpad);
  
  *outN = Npad;
  *outM = Mpad;
}

  
void lbmOutput(const char *fname,
	       const int *nodeType,
	       unsigned char *rgb,
	       unsigned char *alpha,
	       int N,
	       int M,
	       const double *f){
  int n,m,s;
  FILE *bah = fopen(fname, "w");

  float plotMin = .95, plotMax = 1.05;
  for(m=0;m<=M+1;++m){
    for(n=0;n<=N+1;++n){
      int id = idx(N,n,m);

      // over write pixels in fluid region
      if(nodeType[id]==FLUID){
	unsigned char r,g,b,a;

	// reconstruct macroscopic density
	double rho = 0;
	for(s=0;s<NSPECIES;++s)
	  rho += f[id+s*(N+2)*(M+2)];
	rho = ((rho-plotMin)/(plotMax-plotMin)); // rescale
	r = 0;
	g = 255*(1.-rho);
	b = 255*rho;
	a = 255;

	rgb[idx(N,n,m)*3+0] = r;
	rgb[idx(N,n,m)*3+1] = g;
	rgb[idx(N,n,m)*3+2] = b;
	alpha[idx(N,n,m)] = a;
      }
    }
  }
  
  write_png(bah, N+2, M+2, rgb, alpha);

  fclose(bah);
}


void lbmEquilibrium(const double c,
		    const double rho,
		    const double Ux, 
		    const double Uy, 
		    double *  feq){

  // resolve macroscopic velocity into lattice particle velocity directions
  const double U2 = Ux*Ux+Uy*Uy;
  
  const double v0 = 0;
  const double v1 = +Ux/c;
  const double v2 = +Uy/c;
  const double v3 = -Ux/c;
  const double v4 = -Uy/c;
  const double v5 =  (+Ux+Uy)/c;
  const double v6 =  (-Ux+Uy)/c;
  const double v7 =  (-Ux-Uy)/c;
  const double v8 =  (+Ux-Uy)/c;
  
  // weights used to compute equilibrium distribution (post collision)
  const double w0 = 4.f/9.f, w1 = 1.f/9.f, w2 = 1.f/9.f, w3 =  1.f/9.f;
  const double w4 = 1.f/9.f, w5 = 1.f/36.f, w6 = 1.f/36.f, w7 = 1.f/36.f, w8 = 1.f/36.f;

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



// perform lattice streaming and collision steps
void lbmUpdate(const int N,                  // number of nodes in x
	       const int M,                  // number of nodes in y
	       const double c,                // speed of sound
	       const double tau,              // relaxation rate
	       const int   *  nodeType,      // (N+2) x (M+2) node types 
	       const double *  f,             // (N+2) x (M+2) x 9 fields before streaming and collisions
	       double *  fnew){               // (N+2) x (M+2) x 9 fields after streaming and collisions

  // loop counters
  int n,m;

  // number of nodes in whole array including halo
  int Nall = (N+2)*(M+2);
  
  // physics paramaters
  double tauinv = 1.f/tau;

  // loop over all non-halo nodes in lattice
#pragma omp parallel for private(n)
  for(m=1;m<M+1;++m){
    for(n=1;n<=N+1;++n){
      
      // discover type of node (WALL or FLUID)
      const int nt = nodeType[idx(N,n,m)];
      double fnm[NSPECIES];

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
      const double rho = fnm[0]+fnm[1]+fnm[2]+fnm[3]+fnm[4]+fnm[5]+fnm[6]+fnm[7]+fnm[8];

      if(rho<1e-4){ printf("rho(%d,%d)=%g\n", n,m,rho); exit(-1); }
      
      // macroscopic momentum
      const double delta2 = 1e-5;
      const double Ux = (fnm[1] - fnm[3] + fnm[5] - fnm[6] - fnm[7] + fnm[8])*c/sqrt(rho*rho+delta2);
      const double Uy = (fnm[2] - fnm[4] + fnm[5] + fnm[6] - fnm[7] - fnm[8])*c/sqrt(rho*rho+delta2);

      // compute equilibrium distribution
      double feq[NSPECIES];
      lbmEquilibrium(c, rho, Ux, Uy, feq);
      
      // post collision densities
      fnm[0] -= tauinv*(fnm[0]-feq[0]);
      fnm[1] -= tauinv*(fnm[1]-feq[1]);
      fnm[2] -= tauinv*(fnm[2]-feq[2]);
      fnm[3] -= tauinv*(fnm[3]-feq[3]);
      fnm[4] -= tauinv*(fnm[4]-feq[4]);
      fnm[5] -= tauinv*(fnm[5]-feq[5]);
      fnm[6] -= tauinv*(fnm[6]-feq[6]);
      fnm[7] -= tauinv*(fnm[7]-feq[7]);
      fnm[8] -= tauinv*(fnm[8]-feq[8]);

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
}

void lbmCheck(int N, int M, double *f){

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
void lbmInitialConditions(double c, int N, int M, int *nodeType, double *f){
  int n,m,s;
  double feqIC[NSPECIES];
  double feqWALL[NSPECIES];
  double rhoIC = 1.;
  double UxIC = 1.;
  double UyIC = 0.;

  lbmEquilibrium(c, rhoIC, UxIC, UyIC, feqIC);
  lbmEquilibrium(c, rhoIC,    0.,  0., feqWALL);

  for(m=0;m<=M+1;++m){
    for(n=0;n<=N+1;++n){
      int base = idx(N, n, m);
      int s;

      if(n==0){
	//      if(nodeType[idx(N,n,m)]==FLUID){
	for(s=0;s<NSPECIES;++s){
	  f[idx(N,n,m)+s*(N+2)*(M+2)] = feqIC[s];
	}
      }
      else{
	for(s=0;s<NSPECIES;++s){
	  f[idx(N,n,m)+s*(N+2)*(M+2)] = feqWALL[s];
	}
      }
#if 0
      printf("[%g,%g,%g,%g] ",
	     f[idx(N,n,m)+0*(N+2)*(M+2)],
	     f[idx(N,n,m)+1*(N+2)*(M+2)],
	     f[idx(N,n,m)+2*(N+2)*(M+2)],
	     f[idx(N,n,m)+3*(N+2)*(M+2)]);
#endif
    }
    //    printf("\n");
  }
}

int main(int argc, char **argv){

  if(argc!=3){
    printf("usage: ./lbm foo.png threshold\n");
    exit(-1);
  }

  int N, M; // size of lattice
  int n,m;

  // read threshold 
  double threshold = atof(argv[2]);
  char *imageFileName = strdup(argv[1]);

  unsigned char *rgb, *alpha;
  int *nodeType;
  lbmInput(imageFileName, threshold, &N, &M, &rgb, &alpha, &nodeType); 
  
  // physical parameters
  double dx = .01;    // lattice node spacings in x
  double dy = .01;
  double dt = dx*.1; // time step (also determines Mach number)
  double c  = dx/dt; // speed of sound
  double tau = .525; // relaxation rate
  double Reynolds = 2./((tau-.5)*c*c*dt/3.);

  printf("Reynolds number %g\n", Reynolds);

  // create lattice storage
  double *f    = (double*) calloc((N+2)*(M+2)*NSPECIES, sizeof(double));
  double *fnew = (double*) calloc((N+2)*(M+2)*NSPECIES, sizeof(double));

  // set initial flow densities
  lbmInitialConditions(c, N, M, nodeType, f);
  lbmInitialConditions(c, N, M, nodeType, fnew);
  
  int Nsteps = 30000/2, tstep = 0, iostep = 100;

  // time step
  for(tstep=0;tstep<Nsteps;++tstep){

    // perform two updates
    lbmUpdate(N, M, c, tau, nodeType, f, fnew);
    lbmUpdate(N, M, c, tau, nodeType, fnew, f);

    // check for nans
    lbmCheck(N, M, f);
    if(!(tstep%iostep)){
      printf("tstep = %d\n", tstep);
      char fname[BUFSIZ];
      sprintf(fname, "bah%06d.png", tstep);
      lbmOutput(fname, nodeType, rgb, alpha, N, M, f);
    }
  }


  // output result as image
  lbmOutput("bahFinal.png", nodeType, rgb, alpha, N, M, f);

  exit(0);
  return 0;
}
  
