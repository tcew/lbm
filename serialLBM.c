// http://www.caam.rice.edu/~timwar/CAAM210/Flows.html

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "png_util.h"

#define FLUID 0
#define WALL 1

#define NSPECIES 9

// loop up 1D array index from 2D node coordinates
int idx(int N, int n, int m){
  return n + m*(N+2);
}


void lbmEquilibrium(const float c,
		    const float rho,
		    const float Ux, 
		    const float Uy, 
		    float *  feq){

  // resolve macroscopic velocity into lattice particle velocity directions
  const float U2 = Ux*Ux+Uy*Uy;
  
  const float v0 = 0;
  const float v1 = +Ux/c;
  const float v2 = +Uy/c;
  const float v3 = -Ux/c;
  const float v4 = -Uy/c;
  const float v5 =  (+Ux+Uy)/c;
  const float v6 =  (-Ux+Uy)/c;
  const float v7 =  (-Ux-Uy)/c;
  const float v8 =  (+Ux-Uy)/c;
  
  // weights used to compute equilibrium distribution (post collision)
  const float w0 = 4.f/9.f, w1 = 1.f/9.f, w2 = 1.f/9.f, w3 =  1.f/9.f;
  const float w4 = 1.f/9.f, w5 = 1.f/36.f, w6 = 1.f/36.f, w7 = 1.f/36.f, w8 = 1.f/36.f;

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
void update(const int N,                  // number of nodes in x
	    const int M,                  // number of nodes in y
	    const float c,                // speed of sound
	    const float tau,              // relaxation rate
	    const int   *  nodeType, // (N+2) x (M+2) node types 
	    const float *  f,     // (N+2) x (M+2) x 9 fields before streaming and collisions
	    float *  fnew){       // (N+2) x (M+2) x 9 fields after streaming and collisions

  // loop counters
  int n,m;

  // number of nodes in whole array including halo
  int Nall = (N+2)*(M+2);
  
  // physics paramaters
  float cinv = 1.f/c, tauinv = 1.f/tau;

  // loop over all non-halo nodes in lattice
  for(m=1;m<=M;++m){
    for(n=1;n<=N;++n){
      
      // discover type of node (WALL or FLUID)
      const int nt = nodeType[idx(N,n,m)];
      float fnm[NSPECIES];

      // translate particles on lattice
      if(nt != WALL){
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
      const float rho = fnm[0]+fnm[1]+fnm[2]+fnm[3]+fnm[4]+fnm[5]+fnm[6]+fnm[7]+fnm[8];
      
      // macroscopic momentum
      const float Ux = (fnm[1] - fnm[3] + fnm[5] - fnm[6] - fnm[7] + fnm[8])*c/rho;
      const float Uy = (fnm[2] - fnm[4] + fnm[5] + fnm[6] - fnm[7] - fnm[8])*c/rho;

      // compute equilibrium distribution
      float feq[NSPECIES];
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


int main(int argc, char **argv){

  if(argc!=3){
    printf("usage: ./lbm foo.png threshold\n");
    exit(-1);
  }

  int N, M; // size of lattice
  int n,m;

  // read threshold 
  float threshold = atof(argv[2]);
  
  // read png file
  unsigned char *rgb, *alpha;
  read_png(argv[1], &N, &M, &rgb, &alpha);
  
  // threshold walls based on gray scale
  int *ntype = (int*) calloc((N+2)*(M+2), sizeof(int));

  // mark pixels by gray scale intensity
  int wallCount = 0;
  for(m=1;m<=M;++m){
    for(n=1;n<=N;++n){
      int offset = 3*((n-1)+(m-1)*N);
      int r = rgb[offset+0];
      int g = rgb[offset+1];
      int b = rgb[offset+2];

      ntype[idx(N,n,m)] = WALL*(sqrt((double)(r*r+g*g+b*b))<threshold);
      wallCount += ntype[idx(N,n,m)];
      
      if(ntype[idx(N,n,m)]==WALL)
	rgb[offset+0] = 255;
      else
	rgb[offset+0] = 0;

      rgb[offset+1] = 0;
      rgb[offset+2] = 0;
    }
  }
  printf("wallCount = %d out of %d\n", wallCount, N*M);
 
  FILE *fp = fopen("foo.png", "w");
  write_png(fp, N, M, rgb, alpha);
  fclose(fp);

  // parameters
  float dx = 1;    // lattice node spacings in x
  float dy = 1;
  float dt = dx/3; // time step
  float c = dx/dt; // speed of sound
  float tau = .1; // relaxation rate

  // create lattice storage
  float *f    = (float*) calloc((N+2)*(M+2)*NSPECIES, sizeof(float));
  float *fnew = (float*) calloc((N+2)*(M+2)*NSPECIES, sizeof(float));
  
  // set initial conditions (use uniform flow f everywhere)
  float feqIC[NSPECIES];
  float feqWALL[NSPECIES];
  float rhoIC = 1.;
  float UxIC = 1.;
  float UyIC = 1.;

  lbmEquilibrium(c, rhoIC, UxIC, UyIC, feqIC);
  lbmEquilibrium(c, rhoIC, 0., 0., feqWALL);

  for(m=0;m<=M+1;++m){
    for(n=0;n<=N+1;++n){
      int base = idx(N, n, m);
      int s;
      if(ntype[idx(N,n,m)]==WALL){
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

  int Nsteps = 100/2, tstep = 0;

  for(tstep=0;tstep<Nsteps;++tstep){
    update(N, M, c, tau, ntype, f, fnew);
    update(N, M, c, tau, ntype, fnew, f);
  }

  FILE *bah = fopen("bah.png", "w");
  float *rho = (float*) calloc((N+2)*(M+2), sizeof(float));
  for(m=0;m<=M+1;++m){
    for(n=0;n<=N+1;++n){
      float tmp = 0;
      int s;
      for(s=0;s<NSPECIES;++s)
	tmp += f[idx(N,n,m)+s*(N+2)*(M+2)];

      rho[idx(N,n,m)] = tmp;
      if(isnan(tmp))
	printf("rho=%f\n", tmp);
    }
  }

  write_hot_png(bah, N+2, M+2, rho, 0, 1);
  fclose(bah);
}
  
