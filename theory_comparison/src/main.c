///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso, Anze Slosar                        //
//                                                                   //
//                                                                   //
// This file is part of CoLoRe.                                      //
//                                                                   //
// CoLoRe is free software: you can redistribute it and/or modify it //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CoLoRe is distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CoLoRe.  If not, see <http://www.gnu.org/licenses/>.   //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include <ccl.h>
#include "fftlog.h"

typedef struct {
  ccl_cosmology *cosmo;
  double hh;
  gsl_spline *spl_bz;
  double rsm;
  double zmax;
  double dz;
  char prefix_out[256];
  int do_ln;
} ParamLN;

int linecount(FILE *f)
{
  //////
  // Counts #lines from file
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

void *my_malloc(size_t size)
{
  void *ptrout=malloc(size);
  if(ptrout==NULL) {
    fprintf(stderr,"out of memory\n");
    exit(1);
  }
  
  return ptrout;
}

FILE *my_fopen(const char *path, const char *mode)
{
  FILE *fo=fopen(path,mode);
  if(fo==NULL) {
    fprintf(stderr,"Error opening file %s\n",path);
    exit(1);
  }

  return fo;
}

void param_ln_free(ParamLN *par)
{
  gsl_spline_free(par->spl_bz);
  ccl_cosmology_free(par->cosmo);
  free(par);
}    

ParamLN *param_ln_new(double om,double ob,double hh,double ns,double s8,char *fname_bz,
		      double rsm,double zmax,double dz,char *transfer,char *prefix_out,int do_ln)
{
  ParamLN *par=my_malloc(sizeof(ParamLN));
  ccl_configuration config=default_config;
  if(!strcmp(transfer,"eisenstein_hu")) {
    printf("Using EH\n");
    config.transfer_function_method=ccl_eisenstein_hu;
  }

  //Cosmology
  int status=0;
  par->cosmo=ccl_cosmology_create_with_lcdm_params(om-ob,ob,0,hh,s8,ns,config,&status);
  par->hh=hh;

  //Bias spline
  double *zarr,*bzarr;
  FILE *f=my_fopen(fname_bz,"r");
  int nz=linecount(f);
  rewind(f);
  zarr=my_malloc(nz*sizeof(double));
  bzarr=my_malloc(nz*sizeof(double));
  for(int i=0;i<nz;i++) {
    int stat=fscanf(f,"%lf %lf",&(zarr[i]),&(bzarr[i]));
    if(stat!=2) {
      fprintf(stderr,"Error reading file %s, line %d\n",fname_bz,i+1);
      exit(1);
    }
  }
  fclose(f);
  par->spl_bz=gsl_spline_alloc(gsl_interp_cspline,nz);
  gsl_spline_init(par->spl_bz,zarr,bzarr,nz);
  free(zarr);
  free(bzarr);
  
  //Smoothing
  par->rsm=rsm;

  //Z range
  par->zmax=zmax;
  par->dz=dz;

  //Output
  sprintf(par->prefix_out,"%s",prefix_out);  

  //Do LN
  par->do_ln=do_ln;
  printf("Do ln %d\n",par->do_ln);

  return par;
}

void write_predictions(ParamLN *par)
{
  printf("*** Writing predictions\n");
  // first generate k array, sufficiently finely spaced
  // note that we need to sufficiently pad on both ends

  const int Nk=10000;
  const double kmin=1e-3;
  const double kmax=50;;
  const double kminout=kmin;
  const double kmaxout=kmax;
  const double rminout=0.5;
  const double rmaxout=300.;
  double *ka=my_malloc(Nk*sizeof(double));
  double *pk=my_malloc(Nk*sizeof(double));
  double *pklin=my_malloc(Nk*sizeof(double));
  double *ra=my_malloc(Nk*sizeof(double));
  double *xi=my_malloc(Nk*sizeof(double));
  double *xilin=my_malloc(Nk*sizeof(double));
  for (int i=0; i<Nk; i++) ka[i]=kmin*pow((kmax/kmin),i*1.0/(Nk-1));
  FILE *fpk,*fxi;
  char fnamepk[256], fnamexi[256];
  double rsm2=par->rsm*par->rsm;
  gsl_interp_accel *intacc=gsl_interp_accel_alloc();

  // outter loop is over redshifts
  int status=0;
  for (double z=0; z<=par->zmax; z+=par->dz) {
    double a=1./(1+z);
    double bias=gsl_spline_eval(par->spl_bz,z,intacc);
    for (int i=0; i<Nk; i++) pklin[i]=ccl_linear_matter_power(par->cosmo,ka[i]*par->hh,a,&status)*pow(par->hh,3);
    pk2xi(Nk,ka,pklin,ra,xilin);
    // inner loop is over populations, ipop=-1 is the unbiased version
#ifdef _DEBUG
    printf("Writing predictions of redshift %lE: bias %lE\n",z,bias);
#endif //_DEBUG
    for (int i=0; i<Nk; i++) pk[i]=pklin[i]*bias*bias*exp(-rsm2*ka[i]*ka[i]);
    pk2xi(Nk,ka,pk,ra,xi);
    if(par->do_ln) {
      for (int i=0; i<Nk; i++) xi[i]=exp(xi[i])-1;
      xi2pk(Nk,ra,xi,ka,pk);
    }
    // now open the files
    sprintf(fnamepk,"%s_pk_z%.3lf.txt",par->prefix_out,z);
    sprintf(fnamexi,"%s_xi_z%.3lf.txt",par->prefix_out,z);
    fpk=fopen(fnamepk,"w");
    fprintf (fpk, "# k[h/Mpc] P_tt P_tl P_ll\n");
    fxi=fopen(fnamexi,"w");
    fprintf (fxi, "# r[Mpc/h] xi_tt xi_ll*b^2 xi_ll\n");
    for (int i=0; i<Nk; i++) {
      if ((ka[i]>=kminout) && (ka[i]<=kmaxout)) {
	fprintf (fpk,"%lE %lE %lE %lE\n",ka[i],pk[i], 
		 pklin[i]*bias*exp(-rsm2*ka[i]*ka[i]),
		 pklin[i]*exp(-rsm2*ka[i]*ka[i]));
      }
      if ((ra[i]>=rminout) && (ra[i]<=rmaxout))
	fprintf (fxi,"%lE %lE %lE %lE\n",ra[i],xi[i], xilin[i]*bias*bias, xilin[i]);
    }
    fclose(fpk);
    fclose(fxi);
  }

  free(ka);
  free(pk);
  free(pklin);
  free(ra);
  free(xi);
  free(xilin);
  gsl_interp_accel_free(intacc);
}

int main(int argc,char **argv)
{
  int do_ln;
  double om,ob,hh,ns,s8,rsm,zmax,dz;
  char fname_bz[256],prefix_out[256],transfer_function[256];
  if(argc!=13) {
    fprintf(stderr,"Usage : lnpred Om Ob hh ns s8 fname_bz r_smooth zmax delta_z transfer_function prefix_out do_LN\n");
    return 1;
  }
  om=atof(argv[1]);
  ob=atof(argv[2]);
  hh=atof(argv[3]);
  ns=atof(argv[4]);
  s8=atof(argv[5]);
  sprintf(fname_bz,"%s",argv[6]);
  rsm=atof(argv[7]);
  zmax=atof(argv[8]);
  dz=atof(argv[9]);
  sprintf(transfer_function,"%s",argv[10]);
  sprintf(prefix_out,"%s",argv[11]);
  do_ln=atoi(argv[12]);

  ParamLN *par=param_ln_new(om,ob,hh,ns,s8,fname_bz,rsm,zmax,dz,transfer_function,prefix_out,do_ln);
  write_predictions(par);
  param_ln_free(par);
  return 0;
}
