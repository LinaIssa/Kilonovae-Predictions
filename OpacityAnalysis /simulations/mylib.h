
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>


#define PI 3.141592654
#define CLIGHT 2.99792458e10
#define CLIGHT2 8.9875518e20 /*Speed of light squared. */
#define H 6.6260755e-27 /* Planck constant */
#define MSUN 1.989e33 /* Solar mass */
#define DAY 86400
#define MPCTOCM 3.086e+24
#define SIGMA 5.6704e-5 /* Stefan-Boltzmann constant */
#define HCLIGHTOVERFOURPI 1.580764662876770e-17

const gsl_rng *rng;


double escape(double *ttot, int *nsc) ;
    
void chandrafun(double *Il, double *Ir, double *U, double *t1, double *phi1, double *t2, double *phi2) ;

double pdf(double *I, double *Q, double *U, double *phisc, double *M) ;

double maxpdf(double *I, double *Q, double *U) ;

void rejection(double *I, double *Q, double *U, double *phisc, double *M) ;

void stokes_rotation_counterclock(double *I, double *Q, double *U, double *a1, double *b1);

void stokes_rotation_counterclock2(double *I, double *Q, double *U, double *cosi) ;

void stokes_rotation_clock(double *I, double *Q, double *U, double *a1, double *b1);

void stokes_scattered_norm(double *I, double *Q, double *U, double *tsc);

void stokes_norm(double *I, double *Q, double *U, double *tsc, double *a1, double *b1, double *a2, double *b2);

void stokes_norm2(double *I, double *Q, double *U, double *tsc, double *phi1, double *phi2, double *cosi1, double *cosi2);

void stokes_scattered_nonorm(double *I, double *Q, double *U, double *tsc);

double rot_angle(double *n1, double *n2, double *ref1, double *ref2) ;

double rot_angle2(double *phi1, double *phi2, double *t1, double *t2, double *t_scatt) ;

//--------

double dot(double *x, double *y) ;

void cross_prod (double *v1, double *v2, double *v3) ;

void norm(double *x) ;

void get_vel(double *x, double *t, double *v) ;

void aberr(double *n1, double *v, double *n2);

double doppler(double *n,double *vel) ;

int check(double *x, double *ax, double *tcurrent, double *tgrid) ;

int check_grid(double *xt, int *ngrid_tot, double *dgrid, double *tcurrent, double *tgrid) ;

void newdir(double *tsc, double *phisc, double *n, double *n_out) ;

double dist(double *x, double *ax, double *n) ;

double tau_boundary(double *x, int *photosph_flag, double *ax_source, int *nx, double *dgrid, double *ksc_free, double *ksc_rich, double *ksc_int, double *kabs_par_free, double *kabs_par_rich, double *kabs_par_int, int *tdep, double *temission, double *tcurrent, double *tgrid, double *dens, double *Ye, double *Ye_crit, double *n, double *nmax_eps, double *tau_max, double *dbound, double *nu_rf) ;

double tau_cont(double *x, double *ksc_free, double *ksc_rich, double *ksc_int, double *kabs_par_free, double *kabs_par_rich, double *kabs_par_int, int *tdep, double *tcurrent, double *tgrid, int *ngrid_tot, double *dgrid, double *dens, double *Ye, double *Ye_crit, double *n, double *s_l, double *nu_rf) ;

double s_cont(double *x, double *ksc_free, double *ksc_rich, double *ksc_int, double *kabs_par_free, double *kabs_par_rich, double *kabs_par_int, int *tdep, double *tcurrent, double *tgrid, double *n, int *ngrid_tot, double *dgrid, double *dens, double *Ye, double *Ye_crit, double *tau_ev, double *nmax_eps, double *tau_lines, double *nu_rf);

double line_opac(double *tau0, double *x) ;

double abs_coefficient(double *nu, double *kabs_par, double *tcurrent) ;

double abs_coefficient_2(double *nu, double *kabs_par, double *tcurrent) ;

void meridian(double *n, double *ref1, double *ref2);

void frame_transform(double *n_rf, double *Q, double *U, double *v, double *n_cmf);

void lorentz(double *e_rf, double *n_rf, double *v, double *n_cmf) ;

double rho_const(double *x, double *a, double *b, double *c, double *n, double *ind, double *taumax);

void move(double *x, double *Rmin, double *Rmax, double *dir, double *A, double *ind, double *tau_rand);

void source(int *source_shape, int *source_emission, double *ax_source, double *x, double *n_cmf);

void source_radio(int *nx, double *dgrid, double *xgrid, double *ygrid, double *zgrid, double *tgrid, double *t, double *dens, double *msum, double *mtot, double *x, double *n_cmf) ;

double photosphere(int *nx, double *dgrid, double *tgrid, double *t, double *dens, double *ksc, double *alpha) ;

double sample_planck(double *Temp, double tcurrent, double nu_min_r, double nu_max_r);

int search_index(double t, double *time_data, int Nstep_data); 

double temperature(double *tcurrent, double *tgrid, double *dgrid, int *nx, double *x, double *Ye, double *dens, double *eps_data, double *rate_data, double *time_data, int Nstep_data);

double sample_emissivity(double Temp, double *tcurrent, double *tgrid, int *nx, double *dgrid, double *dens, double *Ye, double *Ye_crit, double *x, double *ksc_free, double *ksc_rich, double *ksc_int, double *kabs_par_free, double *kabs_par_rich, double *kabs_par_int, double *emiss, double *nu_cmf);

void read_cell(double *xt, double *t, double *tgrid, int *ngrid_tot, double *dgrid, double *dens, double *Ye, double *dens_cell, double *Ye_cell) ;

void which_obs(double *tobs, double *phiobs, int *ind_obs, double *n_obs) ;
