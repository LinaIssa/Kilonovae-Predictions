#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "mylib.h"


/* --------------------------------------
 
 Escape probability from a uniform sphere
 
 Mattia Bulla, 2013-07-09
 
 ------------------------------------- */



double escape(double *ttot, int *nsc) {
    
    int Nb,N,Nsim,i,j;
    double z,R,a,b,k,t,prob=0;
    double *r,*theta,*d,*te;
    double n;
    
    
    /* NB: For each tau_tot we generate N packets and then for each packet other Nsim packets */
    N=*nsc;
    Nsim=1000;
    Nb=N*Nsim;
    
    R=1.0;      /* sphere radius */
    k=*ttot/R;      /* absorption coefficient */
    
    r=malloc(N*sizeof(int));
    theta=malloc(N*sizeof(int));
    d=malloc(N*sizeof(int));
    te=malloc(N*sizeof(int));
    
    srand((unsigned)time(NULL));
    
    
    /* Randomly place energy bundles in sphere and assign directions */
    
    for(i=0;i<N;i++){
        z=rand()/(RAND_MAX + 1.0);
        r[i]=R*pow(z,1./3.);        /* r=R*z^(1/3) see notes */
    }
    
    /* Randomly assign directions to bundles */
    
    for(i=0;i<N;i++){
        z=rand()/(RAND_MAX + 1.0);
        theta[i]=acos(1-2*z)*(180/PI);      /* direction of the bundles (see notes) */
    }
    
    /* Compute distance (d) to edge of sphere (using the cosine rule) and the relative tau (te) */
    
    for(i=0;i<N;i++){
        a=theta[i]*(PI/180);
        b=pow(R,2)-pow(r[i]*sin(a),2);
        d[i]= -r[i]*cos(a) + sqrt(b) ;  /* distance to edge of sphere (see notes for calculations) */
        te[i]=k*d[i];       /* compute tau edge (see notes) */
    }
    
    
    /* Randomly decide whether bundle is absorbed before it travels the distance d (use interaction pdf) */
    
    n=0;
    
    for(j=0;j<Nsim;j++){
        
        z=rand()/(RAND_MAX + 1.0);
        t=-log(z);      /* compute tau (see notes) */
        
        for(i=0;i<N;i++) if ( t>te[i] ) ++n;        /* n stores the number of escaped bundles */
        
    }
    
    prob=n/Nb;  /* compute escape probability */
    
    return prob;
    
}



/* --------------------------------------------------------------------------------------------------------------------
 
 Stokes Vector Transformation (see Chandrasekhar 1960 p. 38-42).
 However, it's not clear how the scattering angle is linked to the coordinates angles (t1,phi1,t2,phi2) --> stokesfun.c
 
 NB: the Stokes Vector is (Il,Ir,U) rather than (I,Q,U). Transformation: Il=(I+Q)/2 and Ir=(I-Q)/2.
 
 Mattia Bulla, 2013-09-23
 
 -------------------------------------------------------------------------------------------------------------------- */


void chandrafun(double *Il, double *Ir, double *U, double *t1, double *phi1, double *t2, double *phi2) {
    
    double u1,u2,dphi,a,b,c,d,e,f,n;
    double Il0,Il1,Il2;
    double Ir0,Ir1,Ir2;
    double U0,U1,U2;
    
    u1=cos(*t1);
    u2=cos(*t2);
    dphi=(*phi2)-(*phi1);
    a=1-pow(u1,2.);
    b=1-pow(u2,2.);
    c=cos(dphi);
    d=sin(dphi);
    e=cos(2*dphi);
    f=sin(2*dphi);
    
    
    Il0 = (*Il) * ( 2*a*b+pow(u1*u2,2.) )  +  (*Ir) * pow(u1,2.) ;
    Ir0 = (*Il) * pow(u2,2.)  +  (*Ir) ;
    U0 = 0 ;
    
    Il1 = sqrt(a) * sqrt(b) * ( (*Il) * 4 * u1 * u2 * c  +  (*U) * 2 * u1 * d  ) ;
    Ir1 = 0 ;
    U1 = sqrt(a) * sqrt(b) * (  - (*Il) * 2 * u2 * d  +  (*U) * c  ) ;
    
    Il2 = (*Il) * pow(u1,2.) * pow(u2,2.) * e   -   (*Ir) * pow(u1,2.) * e   +   (*U) * pow(u1,2.) * u2 * f  ;
    Ir2 = - (*Il) * pow(u2,2.) * e   +   (*Ir) * e   -   (*U) * u2 * f ;
    U2 = (*Il) * u1 * pow(u2,2.) * f   -   (*Ir) * u1 * f   +   (*U) * u1 * u2 * e ;
    
    
    *Il = 3./4. * ( Il0 + Il1 + Il2 ) ;
    *Ir = 3./4. * ( Ir0 + Ir1 + Ir2 ) ;
    *U = 3./2. * ( U0 + U1 + U2 ) ;
    
    /* Normalization (see Kasen) */
    
    n = sqrt ( pow(*Il,2.) + pow(*Ir,2.) + pow(*U,2.) ) ;
    
    *Il = *Il/n ;
    *Ir = *Ir/n ;
    *U = *U/n ;
    
    
}



/* --------------------------------------------------------------------------------------------------------------------
 
 Scattering (see Code & Whitney 1995 and Whitney 2011) with angular binning
 
 Mattia Bulla, 2013-09-25
 
 -------------------------------------------------------------------------------------------------------------------- */




/* ---------------------- Probability function P(M,i1) --------------------- */

double pdf(double *I, double *Q, double *U, double *phisc, double *M) {
    
    double a,If;
    
    a = pow(*M,2.) ;
    
    // NB: the rotational matrix R here is chosen in the clockwise direction ("+").
    // In Bulla+2015 equation (10) and (12) refer to the specific case shown in Fig.2 where the angle i2
    // is measured in the counter-clockwise direction. Therefore we use the clockwise rotation matrix but
    // with -i1. Here, instead, we calculate the angle in the clockwise direction from 0 to 2PI.
    // For instance, the i1 angle in Fig.2 of Bulla+2015 corresponds to 2PI-i1 here.
    // NB2: the i1 and i2 angles computed in the code (before and after scattering) are instead as in Bulla+2015
    If =  (a+1) * (*I)   +   (a-1) *  (   cos( 2 * (*phisc) ) * (*Q)   +  sin( 2 * (*phisc) ) * (*U)  );
    
    return If;
    
}


/* --------------- Maximum of the probability function P(M,i1) ------------- */

double maxpdf(double *I, double *Q, double *U) {
    
    int l,N;
    double phisc,M,Pmax,tmp;
    
    N=10000;
    Pmax=0; tmp=0; phisc=0; M=0;
    
    for(l=0;l<N;l++){
        
        phisc=2*PI* ( rand()/(RAND_MAX + 1.0) );
        M=1 - 2* ( rand()/(RAND_MAX + 1.0) );
        
        tmp = pdf(I,Q,U,&phisc,&M) ;
        
        if(tmp>Pmax) Pmax=tmp;

    }
    
    /*printf("\nPmax = %3.2f   find at (i1 = %1.3f) and (theta = %1.3f) \n\n",Pmax,imax*(180./PI),acos(mmax)*(180./PI)); */
    
    return Pmax;
    
}




/* --------------------------- Rejection Code ------------------------------ */

void rejection(double *I, double *Q, double *U, double *phisc, double *M) {
    
    double phisct,Mt,If,x;
    
    
    do {
        
        /* Sample i1 and M=cos(theta) from an isotropic distribution */
        
        phisct = 2*PI* gsl_rng_uniform(rng);
        Mt=2* gsl_rng_uniform(rng)-1;
        
        /* Calculate If=pdf(M,i1) (the Stoke parameter If in the reference frame) */
        
        If = pdf(I,Q,U,&phisct,&Mt) ;
        
        /* Compute x3*Pmax=2  */
        
        x = 2 * gsl_rng_uniform(rng) ;
        
    }
    
    while (x>If);  /* If x3*Pmax is less than pdf(M,i1) we accept M and i1. Otherwise, we resample. */
    
    
    *phisc = phisct ;
    *M = Mt ;
    
    
}






/* --------------------------- Stokes Parameter after rotation ------------------------------ */


void stokes_rotation_counterclock(double *I, double *Q, double *U, double *a, double *b) {
    
    double Qt,Ut;
    
    Qt = (*Q) * (*a) - (*U) * (*b);
    Ut = (*Q) * (*b) + (*U) * (*a);
    
    *Q = Qt;
    *U = Ut;
    
}


void stokes_rotation_clock(double *I, double *Q, double *U, double *a, double *b) {
    
    double Qt,Ut;
    
    Qt = (*Q) * (*a) + (*U) * (*b);
    Ut = -(*Q) * (*b) + (*U) * (*a);
    
    *Q = Qt;
    *U = Ut;

}



/* --------------------------- Stokes Parameter after scattering ------------------------------ */

void stokes_scattered_nonorm(double *I, double *Q, double *U, double *tsc) {
    
    double It,Qt,Ut,M,e;
    
    M = cos(*tsc) ;
    e = pow(M,2.) ;
    
    It = 0.75 * ( (*I) * (e+1) + (*Q) * (e-1) );
    Qt = 0.75 * ( (*I) * (e-1) + (*Q) * (e+1) ) ;
    Ut = 0.75 * ( 2 * M * (*U) ) ;
    
    *I = It ;
    *Q = Qt ;
    *U = Ut ;
    
}

void stokes_scattered_norm(double *I, double *Q, double *U, double *tsc) {
    
    double It,Qt,Ut,M,e;
    
    M = cos(*tsc) ;
    e = pow(M,2.) ;
    
    It = 0.75 * ( (*I) * (e+1) + (*Q) * (e-1) );
    Qt = 0.75 * ( (*I) * (e-1) + (*Q) * (e+1) ) ;
    Ut = 0.75 * ( 2 * M * (*U) ) ;
    
    *I = It / It;
    *Q = Qt / It;
    *U = Ut / It;
    
}


/* -------------------- Stokes Parameter Transformation (Hovenier 1983) ------------------------------ */

void stokes_norm(double *I, double *Q, double *U, double *tsc, double *a1, double *b1, double *a2, double *b2) {
    
    stokes_rotation_counterclock(I,Q,U,a1,b1);
    stokes_scattered_norm(I,Q,U,tsc);
    stokes_rotation_counterclock(I,Q,U,a2,b2);
    
}



/* --------------------------------- Rotation angle from the scattering plane ---------------------------------------------*/
/* -------- We need to rotate Stokes Parameters to (or from) the scattering plane from (or to) the meridian frame -------- */
/* ------------------------------- such that Q=1 is in the scattering plane and along ref1 -------------------------------- */


double rot_angle(double *n1, double *n2, double *ref1, double *ref2) {
    
    double ref1_sc[3],cos_stokes_rot_1, cos_stokes_rot_2, i=0;
    
    // ref1_sc is the ref1 axis in the scattering plane ref1 = n1 x ( n1 x n2 )
    ref1_sc[0] = n1[0] * dot(n1,n2) - n2[0];
    ref1_sc[1] = n1[1] * dot(n1,n2) - n2[1];
    ref1_sc[2] = n1[2] * dot(n1,n2) - n2[2];
    norm(ref1_sc);

    cos_stokes_rot_1 = dot(ref1_sc,ref1);
    cos_stokes_rot_2 = dot(ref1_sc,ref2);
    
    if (cos_stokes_rot_1<-1) cos_stokes_rot_1=-1;
    if (cos_stokes_rot_1>1) cos_stokes_rot_1=1;
    
    if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 > 0)) i = acos(cos_stokes_rot_1);
    if ((cos_stokes_rot_1 > 0) && (cos_stokes_rot_2 < 0)) i = 2 * acos(-1.) - acos(cos_stokes_rot_1);
    if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 < 0)) i = acos(-1.) + acos(fabs(cos_stokes_rot_1));
    if ((cos_stokes_rot_1 < 0) && (cos_stokes_rot_2 > 0)) i = acos(-1.) - acos(fabs(cos_stokes_rot_1));
    if (cos_stokes_rot_1 == 0) i = acos(-1.)/2.;
    if (cos_stokes_rot_2 == 0) i = 0.0 ;
    
    if (i!=i ) printf("Warning NaN: %3.6f \t %3.6f \t %3.6f \n",cos_stokes_rot_1,cos_stokes_rot_2,acos(cos_stokes_rot_1));
    
    
    return i;
    
}





/* --------------------------------------------------------------------------------------------------------------------
 
 Line Scattering. Include velocity gradient and use Sobolev approximation
 
 Mattia Bulla, 2013-11-19
 
 -------------------------------------------------------------------------------------------------------------------- */


double dot(double *x, double *y){
    
    double result;
    
    result = (x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]);
    
    return(result);
    
}

/**********************************************************************/
void cross_prod (double *v1, double *v2, double *v3) {
    
    v3[0] = (v1[1]*v2[2]) - (v2[1]*v1[2]);
    v3[1] = (v1[2]*v2[0]) - (v2[2]*v1[0]);
    v3[2] = (v1[0]*v2[1]) - (v2[0]*v1[1]);

}

void norm(double *x){
    
    double n;

    n = sqrt(dot(x,x));
    
    x[0]=x[0]/n;
    x[1]=x[1]/n;
    x[2]=x[2]/n;
    
}


void get_vel(double *x, double *t, double *vel) {
    
    vel[0] = x[0] / (*t) ;
    vel[1] = x[1] / (*t) ;
    vel[2] = x[2] / (*t) ;
    
    //vel[0] = 0  ;
    //vel[1] = 0 ;
    //vel[2] = 0 ;
    
}


void aberr(double *n1, double *v, double *n2) {
    
    double gamma_rel;
    double ndotv, fact2, fact1;
    double vsqr;
    
    vsqr = dot(v,v)/CLIGHT2;
    gamma_rel = 1./(sqrt(1 - vsqr));
    
    ndotv = dot(n1,v);
    fact1 = gamma_rel * (1 - (ndotv/CLIGHT));
    fact2 = (gamma_rel - (gamma_rel*gamma_rel*ndotv/(gamma_rel + 1)/CLIGHT))/CLIGHT;
    
    n2[0] = (n1[0] - (v[0] * fact2))/fact1;
    n2[1] = (n1[1] - (v[1] * fact2))/fact1;
    n2[2] = (n1[2] - (v[2] * fact2))/fact1;

    norm(n2);

    //double a;
    //a=dot(n1,n2);
    //if (n1[0]*n1[1]>0 && n2[0]*n2[1]<0) printf("%g %g %g %g \n",n1[0],n1[1],n2[0],n2[1]);
    
    //n2[0] = n1[0] ;
    //n2[1] = n1[1] ;
    //n2[2] = n1[2] ;

}


double doppler(double *n,double *vel) {
    
    double ndotv, fact1,gamma_rel,vsqr;//
    
    ndotv = dot(n,vel);
    
    vsqr = dot(vel,vel)/CLIGHT2;
    gamma_rel = 1./(sqrt(1 - vsqr));
    //gamma_rel = 1.;

    fact1 = gamma_rel*(1. - (ndotv/CLIGHT));
    
    //if (fabs(fact1-1) > 0.5)
    //{
    //    printf("Doppler factor > 1.05?? Abort.\n");
    //    printf("%g %g %g %g %g\n",sqrt(vsqr),vsqr,vel[0]*1.5*DAY,vel[1]*1.5*DAY,vel[2]*1.5*DAY);
    //    //return -1;
    //}
    
    return(fact1);
    
}


int check(double *x, double *ax, double *tcurrent, double *tgrid) {
    
    double result,ax_time[3];

    ax_time[0] = ax[0] * (*tcurrent) / (*tgrid) ;
    ax_time[1] = ax[1] * (*tcurrent) / (*tgrid) ;
    ax_time[2] = ax[2] * (*tcurrent) / (*tgrid) ;

    result = (x[0] * x[0] / ax_time[0] / ax_time[0]) + 
        (x[1] * x[1]  / ax_time[1] / ax_time[1] ) + (x[2] * x[2] / ax_time[2] / ax_time[2] ) ;

    if (result<1 && fabs(result-1)>1e-10) return 1;
    
    else return 0;
    
    
}

int check_grid(double *xt, int *nx, double *dgrid, double *tcurrent, double *tgrid) {

    double rmax, dgrid_t ;

    dgrid_t = (*dgrid) * (*tcurrent) / (*tgrid) ;

    rmax = 0.5 * dgrid_t * (*nx) ;
    
    if (fabs(xt[0])<=rmax && fabs(xt[1])<=rmax && fabs(xt[2])<=rmax) return 1;

    else return 0 ;

}



/* -------------------------------- New direction after scattering (Kalos & Whitlock 2008) ------------------------------ */

void newdir(double *tsc, double *phisc, double *n, double *n_out) {
    
    double a,b,c;
    
    
    if( fabs(n[2]) < 0.99999 ) {
        
        a = sin((*tsc))/sqrt(1.-pow(n[2],2.)) * ( n[1] * sin(*phisc) - n[0] * n[2] * cos(*phisc) ) + n[0] * cos(*tsc) ;
        b = sin(*tsc)/sqrt(1-pow(n[2],2.)) * ( - n[0] * sin(*phisc) - n[1] * n[2] * cos(*phisc) ) + n[1] * cos(*tsc) ;
        c = sin(*tsc) * cos(*phisc) * sqrt(1-pow(n[2],2.))  +  n[2] * cos(*tsc) ;
        
    }
    
    else {
        
        a = sin(*tsc) * cos(*phisc) ;
        b = sin(*tsc) * sin(*phisc) ;
        
        if( n[2]>0 ) c =  cos(*tsc) ;
        else c = - cos(*tsc) ;
        
    }
    
    n_out[0]=a;
    n_out[1]=b;
    n_out[2]=c;
    
}


double dist(double *x, double *ax, double *n) {
    
    
    int N=0;
    double xt,yt,zt,eps,d;
    
    if (ax[0] < ax[2]) eps=ax[0]/1e3;
    else eps=ax[2]/1e3;
    
    xt=x[0];   yt=x[1];   zt=x[2];
    
    while( pow(xt/ax[0],2.)+pow(yt/ax[1],2.)+pow(zt/ax[2],2.)<=1. ){
        
        xt += eps*n[0];
        yt += eps*n[1];
        zt += eps*n[2];
        
        N++;
        
    }
    
    d = eps * N;
    
    return d;
    
    
}


// Calculate the opacity to the boundary taking into account different shapes

double tau_boundary(double *x, int *photosph_flag, double *ax_source, int *nx, double *dgrid, double *ksc_free, double *ksc_rich, double *ksc_int, 
    double *kabs_par_free, double *kabs_par_rich, double *kabs_par_int, int *tdep, double *temission, double *tcurrent, double *tgrid, double *dens, 
    double *Ye, double *Ye_crit, double *n, double *nmax_eps, double *tau_max, double *dbound, double *nu_rf) {
    
    int inside_inner=0,inside_grid=0;
    double eps=0,xt[3];
    double dens_cell,Ye_cell;
    double tbound=0,kabstmp=0,ksctmp=0;
    //double extrad=0;
    double tfuture;
    double frac_eps,vel[3],nu_cmf;

    // Max path within a cell is l x sqrt(3) ~ l x 1.73 (diagonal)
    frac_eps = sqrt(3) / (*nmax_eps) ; 

    // Set initial quantities
    xt[0]=x[0];   xt[1]=x[1];   xt[2]=x[2];
    tfuture = (*tcurrent);

    // If vpkt has entered the inner region, kill it!
    inside_inner = check(xt,ax_source,&tfuture,temission); // Source is defined at temission 
    if (*photosph_flag==1 && inside_inner==1) return -1;

    // Check you are inside the grid
    inside_grid = check_grid(xt,nx,dgrid,&tfuture,tgrid); // Grid is defined at tgrid

    while( inside_grid == 1 && tbound < *tau_max) {

        eps = (*dgrid) * tfuture / (*tgrid) * frac_eps ;

        xt[0] += eps * n[0];
        xt[1] += eps * n[1];
        xt[2] += eps * n[2];

        if (*tdep==1) tfuture += eps / CLIGHT ;

        // If vpkt has entered the inner region, kill it!
        inside_inner = check(xt,ax_source,&tfuture,temission);
        if (*photosph_flag==1 && inside_inner==1) return -1;

        inside_grid = check_grid(xt,nx,dgrid,&tfuture,tgrid);

        if (inside_grid == 1) {

            // Calculate cmf frequency
            get_vel(xt,&tfuture,vel);
            nu_cmf = *nu_rf * doppler(n,vel) ;

            // Identify which cell you are in and read ejecta parameters
            read_cell(xt,&tfuture,tgrid,nx,dgrid,dens,Ye,&dens_cell,&Ye_cell) ;

            // Lanthanide-rich or lanthanide-free?
            if (Ye_cell<(*Ye_crit)) {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_rich,&tfuture) ;
                ksctmp = (*ksc_rich) * pow(*tgrid/tfuture,kabs_par_rich[4]) ; 
            }

            else if (Ye_cell>(*Ye_crit)) {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_free,&tfuture) ;
                ksctmp = (*ksc_free) * pow(*tgrid/tfuture,kabs_par_free[4]) ; 
            }

            else {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_int,&tfuture) ;
                ksctmp = (*ksc_int) * pow(*tgrid/tfuture,kabs_par_int[4])  ; 
            }

            // Increment opacity
            tbound += eps * (ksctmp + kabstmp) * dens_cell ;

            // Increment path
            *dbound += eps ;

        }
        
    }
    
//    // We need to add the last bit to the boundary (I rejected contributions when packets went out of boundary)
//
//    // Bring packet back inside
//    xt[0] = xt[0] - eps * n[0];
//    xt[1] = xt[1] - eps * n[1];
//    xt[2] = xt[2] - eps * n[2];
//    
//    if (*tdep==1) tfuture -= eps / CLIGHT ;
//    
//    // Identify which cell you are in and read ejecta parameters
//    read_cell(xt,&tfuture,tgrid,nx,dgrid,dens,Ye,&dens_cell,&Ye_cell) ;
//    
//    // Lanthanide-rich or lanthanide-free?
//    if (Ye_cell<(*Ye_crit)) {
//        kabstmp = abs_coefficient(nu,kabs_par_rich,tcurrent) ;
//        ksctmp = (*ksc_rich) * pow(*tgrid/tfuture,kabs_par_rich[4]) ; 
//    }
//    
//    else {
//        kabstmp = abs_coefficient(nu,kabs_par_free,tcurrent) ;
//        ksctmp = (*ksc_free) * pow(*tgrid/tfuture,kabs_par_free[4]) ; 
//    }
//
//    // Coefficient of outer ellipsoid/sphere defining the 
//    r[0] = 0.5 * (*dgrid) * tfuture / (*tgrid) * (*nx) ;
//    r[1] = 0.5 * (*dgrid) * tfuture / (*tgrid) * (*nx) ;
//    r[2] = 0.5 * (*dgrid) * tfuture / (*tgrid) * (*nx) ;
//    
//    extrad = dist(xt,r,n);
//    
//    tbound += extrad * (ksctmp + kabstmp) * dens_cell ;
//
//    *dbound = eps * N + extrad ;

    return tbound;
   
}



// Calculate the continuum opacity at the Sobolev point taking into account different shapes

double tau_cont(double *x, double *ksc_free, double *ksc_rich, double *ksc_int, double *kabs_par_free, double *kabs_par_rich, double *kabs_par_int, 
    int *tdep, double *tcurrent, double *tgrid, int *nx, double *dgrid, double *dens, 
    double *Ye, double *Ye_crit, double *n, double *s_l, double *nu_rf) {

    int inside_grid = 0;
    double xt[3],eps,s=0,ksctmp,kabstmp,tcont=0,dens_cell,Ye_cell,tfuture,nu_cmf,vel[3];

    xt[0]=x[0];   xt[1]=x[1];   xt[2]=x[2];

    eps= (*s_l) / 1e3 ;

    tfuture = (*tcurrent) ;

    while( s < (*s_l) ) {

        xt[0] += eps * n[0];
        xt[1] += eps * n[1];
        xt[2] += eps * n[2];

        if (*tdep==1) tfuture += eps / CLIGHT ;

        inside_grid = check_grid(xt,nx,dgrid,&tfuture,tgrid);
        
        // Calculate cmf frequency
        get_vel(xt,&tfuture,vel);
        nu_cmf = *nu_rf * doppler(n,vel) ;

        if (inside_grid == 1) {

            // Identify which cell you are in and read ejecta parameters
            read_cell(xt,&tfuture,tgrid,nx,dgrid,dens,Ye,&dens_cell,&Ye_cell) ;
            
            // Lanthanide-rich or lanthanide-free?
            if (Ye_cell<(*Ye_crit)) {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_rich,&tfuture) ;
                ksctmp = (*ksc_rich) * pow(*tgrid/tfuture,kabs_par_rich[4]) ; 
            }

            else if (Ye_cell>(*Ye_crit)) {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_free,&tfuture) ;
                ksctmp = (*ksc_free) * pow(*tgrid/tfuture,kabs_par_free[4]) ; 
            }

            else {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_int,&tfuture) ;
                ksctmp = (*ksc_int) * pow(*tgrid/tfuture,kabs_par_int[4])  ; 
            }

            // Increment opacity
            tcont += eps * (ksctmp + kabstmp) * dens_cell ;

        }

        // Packet outside both dynamical ejecta and wind
        // This means the the sobolev point was outside the outer boundary
        // I increment the opacity so that s will be bigger than d_boundary
        // This packet will be then rejected later on in the main function via the check function
        else {

            ksctmp = (*ksc_free) * pow(*tgrid/tfuture,kabs_par_free[4]) ;
            kabstmp = abs_coefficient(&nu_cmf,kabs_par_free,&tfuture);
            tcont += eps * (ksctmp + kabstmp) * dens_cell ;
        }

        // Increment path
        s += eps ;
    }

    return tcont;
    
}



// Calculate the continuum optical depth s 

double s_cont(double *x, double *ksc_free, double *ksc_rich, double *ksc_int, double *kabs_par_free, double *kabs_par_rich, double *kabs_par_int, 
    int *tdep, double *tcurrent, double *tgrid, double *n, int *nx, double *dgrid, 
    double *dens, double *Ye, double *Ye_crit, double *tau_ev, double *nmax_eps, double *tau_lines, double *nu_rf) {
    
    //printf("1 \n");    
    int inside_grid;
    double s=0,tautmp=0,xt[3],eps=0,kabstmp,ksctmp,dens_cell,Ye_cell,tfuture;
    double frac_eps,nu_cmf,vel[3];

    // Max path within a cell is l x sqrt(3) ~ l x 1.73 (diagonal)
    frac_eps = sqrt(3) / (*nmax_eps) ; 

    // Set initial quantities
    xt[0]=x[0];   xt[1]=x[1];   xt[2]=x[2];
    tfuture = (*tcurrent);

    while( tautmp < (*tau_ev) - (*tau_lines) ) {

        eps = (*dgrid) * tfuture / (*tgrid) * frac_eps ;

        xt[0] += eps * n[0];
        xt[1] += eps * n[1];
        xt[2] += eps * n[2];

        if (*tdep==1) tfuture += eps / CLIGHT ;

        inside_grid = check_grid(xt,nx,dgrid,&tfuture,tgrid);

        // Calculate cmf frequency
        get_vel(xt,&tfuture,vel);
        nu_cmf = *nu_rf * doppler(n,vel) ;

        if (inside_grid == 1) {
            
            // Identify which cell you are in and read ejecta parameters
            read_cell(xt,&tfuture,tgrid,nx,dgrid,dens,Ye,&dens_cell,&Ye_cell) ;
        
            // Lanthanide-rich or lanthanide-free?
            if (Ye_cell<(*Ye_crit)) {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_rich,&tfuture) ;
                ksctmp = (*ksc_rich) * pow(*tgrid/tfuture,kabs_par_rich[4]) ; 
            }

            else if (Ye_cell>(*Ye_crit)) {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_free,&tfuture) ;
                ksctmp = (*ksc_free) * pow(*tgrid/tfuture,kabs_par_free[4]) ; 
            }

            else {

                kabstmp = abs_coefficient(&nu_cmf,kabs_par_int,&tfuture) ;
                ksctmp = (*ksc_int) * pow(*tgrid/tfuture,kabs_par_int[4])  ; 
            }

            // Increment opacity
            //printf("%g %g %g %g %g %g \n",frac_eps,tfuture,eps,tautmp,(*tau_ev) - (*tau_lines),(ksctmp + kabstmp) * dens_cell);
            tautmp += eps * (ksctmp + kabstmp) * dens_cell ;
        }

        // Packet outside both dynamical ejecta and wind
        // This means the event was outside the outer boundary
        // This packet will be then rejected later on in the main function via the check function
        else tautmp += 1e45 ;

        //printf("%g %g %g %g\n",dens[ind],eps,tautmp,(*tau_ev) - (*tau_lines));

        if (tautmp < (*tau_ev) - (*tau_lines)) s += eps ;     
        
    }

    // Remove last contribution - that brought tautmp above (tau_ev - tau_lines)
    // And then add last path length to reach the event    

    if (tautmp < 1e44) {

        // Identify which cell you are in and read ejecta parameters
        read_cell(xt,&tfuture,tgrid,nx,dgrid,dens,Ye,&dens_cell,&Ye_cell) ;

        // Lanthanide-rich or lanthanide-free?
        if (Ye_cell<(*Ye_crit)) {
            kabstmp = abs_coefficient(&nu_cmf,kabs_par_rich,&tfuture) ;
            ksctmp = (*ksc_rich) * pow(*tgrid/tfuture,kabs_par_rich[4]) ; 
        }

        else if (Ye_cell>(*Ye_crit)) {
            kabstmp = abs_coefficient(&nu_cmf,kabs_par_free,&tfuture) ;
            ksctmp = (*ksc_free) * pow(*tgrid/tfuture,kabs_par_free[4]) ; 
        }

        else {
            kabstmp = abs_coefficient(&nu_cmf,kabs_par_int,&tfuture) ;
            ksctmp = (*ksc_int) * pow(*tgrid/tfuture,kabs_par_int[4])  ; 
        }

        tautmp = tautmp - eps * (ksctmp + kabstmp) * dens_cell ;

        s = s + ( (*tau_ev) - (*tau_lines) - tautmp ) / (ksctmp + kabstmp) / dens_cell ;
    }

    else s = 1e45;

    
    return s;

}



double rho_const(double *x, double *a, double *b, double *c, double *dir, double *ind, double *taumax){
    
    
    double r,xt,yt,zt,eps,tau,rho0t;
    double ratio;
    
    ratio = *c / (*a) ;
    tau=0;
    
    eps=*a/1e6;
    
    xt=x[0];   yt=x[1];   zt=x[2];
    
    while( sqrt( pow(xt/(*a),2.)+pow(yt/(*b),2.)+pow(zt/(*c),2.) )  <= 1. ){
        
        xt += eps*dir[0];
        yt += eps*dir[1];
        zt += eps*dir[2];
        
        r = sqrt( pow(xt,2.)+pow(yt,2.)+pow(zt/ratio,2.) ) ;
        
        tau +=  1 / pow( r , *ind ) * eps ;
        
    }
    
    
    rho0t = *taumax / tau ;
    
    return rho0t;
    
    
}




void move(double *x, double *Rmin, double *Rmax, double *dir, double *A, double *ind, double *tau_rand){
    
    
    double r,xt,yt,zt,eps_ins, eps,tau;
    double rho, cosbeta;
    
    tau = 0;
    
    eps_ins = *Rmin/1e2;
    eps = *Rmax/5e2;
    
    
    xt=x[0];   yt=x[1];   zt=x[2];
    
    r = sqrt( pow(xt,2.)+pow(yt,2.)+pow(zt,2.) ) ;
    
    // Inside Rmin: fly to Rmin
    if ( r < *Rmin ) {
        
        while(  sqrt( pow(xt,2.)+pow(yt,2.)+pow(zt,2.) ) <= *Rmin ) {

            xt += eps_ins * dir[0];
            yt += eps_ins * dir[1];
            zt += eps_ins * dir[2];

        }
        
    }
    
        

    
    while(  sqrt( pow(xt,2.)+pow(yt,2.)+pow(zt,2.) ) <= *Rmax  && tau < *tau_rand   ){
        
        r = sqrt( pow(xt+eps*dir[0]/2.,2.)+pow(yt+eps*dir[1]/2.,2.)+pow(zt+eps*dir[2]/2.,2.) ) ;
        cosbeta = (zt+eps*dir[2]/2.) / r ;
        
        rho = (*A) * pow ( (*Rmin/r) , *ind ) * (1+10*cosbeta*cosbeta) ;
        
        tau += rho * eps ;
        
        
        if ( tau < (*tau_rand) )  {  // go back to previous position
            
            xt += eps*dir[0];
            yt += eps*dir[1];
            zt += eps*dir[2];
        }
        
    }
    
    x[0]=xt;
    x[1]=yt;
    x[2]=zt;
        
        
    
    
    
}


double line_opac(double *tau0, double *x) {
    
    double tau;
    
    tau = *tau0;
    
    //Clumpy
    
    /*double r,theta,phi,,yc1,zc1,rclumpy1,yc2,zc2,rclumpy2;
    
    r = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ) ;
    theta = acos ( x[2] / r ) ;
    phi = atan2 ( x[1] , x[0] ) ;
    
    zc1 = 2.5e14 ;
    rclumpy1 = sqrt ( x[0]*x[0] + (x[1]-yc1)*(x[1]-yc1) + (x[2]-zc1)*(x[2]-zc1)   ) ;
    yc2 = 2.5e14 ;
    zc2 = 2e14 ;
    rclumpy2 = sqrt ( x[0]*x[0] + (x[1]-yc2)*(x[1]-yc2) + (x[2]-zc2)*(x[2]-zc2)   ) ;
    if ( rclumpy1 < 2.5e14 || rclumpy2 < 2e14)  tau = *tau0 ;
    */
     
    return tau;
    
}



double abs_coefficient(double *nu, double *kabs_par, double *tcurrent) {
    
    double wave, logkabs, kabs;
    double Bt, m, q ;

    wave = CLIGHT / (*nu) ;

    // Calculate B at time tcurrent
    Bt = kabs_par[3] * pow( ( *tcurrent / 1.5 / DAY ), kabs_par[4]) ;

    // Calculate m and q for intermediate region
    m = ( log10(Bt) - log10(kabs_par[1]) ) / ( kabs_par[2] - kabs_par[0] ) ;
    q = ( log10(kabs_par[1]) * kabs_par[2] - log10(Bt) * kabs_par[0] ) / ( kabs_par[2] - kabs_par[0] ) ;

    // Calculate logkabs for different wavelength ranges
    if (wave > kabs_par[2]) logkabs = log10( Bt ) ;
    else logkabs = m * wave + q ;

    kabs = pow(10.,logkabs);

    return kabs;

}




/* ----------------------- Routine to compute the meridian frame axes ref1 and ref2 ----------------------------------------*/

void meridian(double *n, double *ref1, double *ref2){
    
    
    // for ref_1 use (from triple product rule)
        
    ref1[0] = -1. * n[0] * n[2]/ sqrt(n[0]*n[0] + n[1]*n[1]);
    ref1[1] = -1. * n[1] * n[2]/ sqrt(n[0]*n[0] + n[1]*n[1]);
    ref1[2] = (1 - (n[2] * n[2]))/ sqrt(n[0]*n[0] + n[1]*n[1]);
    
    // for ref_2 use vector product of n_cmf with ref1
        
    ref2[0] = n[2] * ref1[1] - n[1] * ref1[2];
    ref2[1] = n[0] * ref1[2] - n[2] * ref1[0];
    ref2[2] = n[1] * ref1[0] - n[0] * ref1[1];
    
    
}


/* ----------------------- Routine to transform the Stokes Parameters from RF to CMF ----------------------------------------*/

void frame_transform(double *n_rf, double *Q, double *U, double *v, double *n_cmf) {
    
    double rot_angle,cos2rot_angle,sin2rot_angle,p,Q0,U0;
    double ref1[3],ref2[3],e_rf[3],e_cmf[3];
    double e_cmf_ref1, e_cmf_ref2, theta_rot=0;

    // Meridian frame in the RF
    meridian(n_rf,ref1,ref2);
    
    Q0 = *Q;
    U0 = *U;
    
    // Compute polarisation (which is invariant)
    p = sqrt(Q0*Q0+U0*U0);
    
    // We want to compute the angle between ref1 and the electric field
    rot_angle=0;
    
    if (p > 0) {
        
        cos2rot_angle = Q0/p;
        sin2rot_angle = U0/p;
        
        if ((cos2rot_angle > 0) && (sin2rot_angle > 0)) rot_angle = acos(Q0 / p) / 2. ;
        if ((cos2rot_angle < 0) && (sin2rot_angle > 0)) rot_angle = (acos(-1.) - acos(fabs(Q0 / p))) / 2. ;
        if ((cos2rot_angle < 0) && (sin2rot_angle < 0)) rot_angle = (acos(-1.) + acos(fabs(Q0 / p))) / 2. ;
        if ((cos2rot_angle > 0) && (sin2rot_angle < 0)) rot_angle = (2. * acos(-1.) - acos(fabs(Q0 / p))) / 2. ;
        if (cos2rot_angle == 0) {
            rot_angle = 0.25 * acos(-1);
            if (U0 < 0) rot_angle = 0.75 * acos(-1);
        }
        if (sin2rot_angle == 0) {
            rot_angle = 0.0;
            if (Q0 < 0) rot_angle = 0.5 * acos(-1);
                
        }
        
    }
    
    // Define electric field by linear combination of ref1 and ref2 (using the angle just computed)
    e_rf[0] =  cos(rot_angle) * ref1[0] - sin(rot_angle) * ref2[0];
    e_rf[1] =  cos(rot_angle) * ref1[1] - sin(rot_angle) * ref2[1];
    e_rf[2] =  cos(rot_angle) * ref1[2] - sin(rot_angle) * ref2[2];
    
    // Aberration
    aberr(n_rf,v,n_cmf);
    
    // Lorentz transformation of E
    lorentz(e_rf,n_rf,v,e_cmf);
    
    // Meridian frame in the CMF
    meridian(n_cmf,ref1,ref2);
    
    // Projection of E onto ref1 and ref2
    e_cmf_ref1 = e_cmf[0] * ref1[0] + e_cmf[1] * ref1[1] + e_cmf[2] * ref1[2];
    e_cmf_ref2 = e_cmf[0] * ref2[0] + e_cmf[1] * ref2[1] + e_cmf[2] * ref2[2];

    // Compute the angle between ref1 and the electric field
    if ((e_cmf_ref1 > 0) && (e_cmf_ref2 < 0)) theta_rot = acos(e_cmf_ref1);
    if ((e_cmf_ref1 < 0) && (e_cmf_ref2 < 0)) theta_rot = acos(-1.) - acos(fabs(e_cmf_ref1));
    if ((e_cmf_ref1 < 0) && (e_cmf_ref2 > 0)) theta_rot = acos(-1.) + acos(fabs(e_cmf_ref1));
    if ((e_cmf_ref1 > 0) && (e_cmf_ref2 > 0)) theta_rot = 2*acos(-1.) - acos(e_cmf_ref1);
    if (e_cmf_ref1 == 0) theta_rot = acos(-1.)/2. ;
    if (e_cmf_ref2 == 0) theta_rot = 0.0 ;
    if (e_cmf_ref1 > 1) theta_rot = 0.0 ;
    if (e_cmf_ref1 < -1) theta_rot = acos(-1.) ;
    
    // Compute Stokes Parameters in the CMF
    *Q = cos(2 * theta_rot ) * p ;
    *U = sin(2 * theta_rot ) * p ;

}


// Lorentz transformations from RF to CMF
void lorentz(double *e_rf, double *n_rf, double *v, double *e_cmf) {
    
    double beta[3],e_par[3], e_perp[3], b_rf[3], b_par[3], b_perp[3], vsqr, gamma_rel, v_cr_b[3], v_cr_e[3], b_cmf[3];
    
    beta[0] = v[0] / CLIGHT ;
    beta[1] = v[1] / CLIGHT ;
    beta[2] = v[2] / CLIGHT ;
    vsqr = dot(beta,beta);
    
    gamma_rel = 1./(sqrt(1 - vsqr));
    
    
    e_par[0] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[0] / (vsqr);
    e_par[1] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[1] / (vsqr);
    e_par[2] = (e_rf[0]*beta[0] + e_rf[1]*beta[1] + e_rf[2]*beta[2]) * beta[2] / (vsqr);
 
    e_perp[0] = e_rf[0] - e_par[0];
    e_perp[1] = e_rf[1] - e_par[1];
    e_perp[2] = e_rf[2] - e_par[2];
    
    b_rf[0]=n_rf[1]*e_rf[2] - n_rf[2]*e_rf[1];
    b_rf[1]=n_rf[2]*e_rf[0] - n_rf[0]*e_rf[2];
    b_rf[2]=n_rf[0]*e_rf[1] - n_rf[1]*e_rf[0];
    
    b_par[0] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[0] / (vsqr);
    b_par[1] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[1] / (vsqr);
    b_par[2] = (b_rf[0]*beta[0] + b_rf[1]*beta[1] + b_rf[2]*beta[2]) * beta[2] / (vsqr);
    
    b_perp[0] = b_rf[0] - b_par[0];
    b_perp[1] = b_rf[1] - b_par[1];
    b_perp[2] = b_rf[2] - b_par[2];
    
    
    
    v_cr_b[0]=beta[1]*b_rf[2] - beta[2]*b_rf[1];
    v_cr_b[1]=beta[2]*b_rf[0] - beta[0]*b_rf[2];
    v_cr_b[2]=beta[0]*b_rf[1] - beta[1]*b_rf[0];
    
    v_cr_e[0]=beta[1]*e_rf[2] - beta[2]*e_rf[1];
    v_cr_e[1]=beta[2]*e_rf[0] - beta[0]*e_rf[2];
    v_cr_e[2]=beta[0]*e_rf[1] - beta[1]*e_rf[0];
    
    
    e_cmf[0] = e_par[0] + gamma_rel * (e_perp[0] + v_cr_b[0]);
    e_cmf[1] = e_par[1] + gamma_rel * (e_perp[1] + v_cr_b[1]);
    e_cmf[2] = e_par[2] + gamma_rel * (e_perp[2] + v_cr_b[2]);
    
    b_cmf[0] = b_par[0] + gamma_rel * (b_perp[0] - v_cr_e[0]);
    b_cmf[1] = b_par[1] + gamma_rel * (b_perp[1] - v_cr_e[1]);
    b_cmf[2] = b_par[2] + gamma_rel * (b_perp[2] - v_cr_e[2]);
    
    norm(e_cmf);
    norm(b_cmf);
   
    
    
}



// Source emission
void source(int *source_shape, int *source_emission, double *ax_source, double *x, double *n_cmf) {

    double u,v,ri,mu,ti,phii,n_cmf_xyz[3],n_pos[3],n_sum[3],z_axis[3];
    double A,B,C,t,Rinner;
    int inside_source;
    double govergmax;
    double Nell[3], Nsph[3], Nperp[3];
    double a11,a12,a13,a21,a22,a23,a31,a32,a33;
    double x_test[3],eps;
    int ins_source;
    double t_unused;

    // -----------------------------------------------------
    // Sphere (surface) - Costant surface brightness
    // -----------------------------------------------------
    
    if (*source_shape==0 && *source_emission==0) {
        
        Rinner = ax_source[0] ;
        
        /* -------- Position ---------- */
        ri = Rinner ;
        mu = 2. *  gsl_rng_uniform(rng) -1.;
        ti = acos(mu);
        phii = 2. * PI * gsl_rng_uniform(rng);
        
        x[0] = ri * sin(ti) * cos(phii) ;
        x[1] = ri * sin(ti) * sin(phii) ;
        x[2] = ri * cos(ti) ;
        
        /* -------- Direction --------- */
        
        /* Constant brightness in the xyz system */
        u = sqrt ( gsl_rng_uniform(rng) ) ;
        v = 2. * PI * gsl_rng_uniform(rng);
        
        n_cmf_xyz[0] = sqrt(1-pow(u,2.)) * cos(v);
        n_cmf_xyz[1] = sqrt(1-pow(u,2.)) * sin(v);
        n_cmf_xyz[2] = u;
        
        /* Normal to the sphere surface at x */
        Nsph[0] = x[0] ;
        Nsph[1] = x[1] ;
        Nsph[2] = x[2] ;
        norm(Nsph);
        
        /* z_axis */
        z_axis[0] = 0 ;
        z_axis[1] = 0 ;
        z_axis[2] = 1 ;
        
        /* Angle between z and Nsph and orthoghonal axis */
        mu = dot(z_axis,Nsph) ;
        cross_prod(z_axis,Nsph,Nperp) ;
        norm(Nperp);
        
        /* Rotation matrix. See 5.2 http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ */
        a11 = pow(Nperp[0],2.) + (1 - pow(Nperp[0],2.)) * mu ;
        a12 = Nperp[0] * Nperp[1] * (1 - mu) - Nperp[2] * sqrt(1-pow(mu,2.)) ;
        a13 = Nperp[0] * Nperp[2] * (1 - mu) + Nperp[1] * sqrt(1-pow(mu,2.)) ;
        
        a21 = Nperp[0] * Nperp[1] * (1 - mu) + Nperp[2] * sqrt(1-pow(mu,2.)) ;
        a22 = pow(Nperp[1],2.) + (1 - pow(Nperp[1],2.)) * mu ;
        a23 = Nperp[1] * Nperp[2] * (1 - mu) - Nperp[0] * sqrt(1-pow(mu,2.)) ;
        
        a31 = Nperp[0] * Nperp[2] * (1 - mu) - Nperp[1] * sqrt(1-pow(mu,2.)) ;
        a32 = Nperp[1] * Nperp[2] * (1 - mu) + Nperp[0] * sqrt(1-pow(mu,2.)) ;
        a33 = pow(Nperp[2],2.) + (1 - pow(Nperp[2],2.)) * mu ;
        
        /* Rotate n_cmf_xyz */
        n_cmf[0] = a11 * n_cmf_xyz[0] + a12 * n_cmf_xyz[1] + a13 * n_cmf_xyz[2] ;
        n_cmf[1] = a21 * n_cmf_xyz[0] + a22 * n_cmf_xyz[1] + a23 * n_cmf_xyz[2] ;
        n_cmf[2] = a31 * n_cmf_xyz[0] + a32 * n_cmf_xyz[1] + a33 * n_cmf_xyz[2] ;
        
        
        // Check that the direction does not point inside the source
        eps=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])*0.1;
        
        x_test[0] = x[0] + n_cmf[0] * eps ;
        x_test[1] = x[1] + n_cmf[1] * eps ;
        x_test[2] = x[2] + n_cmf[2] * eps ;
        
        ins_source = check(x_test,ax_source,&t_unused,&t_unused);
        
        if (ins_source==1) {
            
            printf("ERROR: Constant surface brightness. Direction points inside the inner source \n");
            printf("ERROR: a_source %.3g b_source %.3g c_source %.3g \n",ax_source[0],ax_source[1],ax_source[2]);
            printf("ERROR: x[0] %.3g x[1] %.3g x[2] %.3g \n",x[0],x[1],x[2]);
            printf("ERROR: n[0] %.3g n[1] %.3g n[2] %.3g \n",n_cmf[0],n_cmf[1],n_cmf[2]);
            
            exit(0);
        }
        
    }
    
    
    // -----------------------------------------------------
    // Ellipsoid (surface) - Costant surface brightness
    // -----------------------------------------------------

    else if (*source_shape==1 && *source_emission==0) {
        
        /* -------- Position ---------- */
        /* We use the approach of Chen & Glotzer 2007 (see also Appendix for derivation of oblate) */
        do {
            
            mu = 2. *  gsl_rng_uniform(rng) -1.;
            phii = 2. * PI * gsl_rng_uniform(rng);
            
            x[0] = ax_source[0] * sqrt(1-pow(mu,2.)) * cos(phii) ;
            x[1] = ax_source[1] * sqrt(1-pow(mu,2.)) * sin(phii) ;
            x[2] = ax_source[2] * mu ;
            
            if ( ax_source[0] > ax_source[2] && ax_source[0] == ax_source[1])  
                govergmax = ax_source[2] * sqrt( ( pow(x[0],2.) + pow(x[1],2.) ) / pow(ax_source[0],4.) + pow(x[2],2.) / pow(ax_source[2],4.) ) ;
            else
                govergmax = ax_source[0] * sqrt( ( pow(x[0],2.) + pow(x[1],2.) ) / pow(ax_source[0],4.) + pow(x[2],2.) / pow(ax_source[2],4.) ) ;

            u = gsl_rng_uniform(rng);
            
        }
        
        while ( govergmax < u );
        
        
        /* -------- Direction --------- */
        
        /* Constant brightness in the xyz system */
        u = sqrt ( gsl_rng_uniform(rng) ) ;
        v = 2. * PI * gsl_rng_uniform(rng);
        
        n_cmf_xyz[0] = sqrt(1-pow(u,2.)) * cos(v);
        n_cmf_xyz[1] = sqrt(1-pow(u,2.)) * sin(v);
        n_cmf_xyz[2] = u;
        
        /* Normal to the ellipsoid surface at x */
        Nell[0] = 2 * x[0] / pow(ax_source[0],2.) ;
        Nell[1] = 2 * x[1] / pow(ax_source[1],2.) ;
        Nell[2] = 2 * x[2] / pow(ax_source[2],2.) ;
        norm(Nell);
        
        /* z_axis */
        z_axis[0] = 0 ;
        z_axis[1] = 0 ;
        z_axis[2] = 1 ;
        
        /* Angle between z and Nsph and orthoghonal axis */
        mu = dot(z_axis,Nell) ;
        cross_prod(z_axis,Nell,Nperp) ;
        norm(Nperp);
        
        /* Rotation matrix. See 5.2 http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ */
        a11 = pow(Nperp[0],2.) + (1 - pow(Nperp[0],2.)) * mu ;
        a12 = Nperp[0] * Nperp[1] * (1 - mu) - Nperp[2] * sqrt(1-pow(mu,2.)) ;
        a13 = Nperp[0] * Nperp[2] * (1 - mu) + Nperp[1] * sqrt(1-pow(mu,2.)) ;
        
        a21 = Nperp[0] * Nperp[1] * (1 - mu) + Nperp[2] * sqrt(1-pow(mu,2.)) ;
        a22 = pow(Nperp[1],2.) + (1 - pow(Nperp[1],2.)) * mu ;
        a23 = Nperp[1] * Nperp[2] * (1 - mu) - Nperp[0] * sqrt(1-pow(mu,2.)) ;
        
        a31 = Nperp[0] * Nperp[2] * (1 - mu) - Nperp[1] * sqrt(1-pow(mu,2.)) ;
        a32 = Nperp[1] * Nperp[2] * (1 - mu) + Nperp[0] * sqrt(1-pow(mu,2.)) ;
        a33 = pow(Nperp[2],2.) + (1 - pow(Nperp[2],2.)) * mu ;
        
        /* Rotate n_cmf_xyz */
        n_cmf[0] = a11 * n_cmf_xyz[0] + a12 * n_cmf_xyz[1] + a13 * n_cmf_xyz[2] ;
        n_cmf[1] = a21 * n_cmf_xyz[0] + a22 * n_cmf_xyz[1] + a23 * n_cmf_xyz[2] ;
        n_cmf[2] = a31 * n_cmf_xyz[0] + a32 * n_cmf_xyz[1] + a33 * n_cmf_xyz[2] ;
        
        
        // Check that the direction does not point inside the source
        eps=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])*0.1;
        
        x_test[0] = x[0] + n_cmf[0] * eps ;
        x_test[1] = x[1] + n_cmf[1] * eps ;
        x_test[2] = x[2] + n_cmf[2] * eps ;
        
        ins_source = check(x_test,ax_source,&t_unused,&t_unused);

        if (ins_source==1) {
         
            printf("ERROR: Constant surface brightness. Direction points inside the inner source \n");

            exit(0);
        }
        
        
    }
    

    
    // -----------------------------------------------------
    // Sphere (throughout the volume) - Isotropic emission
    // -----------------------------------------------------
    
    else if (*source_shape==0 && *source_emission==1) {
        
        Rinner = ax_source[0] ;
        
        /* Position */
        u = gsl_rng_uniform(rng);
        ri = Rinner * pow(u,1./3.);
        mu = 2. *  gsl_rng_uniform(rng) -1.;
        ti = acos(mu);
        phii = 2. * PI * gsl_rng_uniform(rng);
        
        x[0] = ri * sin(ti) * cos(phii) ;
        x[1] = ri * sin(ti) * sin(phii) ;
        x[2] = ri * cos(ti) ;
        
        /* Direction */
        
        u = 2. * gsl_rng_uniform(rng) - 1.;
        v = 2. * PI * gsl_rng_uniform(rng);
        
        n_cmf[0] = sqrt(1-pow(u,2.)) * cos(v);
        n_cmf[1] = sqrt(1-pow(u,2.)) * sin(v);
        n_cmf[2] = u;
        
        /* Move packet to the surface of the sphere (do not want it to interact inside the sphere) */
        
        /* 1st find unit vector to the position */
        norm(x) ;
        n_pos[0] = x[0] ;
        n_pos[1] = x[1] ;
        n_pos[2] = x[2] ;
        
        /* 2nd sum the position and the direction vectors */
        n_sum[0] = n_pos[0] + n_cmf[0] ;
        n_sum[1] = n_pos[1] + n_cmf[1] ;
        n_sum[2] = n_pos[2] + n_cmf[2] ;
        norm(n_sum) ;
        
        /* 3rd find intersection between n_sum and sphere */
        x[0] = Rinner * n_sum[0] ;
        x[1] = Rinner * n_sum[1] ;
        x[2] = Rinner * n_sum[2] ;
        
        
        // Check that the direction does not point inside the source
        eps=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])*0.1;
        
        x_test[0] = x[0] + n_cmf[0] * eps ;
        x_test[1] = x[1] + n_cmf[1] * eps ;
        x_test[2] = x[2] + n_cmf[2] * eps ;
        
        ins_source = check(x_test,ax_source,&t_unused,&t_unused);
        
        if (ins_source==1) {
            
            printf("ERROR: Isotropic emission. Direction points inside the inner source \n");
            
            exit(0);
        }
        
        
    }
    
    
    
    //--------------------------------------------
    // Ellipsoid (volume) - Isotropic emission
    //--------------------------------------------
    
    else if (*source_shape==1 && *source_emission==1) {
        
        /* Position */
        do {
            
            u = 2 * gsl_rng_uniform(rng) -1;
            x[0] = ax_source[0] * u ;
            u = 2 * gsl_rng_uniform(rng) -1;
            x[1] = ax_source[1] * u ;
            u = 2 * gsl_rng_uniform(rng) -1;
            x[2] = ax_source[2] * u ;
            
            inside_source = check(x,ax_source,&t_unused,&t_unused);
            
        }
        
        while ( inside_source==0 );
        
        
        /* Direction */
        
        u = 2. * gsl_rng_uniform(rng) - 1.;
        v = 2. * PI * gsl_rng_uniform(rng);
        
        n_cmf[0] = sqrt(1-pow(u,2.)) * cos(v);
        n_cmf[1] = sqrt(1-pow(u,2.)) * sin(v);
        n_cmf[2] = u;
        
        /* Move packet to the surface of the ellipsoid (do not want it to interact inside the ellipsoid) */
        
        /*  (i) write line in parametric form: x = x0 + t * nx, y = y0 + t * ny, z = z0 + t * nz
           (ii) put x, y and z into the ellipsoid equation and solve for t 
            NB  take the positive root, since the sign is already given by nx in the product t * nx */
        
        A = pow( n_cmf[0] / ax_source[0] ,2.) + pow( n_cmf[1] / ax_source[1] ,2.) + pow( n_cmf[2] / ax_source[2] ,2.)   ;
        B = 2 * ( x[0] * n_cmf[0] / pow(ax_source[0],2.) + x[1] * n_cmf[1] / pow(ax_source[1],2.) + x[2] * n_cmf[2] / pow(ax_source[2],2.) ) ;
        C = pow( x[0] / ax_source[0] ,2.) + pow( x[1] / ax_source[1] ,2.) + pow( x[2] / ax_source[2] ,2.) - 1 ;
        
        t = ( - B + sqrt( pow( - B ,2.) - 4 * A * C ) ) / (2 * A)  ;
        
        x[0] = x[0] + t * n_cmf[0] ;
        x[1] = x[1] + t * n_cmf[1] ;
        x[2] = x[2] + t * n_cmf[2] ;
        
        
        // Check that the direction does not point inside the source
        eps=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])*0.1;
        
        x_test[0] = x[0] + n_cmf[0] * eps ;
        x_test[1] = x[1] + n_cmf[1] * eps ;
        x_test[2] = x[2] + n_cmf[2] * eps ;
        
        ins_source = check(x_test,ax_source,&t_unused,&t_unused);
        
        if (ins_source==1) {
            
            printf("ERROR: Isotropic emission. Direction points inside the inner source \n");
            
            exit(0);
        }
        
        
    }
    
    
    
    else {
    
        printf("Source not recognized \n");
        
        exit(0);
       
    }
    
}


// Localize photosphere given a density map and ksc coefficient

double photosphere(int *nx, double *dgrid, double *tgrid, double *t, double *dens, double *ksc, double *alpha) {

    double xt[3], nt[3];
    double dgrid_t,rmax,dens_cell,vphot;
    double tau_phot,tautmp = 0;
    int indx,indy,indz=0,ind_cell;
    double ksc_t;

    // Optical depth at photosphere
    tau_phot = 2/3. ; 
    // Scattering coefficient at time t
    ksc_t = (*ksc) * pow(*tgrid/(*t),*alpha) ;

    // Calculate rmax
    dgrid_t = (*dgrid) * (*t) / (*tgrid) ;
    rmax = 0.5 * dgrid_t * (*nx) ;

    // Locate test particle at rmax along z
    xt[0] = 0 ;
    xt[1] = 0 ;
    xt[2] = rmax ;

    // Direction of test particle along -z
    nt[0] = 0 ;
    nt[1] = 0 ;
    nt[2] = -1 ;

    // Move test particle until you reach tau_phot
    while (tautmp < tau_phot && indz < (*nx)/2. ) {

        indx = (rmax - xt[0]) / dgrid_t ;
        indy = (rmax - xt[1]) / dgrid_t ;
        indz = (rmax - xt[2]) / dgrid_t ;

        ind_cell = pow((*nx),2) * indx + (*nx) * indy + indz ;

        dens_cell = dens[ind_cell] * pow((*tgrid)/(*t),3.) ;

        tautmp += dens_cell * ksc_t * dgrid_t;

        if (tautmp < tau_phot && indz < (*nx)/2. ) xt[2] = xt[2] + nt[2] * dgrid_t ;

    }

    // Bring packet back to last cell where optical depth was smaller than tau_phot
    xt[2] = xt[2] - nt[2] * dgrid_t ;

    vphot = sqrt(dot(xt,xt)) / (*t) ;

    return vphot;

}


// Source radioactive. Place packets according to the radioactive element distribution
void source_radio(int *nx, double *dgrid, double *xgrid, double *ygrid, double *zgrid, 
     double *tgrid, double *t, double *dens, double *msum, double *mtot, double *x, double *n_cmf) {

    // -------------------------- //
    // -------- Position -------- //
    // -------------------------- //

    int ngrid_tot;
    double zrand;
    int ind_below, ind_above, ind;
    
    ngrid_tot = pow(*nx,3.);
    
    // Place packet
    ind_above = ngrid_tot ;
    ind_below = 0;
    zrand = gsl_rng_uniform(rng);

    while (ind_above != (ind_below+1)) {

        ind = (ind_above + ind_below) / 2 ;

        // Accumulated mass is larger than the one drawn, cell is below
        if (msum[ind] > (zrand * (*mtot)) ) ind_above = ind;

        // Accumulated mass is smaller than the ose drawn, cell is above
        else ind_below = ind;

    }
    
    // Some sanity check
    if (msum[ind_below] > (zrand * (*mtot))) {

        printf("ind_below %d msum[ind_below] %g zrand*mtot %g\n", ind_below, msum[ind_below], zrand * (*mtot));
        exit(0);
    }

    if ((msum[ind_above] < (zrand * (*mtot))) && (ind_above != ngrid_tot)) {

        printf("ind_above %d msum[ind_above] %g zrand*mtot %g\n", ind_above, msum[ind_above], zrand * (*mtot));
        exit(0);
    }

    // Identify initial position
    ind = ind_below;

    x[0] = xgrid[ind];
    x[1] = ygrid[ind];
    x[2] = zgrid[ind];

    // Take into account that grid had expanded at time t
    x[0] = x[0] * (*t) / (*tgrid) ;
    x[1] = x[1] * (*t) / (*tgrid) ;
    x[2] = x[2] * (*t) / (*tgrid) ;

    // -------------------------- //
    // -------- Direction ------- //
    // -------------------------- //    

    double u,v;

    // Isotropic direction

    u = 2. * gsl_rng_uniform(rng) - 1.;
    v = 2. * PI * gsl_rng_uniform(rng);
    
    n_cmf[0] = sqrt(1-pow(u,2.)) * cos(v);
    n_cmf[1] = sqrt(1-pow(u,2.)) * sin(v);
    n_cmf[2] = u;

}



///****************************************************************************
double sample_planck(double *Temp, double tcurrent, double nu_min_r, double nu_max_r)
{
  double zrand,zrand2;
  double nu,nu_peak;
  double B_peak;
  int endloop;
  double T;

  double HOVERKB = 4.799243681748932e-11;
  double TWOHOVERCLIGHTSQUARED = 1.4745007e-47;
  
  // From Kasliwal+2017
  //T = T0 * pow( (tcurrent/DAY),alphaT ) ;
  T = (*Temp) ;
  nu_peak = 5.879e10 * T;
  B_peak = TWOHOVERCLIGHTSQUARED*pow(nu_peak,3) / (exp(HOVERKB*nu_peak/T) - 1);

  endloop = 0;
  while (endloop == 0)
  {
    zrand = gsl_rng_uniform(rng);
    zrand2 = gsl_rng_uniform(rng);
    nu = nu_min_r + zrand*(nu_max_r - nu_min_r);
    if (zrand2*B_peak <= TWOHOVERCLIGHTSQUARED*pow(nu,3) / (exp(HOVERKB*nu/T) - 1)) endloop = 1;
  }

  return nu;
}

///****************************************************************************

int search_index( double t,  double *time_data, int Nstep_data) {

    double min ;
    int index, i ;
    index = 0 ;
    min   = 1e30 ;
    for (i = 0; i<Nstep_data; i++) {
        if (fabs(t-time_data[i]) < min) {
            min   = fabs(t-time_data[i]) ;
            index = i ;
        }
    }
    return index ;
}

double temperature(double *tcurrent, double *tgrid, double *dgrid, int *nx, double *x, double *Ye, double *dens, double *eps_data, double *rate_data, double *time_data, int Nstep_data){
    
    double Temp, u_rad, Em=0, t ;
    double dens_cell, Ye_cell;
    int  i, index ;

    t = (*tcurrent) ;
    // Density of the cell 
    read_cell(x,tcurrent,tgrid,nx,dgrid,dens,Ye,&dens_cell,&Ye_cell);

    // Energy per mass 
    index = search_index(t,time_data,Nstep_data) ;
    for (i=0; i<index; i++){
        if (i==0) Em += eps_data[i] * rate_data[i] * time_data[i] * (time_data[i] / *tcurrent) ;
        else Em += eps_data[i] * rate_data[i] * (time_data[i] - time_data[i-1]) * (time_data[i] / *tcurrent) ;
    }

    // Temperature of the cell 
    u_rad = Em * dens_cell ;
    Temp = pow(CLIGHT*u_rad/(4*SIGMA),0.25);

    return Temp ;
}

double sample_emissivity(double T, double *tcurrent, double *tgrid, int *nx, double *dgrid, double *dens, double *Ye, double *Ye_crit, 
    double *x, double *ksc_free, double *ksc_rich, double *ksc_int, double *kabs_par_free, double *kabs_par_rich, double *kabs_par_int, double *emiss, double *nu_cmf) {

  double zrand1,zrand2,zrand3;
  int endloop;
  
  int inu,Nnu;
  double nu_min,nu_max;
  double nu,dnu,Bnu,snu_peak,snu;
  double dens_cell,Ye_cell,knu;
  double HOVERKB = 4.799243681748932e-11;
  double TWOHOVERCLIGHTSQUARED = 1.4745007e-47;

  // Select how to sample frequency
  zrand1 = gsl_rng_uniform(rng);

  // nu = nu_cmf (resonant)
  if ( zrand1 > (*emiss) ) return *nu_cmf ;

  // nu from Bnu x kappa_nu (emissivity)
  else {

      // From Kasliwal+2017
      //T = (*T0) * pow( (*tcurrent/DAY),*alphaT ) ;
      nu_min = 1e12 ;
      nu_max = 5e15 ;
      Nnu = 1e2 ;
      dnu = (nu_max - nu_min) / (double) Nnu ;

      snu_peak = 0;
      
      // Density and Ye for given cell
      read_cell(x,tcurrent,tgrid,nx,dgrid,dens,Ye,&dens_cell,&Ye_cell) ;

      // Find maximum emissivity
      for (inu=0; inu<Nnu; inu++) {

        // Frequency 
        nu = nu_min + inu * dnu ;

        // Planck funciton
        Bnu = TWOHOVERCLIGHTSQUARED*pow(nu,3) / (exp(HOVERKB*nu/T) - 1) ;

        // Opacity
        if (Ye_cell<(*Ye_crit)) knu = abs_coefficient(&nu,kabs_par_rich,tcurrent) ;
        else if (Ye_cell>(*Ye_crit)) knu = abs_coefficient(&nu,kabs_par_free,tcurrent) ;
        else knu = abs_coefficient(&nu,kabs_par_int,tcurrent) ;

        // Find peak
        if (Bnu * knu > snu_peak) snu_peak = Bnu * knu ;

        //printf("%g %g %g %g \n",nu,Bnu,knu,Bnu*knu);
      }

      //printf("%g %g \n",Ye_cell,snu_peak);


      // Do the frequency sampling
      endloop = 0;
      while (endloop == 0)
      {
        zrand2 = gsl_rng_uniform(rng);
        zrand3 = gsl_rng_uniform(rng);
        
        nu = nu_min + zrand2 * (nu_max - nu_min);
        
        // Snu        
        Bnu = TWOHOVERCLIGHTSQUARED*pow(nu,3) / (exp(HOVERKB * nu/T) - 1) ;
        if (Ye_cell<(*Ye_crit)) knu = abs_coefficient(&nu,kabs_par_rich,tcurrent) ;
        else if (Ye_cell>(*Ye_crit)) knu = abs_coefficient(&nu,kabs_par_free,tcurrent) ;
        else knu = abs_coefficient(&nu,kabs_par_int,tcurrent) ;       
        snu = Bnu * knu ;

        if (zrand3 * snu_peak <= snu) endloop = 1;
      }

      return nu;

  }

}



// Read into input grid and save density and Ye for the required cell
void read_cell(double *xt, double *t, double *tgrid, int *nx, double *dgrid, double *dens, double *Ye, double *dens_cell, double *Ye_cell) {

    int indx,indy,indz,ind_cell;
    double rmax, dgrid_t;

    // size of cells changes with time according to homologous expansion
    dgrid_t = (*dgrid) * (*t) / (*tgrid) ;

    rmax = 0.5 * dgrid_t * (*nx) ;

    indx = (rmax - xt[0]) / dgrid_t ;
    indy = (rmax - xt[1]) / dgrid_t ;
    indz = (rmax - xt[2]) / dgrid_t ;

    ind_cell = pow((*nx),2) * indx + (*nx) * indy + indz ;
    //printf("%g %g %g %g %g %d %d %d %d \n",*t,rmax,xt[0],xt[1],xt[2],indx,indy,indz, ind_cell);

    // density changes with time according to homologous expansion
    *dens_cell = dens[ind_cell] * pow((*tgrid)/(*t),3.)  ;
    
    *Ye_cell = Ye[ind_cell] ;

}


// Select which observer

void which_obs(double *tobs, double *phiobs, int *ind_obs, double *n_obs) {

    n_obs[0]=sin(tobs[*ind_obs])*cos(phiobs[*ind_obs]);
    n_obs[1]=sin(tobs[*ind_obs])*sin(phiobs[*ind_obs]);
    n_obs[2]=cos(tobs[*ind_obs]);
}


