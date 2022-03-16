 /*
 
 Polarization toymodel for kilonovae
 
 Mattia Bulla, 2017-11-15
 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "mylib.h"


int main() {
    
    FILE *fsetup;
    char *filesetup;
    int unused __attribute__((unused));
    int Npar,Nobs,Nsim,isim,ipar;              
    double dum1,dum2,dum3;
    double Nph;
    int lf_contr;

    Nph = 1e4;                                  // # of photons packets
    
    lf_contr = 0; // If lf_contr = 1, print out flux contribution from lf region

    // Load parameters from simulation

    filesetup = "setup.txt";

    fsetup=fopen(filesetup,"r");                                        // open the setup file 
    unused = fscanf(fsetup,"%lf %lf %lf",&dum1,&dum2,&dum3);
    Npar = dum1;    // number of parameters
    Nobs = dum2;    // number of observers
    Nsim = dum3;    // number of simulations

    double sim_par[Nsim][Npar] ;   //2D arrays stocke les paramètres pour chaque simulation 

    for (isim=0;isim<Nsim;isim++) {
        for (ipar=0;ipar<Npar;ipar++) {
            
            unused = fscanf(fsetup,"%lf ",&dum1);
            sim_par[isim][ipar] = dum1;   
        }
    }

    // Start simulation
    for (isim=0;isim<Nsim;isim++) {

        // Check if simulation is present already       // à la fin flag_sim doit être égal à 1 pour que la simultaion commence  
        FILE *fq, *finfo, *fgrid, *fqdirect, *fqother;
        char filespec[128],fileinfo[128],filegrid[128], filedirec[128], fileother[128];
        char str[2];
        int flag_sim=1;

        if (lf_contr == 1) {
         
            snprintf(filespec,sizeof filespec,"outputs/nph%.1e_mej%.2f_phi%.0f_T%.1e_aopac%.1f_lfcontr_spec.txt",
                 Nph,sim_par[isim][0],sim_par[isim][1],sim_par[isim][2],sim_par[isim][3]) ;
            snprintf(fileinfo,sizeof fileinfo,"outputs/nph%.1e_mej%.2f_phi%.0f_T%.1e_aopac%.1f_lfcontr_info.txt",
                 Nph,sim_par[isim][0],sim_par[isim][1],sim_par[isim][2],sim_par[isim][3]) ;
            snprintf(filedirec,sizeof filedirec,"outputs/nph%.1e_mej%.2f_phi%.0f_T%.1e_aopac%.1f_direct_photons.txt",
                 Nph,sim_par[isim][0],sim_par[isim][1],sim_par[isim][2],sim_par[isim][3]) ;
            snprintf(fileother,sizeof fileother,"outputs/nph%.1e_mej%.2f_phi%.0f_T%.1e_aopac%.1f_other_photons.txt",
                 Nph,sim_par[isim][0],sim_par[isim][1],sim_par[isim][2],sim_par[isim][3]) ;
        }

        else {
         
            snprintf(filespec,sizeof filespec,"outputs/nph%.1e_mej%.2f_phi%.0f_T%.1e_aopac%.1f_spec.txt",
                 Nph,sim_par[isim][0],sim_par[isim][1],sim_par[isim][2],sim_par[isim][3]) ;
            snprintf(fileinfo,sizeof fileinfo,"outputs/nph%.1e_mej%.2f_phi%.0f_T%.1e_aopac%.1f_info.txt",
                 Nph,sim_par[isim][0],sim_par[isim][1],sim_par[isim][2],sim_par[isim][3]) ;
            snprintf(filedirec,sizeof filedirec,"outputs/nph%.1e_mej%.2f_phi%.0f_T%.1e_aopac%.1f_direct_photons.txt",
                 Nph,sim_par[isim][0],sim_par[isim][1],sim_par[isim][2],sim_par[isim][3]) ;
            snprintf(fileother,sizeof fileother,"outputs/nph%.1e_mej%.2f_phi%.0f_T%.1e_aopac%.1f_other_photons.txt",
                 Nph,sim_par[isim][0],sim_par[isim][1],sim_par[isim][2],sim_par[isim][3]) ;
        }

        finfo=fopen(fileinfo,"a+");
        fseek(finfo, 0, SEEK_END);
        if (ftell(finfo) != 0) {

            printf("\n%s already exists!!\n\n",fileinfo);
            printf("Should I overwrite this [Y/N]? ");
            if (fgets(str, sizeof str, stdin)) ;
            if (strcmp(str,"Y") != 0) flag_sim = 0;
            else {
                flag_sim = 1 ; 
                fq=fopen(filespec,"w");
                finfo=fopen(fileinfo,"w");
                fclose(fq);
                fclose(finfo);
            }
        }

        // Start        
        if (flag_sim == 1) {
        
            int l;
            double costobs,tobs[Nobs],phiobs[Nobs];
            double vphot,ax_photosph[3];
            double gamma,cos2gamma,sin2gamma;
            double ksc_free,ksc_rich,kabs,ksc,alb=0;
            double n_obs[3];
            double wavei,wavef,step_wave,tau_sob0,t;
            double eps_emiss,eps_emisstmp;
            int Nwave;
            int i,NLINE;
            double N_sc[Nobs],N_line[Nobs],N_abs[Nobs],N_direct[Nobs],N_inner[Nobs],N_real_inner,N_real_highop;
            double dgrid,dum1,dum2,dum3,dum4,dum5;
            int nx,ngrid_tot;
            int photosph_shape, photosph_emission,photosph_flag ;
            int inside;
            double T0,alphaT,mout,lbol,li,t0,sigma,alpha,eps_th,eps_0,eps_dot,factor,dMpc,Ye_crit,dens_cell,Ye_cell;
            double *msum, *xgrid, *ygrid, *zgrid, *dens, *Ye;
            double tcurrent, tgrid, mej_grid, tarrive_vpkt, ti, tf, step_t;
            int tdep,Ntime,itime;
            double kabs_par_free[5], kabs_par_rich[5];
            double tau_rpkt_max,tau_vpkt_max,nmax_eps;

            clock_t begin, end;
            double time_spent;
            begin = clock();

            snprintf(filegrid,sizeof filegrid,"models/grid_input_opang%.0f.txt",sim_par[isim][1]);          //input files
            //snprintf(filegrid,sizeof filegrid,"models/grid_input.txt");
            long seed;
            //seed = time(NULL);              // each timestep a new list of randomly generated nubers
            seed = 1;
            rng = gsl_rng_alloc(gsl_rng_ran3);
            gsl_rng_set (rng, seed);
            
            // ------------  Relevant parameters and flags ------------
            photosph_emission = 0 ;   // 0 constant surface brightness, 1 isotropic emission
            photosph_shape = 0 ;  // 0 sphere, 1 ellipsoid
            dMpc = 1;
            tdep = 1;   // Time dependence (1) or snapshot (0)
            photosph_flag = 0;   // Photons created at photosphere (1) or according to mass distr (0)
            Ye_crit = 0.25 ; // Critical Ye (lanth-free vs lanth-rich)
            tau_rpkt_max = 100; // Maximum optical depth above which we kill real packets
            tau_vpkt_max = 10; // Maximum optical depth above which we kill virtual packets
            nmax_eps = 10.; // Maximum number of steps within a grid cell (when moving photons)     

            // ------------ Temperature (From Kasliwal+ 2017) ------------
            T0 = sim_par[isim][2] ;
            alphaT = -0.4;

            // ----- Time quantities -----------

            if (tdep == 1) {                // if time dependancy

                // Time parameters for grid
                Ntime = 30 ;
                ti = 0.25 * DAY;
                tf = 15.25 * DAY;         
            }

            else {                          // in case of snapshots 

                Ntime = 1 ;
                ti = 1 * DAY ;
                tf = 1e30 * DAY ;
            }

            step_t = (tf - ti) / (double)Ntime;

            // ------------  Opacity coefficients ------------

            // Scattering opacity
            ksc_free = 0.01 ;
            ksc_rich = 0.01 ;

            // Analytic function of line opacity (see notes 25/01/19)

            // --- lanthanide-free in cgts
            kabs_par_free[0] = 0 ; // lambda1 in cm (constant with time)
            kabs_par_free[1] = 1e2 ; // kappa for lambda < lambda1 (constant with time)
            kabs_par_free[2] = 1e4 * 1e-8 ; // lambda2 in cm (constant with time)
            kabs_par_free[3] = 5e-3 ; // kappa for lambda > lambda2 (changing with time)
            if (tdep == 1) kabs_par_free[4] = sim_par[isim][3] ; // index of time-dependence opacity
            else kabs_par_free[4] = 0 ; // no time-dependence

            // --- lanthanide-rich
            kabs_par_rich[0] = 0 ; // lambda1 in cm (constant with time)
            kabs_par_rich[1] = 1e2 ; // kappa for lambda < lambda1 (constant with time)
            kabs_par_rich[2] = 1e4 * 1e-8 ; // lambda2 in cm (constant with time)
            kabs_par_rich[3] = 1 ; // kappa for lambda > lambda2 (changing with time)
            if (tdep == 1) kabs_par_rich[4] = sim_par[isim][3] ; // index of time-dependence opacity
            else kabs_par_rich[4] = 0 ; //o time-dependence

            // ------------  Viewing angle ------------

            for (i=0;i<Nobs;i++) {
                if (Nobs!=1) costobs = 1. / (Nobs-1) * i  ;
                else costobs = 1  ;                 
                tobs[i] = acos(costobs);
                phiobs[i] = 0*PI/180;  
            }


            // Rotation angle of the ellipsoid (with respect to z)
            gamma=0*PI/180; cos2gamma=cos(2*gamma); sin2gamma=sin(2*gamma);
            
            // ------------  Frequencies ------------

            // Resulting spectra
            wavei = 1000;                  // Initial wavelength
            wavef = 23000 ;                 // Final wavelength
            Nwave = 100 ;                     // Frequency binning.
            step_wave=(wavef-wavei)/(double)Nwave;

            // Emissivity coefficient (eps=0 resonance scatter, eps=1 complete redistribution, Kasen 2006 or Magee 2018)
            eps_emiss = 0.9 ;
        
            // ------------  Line scattering ------------

            NLINE = 1;
            double wave0[NLINE],nu0[NLINE];

            tau_sob0 = 0 ;              // Sobolev optical depth (parameter)
            
            wave0[0] = 6000 ;
            //wave0[1] = 5900 ;
            //wave0[2] = 6200 ;
            
            for (i=0;i<NLINE;i++) {
                
                nu0[i] = CLIGHT / wave0[i] * 1e8 ;      // clight en cm/s
                
                if (i>0 && wave0[i]<wave0[i-1]) {
                    
                    printf("Error: Lines are not ordered from blue to red \n");
                    return -1;
                }
            }

            // ------------ Set up grid ------------

            fgrid=fopen(filegrid,"r");
            unused=fscanf(fgrid,"%lf %lf %lf",&dum1,&dum2,&dum3);
            nx = dum1;
            ngrid_tot = pow(nx,3.);
            dgrid = dum2 * 2. / nx ;
            tgrid = dum3 ;

            xgrid = (double *)calloc(ngrid_tot,sizeof(double));
            ygrid = (double *)calloc(ngrid_tot,sizeof(double));
            zgrid = (double *)calloc(ngrid_tot,sizeof(double));
            dens = (double *)calloc(ngrid_tot,sizeof(double));
            Ye = (double *)calloc(ngrid_tot,sizeof(double));
            msum = (double *)calloc(ngrid_tot+1,sizeof(double));        // cumulative ejecta mass

            mej_grid = 0;
            for(i=0;i<ngrid_tot;i++){

                unused=fscanf(fgrid,"%lf %lf %lf %lf %lf",&dum1,&dum2,&dum3,&dum4,&dum5);
                xgrid[i] = dum1;
                ygrid[i] = dum2;
                zgrid[i] = dum3;
                dens[i] = dum4;
                Ye[i] = dum5;
                mej_grid += pow(dgrid,3.) * dens[i] ; // mass_grid
            }

            // Scale to required mass
            for(i=0;i<ngrid_tot;i++) dens[i] = dens[i] * (sim_par[isim][0] * MSUN/mej_grid) ;  // scaling_density to get m_ej required

            // ------------ Parameter for energy deposition (Korobkin+2012) ------------
            eps_0 = 2e18;
            eps_th = 0.5;
            t0 = 1.3;
            sigma = 0.11;
            alpha = 1.3;

            int inside_inner;
            double u=0,v=0,g=0;
            double tsc,phisc,M;
            double I,Q,U,I_virt,Q_virt,U_virt;
            double x[3],xl[3],vel[3],vel_rev[3],n_cmf[3],n_rf[3],n_out_cmf[3],n_obs_cmf[3],ref1[3],ref2[3];
            double nu_cmf,nu_rf;
            double s,tau_ev;
            double cos_scatt,t_scatt,tausob,tesc,pn,prob;
            double *Iwave_direct,*Qwave_direct,*Uwave_direct,*Iwave_sc,*Qwave_sc,*Uwave_sc,
                *Iwave_abs,*Qwave_abs,*Uwave_abs,*Iwave_line,*Qwave_line,*Uwave_line,*Iwave,*Qwave,*Uwave;
            double i1=0, cos2i1, sin2i1, i2=0, cos2i2, sin2i2;
            double color[NLINE],deltanu[NLINE],s_l[NLINE],tau_tot[NLINE],tau_sob[NLINE],tau_c[NLINE];
            double tau_sob_cum;
            int j,h,ind,flag_cont,flag_line,which_line;
            double dbound;
            double Narrived;
            int ind_obs,ind_wave,ind_t,ind_spec;
            double *Iwave_lf_direct,*Qwave_lf_direct,*Uwave_lf_direct,*Iwave_lf_sc,*Qwave_lf_sc,
                *Uwave_lf_sc,*Iwave_lf_abs,*Qwave_lf_abs,*Uwave_lf_abs,*Iwave_lf_line,*Qwave_lf_line,
                *Uwave_lf_line,*Iwave_lf,*Qwave_lf,*Uwave_lf;

            Iwave_direct=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Qwave_direct=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Uwave_direct=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Iwave_sc=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Qwave_sc=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Uwave_sc=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Iwave_abs=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Qwave_abs=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Uwave_abs=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Iwave_line=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Qwave_line=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Uwave_line=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            Iwave=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));    // total
            Qwave=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));    // total
            Uwave=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));    // total

            if (lf_contr == 1) {          
                
                Iwave_lf_direct=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Qwave_lf_direct=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Uwave_lf_direct=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Iwave_lf_sc=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Qwave_lf_sc=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Uwave_lf_sc=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Iwave_lf_abs=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Qwave_lf_abs=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Uwave_lf_abs=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Iwave_lf_line=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Qwave_lf_line=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Uwave_lf_line=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Iwave_lf=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Qwave_lf=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
                Uwave_lf=(double *)calloc(Nobs*Nwave*Ntime,sizeof(double));
            }


            // Virtual packet flags
            for (l=0;l<Nobs;l++) {

                N_direct[l] = 0 ;      // informations shown in the output file
                N_sc[l] = 0 ;
                N_line[l] = 0 ;
                N_abs[l] = 0 ;
                N_inner[l] = 0 ;
            }

            // Real packet flags            // informations shown in the output file
            N_real_inner = 0;
            N_real_highop = 0;

            FILE * fichier = NULL;
            fichier = fopen("outputs/map_photons_origin.txt", "w+");
            fprintf(fichier, "x \t \t y \t \t z \t \t tarrive \t time  \t tesc \t \t wavelength \t observer \n");

            for (itime=0;itime<Ntime;itime++) {     // first loop is time loop

                finfo=fopen(fileinfo,"a+");
                
                // ---- Time when packets are created -----
                t = ti + itime * step_t ;

                // ---- Location of photosphere -----
                vphot = photosphere(&nx,&dgrid,&tgrid,&t,dens,&ksc_free,&kabs_par_free[4]) ;

                if (vphot > 0.5 * dgrid * nx / tgrid) vphot = 0.5 * dgrid * nx / tgrid;

                ax_photosph[0] = vphot * t ;
                ax_photosph[1] = vphot * t ;
                ax_photosph[2] = vphot * t ;

                mout = 0;
                // ---- Mass above the photosphere -----
                for(i=0;i<ngrid_tot;i++){
                    msum[i] = mout ;
                    
                    if (photosph_flag==1) {
                        if ( sqrt( pow(xgrid[i],2.) + pow(ygrid[i],2.) + pow(zgrid[i],2.)) / tgrid > vphot ) 
                            mout += pow(dgrid,3.) * dens[i] ;
                    }

                    else mout += pow(dgrid,3.) * dens[i] ;
                }

                msum[ngrid_tot] = mout;

                // ---- Energy deposition -----
                eps_dot = eps_0 * pow( (0.5 - 1/PI * atan( (t-t0)/sigma) ), alpha) * (eps_th / 0.5);

                lbol = eps_dot * mout ;         // bolometric luminosity
                li = lbol / Nph ;               // each packet carried out the same amount of energy 

                // Scaling factor to transfor I stokes parameter to Flux density (1/4pi already included in prob)
                factor = li / step_wave / pow( (dMpc * MPCTOCM),2.) ;

                //printf("v > vph = %g c : M = %g Msun \n",vphot/CLIGHT,mout/MSUN);

                // Start propagation
                for(i=0;i<Nph;i++) {

                    tcurrent = t ;
                    
                    I=1;        // unpolarised photons
                    Q=0;
                    U=0;
                    
                    inside_inner = 0 ;                    
  
                    // Generate packets
                    if (photosph_flag == 1) source(&photosph_shape,&photosph_emission,ax_photosph,x,n_cmf);
                    else source_radio(&nx,&dgrid,xgrid,ygrid,zgrid,&tgrid,&tcurrent,dens,msum,&mout,x,n_cmf);       // location and the direction in the comoval frame

                    //printf("%g %g %g \n",x[0],x[1],x[2]);

                    // Initial frequencies (between nui and nuf) and energies in the CMF
                    //nu_cmf = sample_planck(T0,alphaT,tcurrent,CLIGHT/wavef*1e8,CLIGHT/wavei*1e8);
                    eps_emisstmp = 1;
                    nu_cmf = sample_emissivity(&T0,&alphaT,&tcurrent,&tgrid,&nx,&dgrid,dens,Ye,&Ye_crit,        //  gives the frequency by the rejection technique
                        x,&ksc_free,&ksc_rich,kabs_par_free,kabs_par_rich,&eps_emisstmp,&nu_cmf);

                    for (ind=0;ind<NLINE;ind++) {
                        
                        deltanu[ind] = nu_cmf - nu0[ind] ;
                        if (deltanu[ind] > 0) color[ind] = 1 ;
                        else color[ind] = -1;

                    }
                   
                    // ----- Compute values in the RF -------
                    
                    get_vel(x,&tcurrent,vel);       // getting the velocity
                    
                    if (lf_contr == 1) read_cell(x,&tcurrent,&tgrid,&nx,&dgrid,dens,Ye,&dens_cell,&Ye_cell);    //

                    // ------- Compute probability and total Stokes vector for virtual-packets
                    
                    for (l=0;l<Nobs;l++) {

                        which_obs(tobs,phiobs,&l,n_obs) ;
                        nu_rf = nu_cmf / doppler(n_obs,vel) ;       // rest frame frequency
                        tarrive_vpkt = tcurrent - (dot(x,n_obs)/CLIGHT);

                        if (fabs(tarrive_vpkt) < tf) {

                            N_abs[l] +=1;

                            for (ind=0;ind<NLINE;ind++) s_l[ind] = fabs( CLIGHT * deltanu[ind] * tcurrent / nu_rf ) ;  // Sobolev Point

                            tesc = tau_boundary(x,&photosph_flag,ax_photosph,&nx,&dgrid,&ksc_free,&ksc_rich,kabs_par_free,
                            kabs_par_rich,&tdep,&t,&tcurrent,&tgrid,dens,Ye,&Ye_crit,n_obs,&nmax_eps,&tau_vpkt_max,&dbound,&nu_rf);     // tau escape

                            fprintf(fichier, "%g \t %g \t %g \t %g \t %g \t %6g \t %3.3e \t %i \n", x[0], x[1], x[2], tarrive_vpkt/DAY, tcurrent/DAY, tesc, CLIGHT * 1e8 / nu_cmf, l);

                            if (tesc == -1) N_inner[l] += 1 ;
                        
                            else {
                            
                                N_direct[l] += 1 ;
                                pn=1./(4.*PI);      // isotropically generated
                            
                            // ------ Virtual packets can also come into resonance with the line
                            
                                tausob = 0 ;
                            
                            // Calculate tau_sob (Sum over NLINE lines)
                                for (ind=0;ind<NLINE;ind++) {
                                
                                    if(color[ind] == 1 && dbound > s_l[ind] ) {
                                    
                                    xl[0] = x[0] + s_l[ind] * n_obs[0] ;
                                    xl[1] = x[1] + s_l[ind] * n_obs[1] ;
                                    xl[2] = x[2] + s_l[ind] * n_obs[2] ;
                                    tausob += line_opac(&tau_sob0,xl) ;
                                    }
                                }

                        // If vpkt has entered the inner boundary, do not count its contribution
                            
                                tesc = tesc + tausob ;      // all the contributions of the source opacities
                            
                                prob=pn*exp(-tesc);

                                ind_obs = l ;
                                ind_wave = (CLIGHT * 1e8 / nu_rf - wavei) / step_wave ;     // integers
                                ind_t = (tarrive_vpkt - ti) / step_t ;
                                ind_spec = ind_wave + ind_t * Nwave + ind_obs * Nwave * Ntime ; // total index

                                if ( ind_wave >= 0 && ind_wave < Nwave && ind_t >= 0 && ind_t < Ntime ) {
                                
                                    Iwave_direct[ind_spec] += I * prob * factor;
                                    Qwave_direct[ind_spec] += Q * prob * factor;
                                    Uwave_direct[ind_spec] += U * prob * factor;

                                    if (lf_contr == 1 && Ye_cell > Ye_crit ) {

                                        Iwave_lf_direct[ind_spec] += I * prob * factor;
                                        Qwave_lf_direct[ind_spec] += Q * prob * factor;
                                        Uwave_lf_direct[ind_spec] += U * prob * factor;
                                    }
                                }
                            
                            }
                        }
                    
                    }

                    // -------------------------------------------------------------------
                    
                    
                    // RF frequency and Sobolev Point of the REAL packet
                    
                    vel_rev[0] = - vel[0] ;
                    vel_rev[1] = - vel[1] ;
                    vel_rev[2] = - vel[2] ;
                    aberr(n_cmf,vel_rev,n_rf);      // aberration relat
                    nu_rf = nu_cmf / doppler(n_rf,vel) ;
                    
                    for (ind=0;ind<NLINE;ind++) s_l[ind] = fabs( CLIGHT * deltanu[ind] * tcurrent / nu_rf ) ;            // Sobolev Point
                    
                    // Propagation

                    while ( tcurrent < tf + 1 * DAY ) {   // Kill real packet if tcurrent larger than tf + N*day (tarrive != tcurrent, see below)
                        
                        // -------------------- Which event?? (See Fig.1 in Mazzali and Lucy 1993) ---------------------        select the nature of the interaction: line or continuum opcacity
                        
                        tau_sob_cum = 0 ;
                        
                        for (ind=0;ind<NLINE;ind++) {

                            if (color[ind]==1) {
                                
                                xl[0] = x[0] + s_l[ind] * n_rf[0] ;
                                xl[1] = x[1] + s_l[ind] * n_rf[1] ;
                                xl[2] = x[2] + s_l[ind] * n_rf[2] ;
                                
                                // Line opacity
                            
                                tau_sob[ind] = line_opac(&tau_sob0,xl) ;
                                
                                tau_sob_cum += tau_sob[ind] ; // Cumulative line opacity (see below)
                                
                                // "Continuum" opacity at Sobolev point: tau_c
                                // tau_c is the opacity at the "ind" sobolev point excluding the line opacity of the "ind" line
                                // i.e. the continuum opacity up to that point + all the opacities given by the (ind-1) previous lines
                               
                                tau_c[ind] = tau_cont(x,&ksc_free,&ksc_rich,kabs_par_free,kabs_par_rich,
                                    &tdep,&tcurrent,&tgrid,&nx,&dgrid,dens,Ye,&Ye_crit,n_rf,&s_l[ind],&nu_rf);
            
                                if (ind!=0) tau_c[ind] = tau_c[ind] + tau_sob_cum - tau_sob[ind] ;   // need to count previous line opacities
                                
                                // Total opacity at Sobolev point: tau_tot
                                
                                tau_tot[ind] = tau_c[ind] + tau_sob_cum ;
                                
                            }
                            
                            else {
                                
                                // Setting all the values to 0 means that in the while loop below I will skip the "ind" line (interaction after the "ind" line)
                                tau_sob[ind] = 0 ;
                                tau_c[ind] = 0 ;
                                tau_tot[ind] = 0 ;
                            }
                            
            
                        }
                        

                        // Draw a random tau
                        tau_ev = - log( 1- gsl_rng_uniform(rng) );  // from the 1st approach
                        
                        // Decide which event
                        
                        flag_cont = 0 ;
                        flag_line = 0 ;
                        ind=0;
                        
                        while (ind < NLINE) {
                                            
                            // Continuum interaction before the line
                            if (tau_c[ind] > tau_ev) {

                                flag_cont = 1 ;
                                
                                if (ind==0) {
                                    
                                    tau_sob_cum = 0 ;
                                    s = s_cont(x,&ksc_free,&ksc_rich,kabs_par_free,kabs_par_rich,
                                    &tdep,&tcurrent,&tgrid,n_rf,&nx,&dgrid,dens,Ye,&Ye_crit,&tau_ev,&nmax_eps,&tau_sob_cum,&nu_rf);
                                }
                                
                                // Need to subtract PREVIOUS lines
                                else {

                                    tau_sob_cum=0;
                                    for (j=0;j<ind;j++) tau_sob_cum+=tau_sob[j];
                                    s = s_cont(x,&ksc_free,&ksc_rich,kabs_par_free,kabs_par_rich,           // s is the path
                                    &tdep,&tcurrent,&tgrid,n_rf,&nx,&dgrid,dens,Ye,&Ye_crit,&tau_ev,&nmax_eps,&tau_sob_cum,&nu_rf);
                                }
                                
                                break; // Exit while loop
                            }
                                
                            // Line scattering
                            else if (tau_tot[ind] > tau_ev) {
                                
                                flag_line = 1 ;
                                which_line = ind;
                                s = s_l[ind] ;
                                
                                break; // Exit while loop
                            }
            
                            // Continuum interaction after all the lines
                            // You can enter this only if ind == NLINE-1 and you have never interact ("break")
                            else if (ind==NLINE-1) {
                                
                                flag_cont = 1 ;
                                tau_sob_cum=0;
                                            
                                for (h=0;h<NLINE;h++) tau_sob_cum+=tau_sob[h];

                                s = s_cont(x,&ksc_free,&ksc_rich,kabs_par_free,kabs_par_rich,
                                    &tdep,&tcurrent,&tgrid,n_rf,&nx,&dgrid,dens,Ye,&Ye_crit,&tau_ev,&nmax_eps,&tau_sob_cum,&nu_rf);                        
                      
                                break; // Exit while loop
                            }

                            ind++ ;
                            
                        }   

                        // ----------------------------- Continuum event --------------------------------
                        
                        //printf("%g %g %g %g %g \n",x[0],x[1],x[2],s/t/CLIGHT,sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/t/CLIGHT);
                        
                        if ( flag_cont == 1 ) {     // always the case
                            
                            // ---- Move packet to the new RF position, check it's between inner and outer regions and update velocity ---      // first path
                       
                            x[0] = x[0] + s * n_rf[0] ;
                            x[1] = x[1] + s * n_rf[1] ;
                            x[2] = x[2] + s * n_rf[2] ;
                            
                            if (tdep==1) tcurrent += s / CLIGHT ;
                                    
                            inside_inner = check(x,ax_photosph,&tcurrent,&t); // Photosphere is created at time t
                            if (photosph_flag==1 && inside_inner == 1 ) {
                                
                                //printf("Inner \n");

                                N_real_inner += 1 ;
                                break;
                            }

                            inside = check_grid(x,&nx,&dgrid,&tcurrent,&tgrid);                    
                            if (inside == 0) break ;
                            
                            get_vel(x,&tcurrent,vel);
                            
                            // ---- Interaction ----
                            nu_cmf = nu_cmf - s * nu_rf / ( CLIGHT * tcurrent ) ;      // Update comoving frequency at interaction point see thesis 3.2.2 => comoving frequencies always redshifted
                            
                            // ----- Calculate absorption and scattering coefficients

                            // Identify which cell you are in and read ejecta parameters
                            read_cell(x,&tcurrent,&tgrid,&nx,&dgrid,dens,Ye,&dens_cell,&Ye_cell) ;

                            // Lanthanide-rich or lanthanide-free?
                            if (Ye_cell<Ye_crit) {

                                kabs = abs_coefficient(&nu_cmf,kabs_par_rich,&tcurrent) ;
                                ksc = ksc_rich * pow((tgrid/tcurrent),kabs_par_rich[4]) ; 
                            }

                            else {

                                kabs = abs_coefficient(&nu_cmf,kabs_par_free,&tcurrent) ;
                                ksc = ksc_free * pow((tgrid/tcurrent),kabs_par_free[4]) ; 
                            }

                            // Kill packet if optical depth in the cell is too high
                            if ( dens_cell * ( ksc + kabs ) * dgrid * tcurrent / tgrid > tau_rpkt_max ) {

                                //printf("%g \n",dens_cell * ( ksc + kabs ) * dgrid * tcurrent / tgrid);
                                N_real_highop +=1 ;
                                break ;
                            }

                            alb = ksc / ( ksc + kabs ) ;        // once we get to the interaction point

                            g=gsl_rng_uniform(rng);

                            // Absorption
                            if(g>alb) {

                                // Isotropically re-emitt packet with Q=0 and U=0 (like a grey opacity)
                                I_virt=1;  Q_virt=0;  U_virt=0;

                                // Sample new frequency using two-level approximation from Kasen et al 2006 (see also Magee et al. 2018)
                                nu_cmf = sample_emissivity(&T0,&alphaT,&tcurrent,&tgrid,&nx,&dgrid,dens,Ye,&Ye_crit,
                                    x,&ksc_free,&ksc_rich,kabs_par_free,kabs_par_rich,&eps_emiss,&nu_cmf);
                                            
                                // -------- Virtual packets -----
                                
                                for (l=0;l<Nobs;l++) {

                                    //break ;   // Destroy packet

                                    which_obs(tobs,phiobs,&l,n_obs) ;                             
                                    
                                    // Compute new frequency in the RF (packet direction is the observer one)
                                    nu_rf = nu_cmf / doppler(n_obs,vel);
                                
                                    // Remove packets arriving to late at the observer
                                    tarrive_vpkt = tcurrent - (dot(x,n_obs)/CLIGHT);
                                    
                                    if (fabs(tarrive_vpkt) < tf) {
                                        
                                        N_abs[l] += 1;

                                        tesc = tau_boundary(x,&photosph_flag,ax_photosph,&nx,&dgrid,&ksc_free,&ksc_rich,kabs_par_free,
                                            kabs_par_rich,&tdep,&t,&tcurrent,&tgrid,dens,Ye,&Ye_crit,n_obs,&nmax_eps,&tau_vpkt_max,&dbound,&nu_rf);

                                        if (tesc == -1) { // If vpkt has entered the inner boundary (tesc==-1), kill it!

                                            N_inner[l] += 1 ;
                                            break ;
                                        }
                                        
                                        pn=1./(4.*PI);
                                        
                                        // -- Virtual packets can also come into resonance with the line
                                        
                                        tausob = 0 ;
                                        
                                        // Calculate tau_sob (Sum over NLINE lines)
                                        for (ind=0;ind<NLINE;ind++) {
                                            
                                            if(color[ind] == 1 && dbound > (s_l[ind]-s) ) {
                                                
                                                xl[0] = x[0] + (s_l[ind]-s) * n_obs[0] ;
                                                xl[1] = x[1] + (s_l[ind]-s) * n_obs[1] ;
                                                xl[2] = x[2] + (s_l[ind]-s) * n_obs[2] ;
                                                tausob += line_opac(&tau_sob0,xl) ;
                                            }
                                        }
                                        
                                        tesc = tesc + tausob ;

                                        prob=pn*exp(-tesc);
                                        
                                        if (I_virt==I_virt && Q_virt==Q_virt && U_virt==U_virt) { // to avoid NaN

                                            ind_obs = l ;
                                            ind_wave = (CLIGHT * 1e8 / nu_rf - wavei) / step_wave ;
                                            ind_t = (tarrive_vpkt - ti) / step_t ;
                                            ind_spec = ind_wave + ind_t * Nwave + ind_obs * Nwave * Ntime ;
                                            
                                            if ( ind_wave >= 0 && ind_wave < Nwave && ind_t >= 0 && ind_t < Ntime ) {
                                                
                                                Iwave_abs[ind_spec] += I_virt * prob * factor;
                                                Qwave_abs[ind_spec] += Q_virt * prob * factor;
                                                Uwave_abs[ind_spec] += U_virt * prob * factor;

                                                if (lf_contr == 1 && Ye_cell > Ye_crit ) {

                                                        Iwave_lf_abs[ind_spec] += I * prob * factor;
                                                        Qwave_lf_abs[ind_spec] += Q * prob * factor;
                                                        Uwave_lf_abs[ind_spec] += U * prob * factor;
                                                }          
                                            }                                 
                                        }

                                    }

                                }

                                // ------- Real-packets -----
                                
                                // New Stokes Parameter (unpolarized!) and direction in the CMF (isotropic!)
                                I=1; Q=0; U=0;
                                u = 2 * gsl_rng_uniform(rng) - 1;
                                v = 2 * PI * gsl_rng_uniform(rng);
                                n_cmf[0] = sqrt(1-pow(u,2.)) * cos(v);
                                n_cmf[1] = sqrt(1-pow(u,2.)) * sin(v);
                                n_cmf[2] = u;
                                
                                // Compute the new packet direction in the RF
                                vel_rev[0] = - vel[0] ;
                                vel_rev[1] = - vel[1] ;
                                vel_rev[2] = - vel[2] ;
                                aberr(n_cmf,vel_rev,n_rf);
                                
                                // Compute new RF frequency
                                nu_rf = nu_cmf / doppler(n_rf,vel);
                                
                            }
                            
                            // Electron scattering
                            else {
                                
                                // ----------------------------------- Virtual packets -------------------------------------
                                
                                for (l=0;l<Nobs;l++) {

                                    I_virt=I;  Q_virt=Q;  U_virt=U;

                                    which_obs(tobs,phiobs,&l,n_obs) ;   // n_obs in rest frame

                                    // Compute new frequency in the RF (packet direction is the observer one)
                                    nu_rf = nu_cmf / doppler(n_obs,vel);

                                    // Remove packets arriving to late at the observer
                                    tarrive_vpkt = tcurrent - (dot(x,n_obs)/CLIGHT);    // delay nduced by light travel
                                    
                                    if (fabs(tarrive_vpkt) < tf) {
                                        
                                        N_sc[l] += 1 ;

                                        // Transform Stokes Parameters from the RF to the CMF // see M. Bulla 2015 + thesis 3.3
                                        
                                        frame_transform(n_rf,&Q_virt,&U_virt,vel,n_cmf);
                                        
                                        // Need to rotate Stokes Parameters in the scattering plane
                                        
                                        aberr(n_obs,vel,n_obs_cmf);
                                        
                                        meridian(n_cmf,ref1,ref2);
                                        
                                        /* This is the i1 angle of Bulla+2015, obtained by computing the angle between the
                                           reference axes ref1 and ref2 in the meridian frame and the corresponding axes
                                           ref1_sc and ref2_sc in the scattering plane. */
                                        i1 = rot_angle(n_cmf,n_obs_cmf,ref1,ref2);
                                        cos2i1 = cos(2 * i1) ;
                                        sin2i1 = sin(2 * i1) ;
                                        
                                        stokes_rotation_counterclock(&I_virt,&Q_virt,&U_virt,&cos2i1,&sin2i1);      // simple rotation
                                        
                                        // Scattering
                                        
                                        cos_scatt = dot(n_cmf,n_obs_cmf);
                                        t_scatt=acos(cos_scatt);
                                        
                                        pn=3./(16.*PI)*( 1+pow(cos_scatt,2.) + ( pow(cos_scatt,2.) - 1 ) * Q_virt ); // see matrix in paper 15 I= 1
                                        
                                        stokes_scattered_norm(&I_virt,&Q_virt,&U_virt,&t_scatt);
                                        
                                        // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame
                                        
                                        meridian(n_obs_cmf,ref1,ref2);
                                        
                                        /* This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
                                           reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
                                           meridian frame. NB: we need to add PI to transform THETA to i2 */
                                        i2 = PI+rot_angle(n_obs_cmf,n_cmf,ref1,ref2);
                                        cos2i2 = cos(2 * i2) ;
                                        sin2i2 = sin(2 * i2) ;
                                        
                                        stokes_rotation_clock(&I_virt,&Q_virt,&U_virt,&cos2i2,&sin2i2);

                                        // Transform Stokes Parameters from the CMF to the RF
                                        
                                        vel_rev[0] = - vel[0] ;
                                        vel_rev[1] = - vel[1] ;
                                        vel_rev[2] = - vel[2] ;
                                        
                                        frame_transform(n_obs_cmf,&Q_virt,&U_virt,vel_rev,n_obs);
                                        
                                        // Compute the opacity to the boundary
                                        
                                        tesc = tau_boundary(x,&photosph_flag,ax_photosph,&nx,&dgrid,&ksc_free,&ksc_rich,kabs_par_free,
                                            kabs_par_rich,&tdep,&t,&tcurrent,&tgrid,dens,Ye,&Ye_crit,n_obs,&nmax_eps,&tau_vpkt_max,&dbound,&nu_rf);

                                        if (tesc == -1) { // If vpkt has entered the inner boundary (tesc==-1), kill it!
                                            
                                            N_inner[l] += 1 ;
                                            break ;
                                        }
                                        
                                        // Virtual packets can also come into resonance with the line
                                        
                                        tausob = 0 ;
                                        
                                        // Calculate tau_sob (Sum over NLINE lines)
                                        for (ind=0;ind<NLINE;ind++) {
                                            
                                            if(color[ind] == 1 && dbound > (s_l[ind]-s) ) {
                                                
                                                xl[0] = x[0] + (s_l[ind]-s) * n_obs[0] ;
                                                xl[1] = x[1] + (s_l[ind]-s) * n_obs[1] ;
                                                xl[2] = x[2] + (s_l[ind]-s) * n_obs[2] ;
                                                tausob += line_opac(&tau_sob0,xl) ;
                                            }
                                        }
                                        
                                        tesc = tesc + tausob ;
                                        
                                        prob=pn*exp(-tesc);
                                        
                                        // Store values
                                        
                                        stokes_rotation_counterclock(&I_virt,&Q_virt,&U_virt,&cos2gamma,&sin2gamma);

                                        
                                        if (I_virt==I_virt && Q_virt==Q_virt && U_virt==U_virt) { // to avoid NaN

                                            ind_obs = l ;
                                            ind_wave = (CLIGHT * 1e8 / nu_rf - wavei) / step_wave ;
                                            ind_t = (tarrive_vpkt - ti) / step_t ;
                                            ind_spec = ind_wave + ind_t * Nwave + ind_obs * Nwave * Ntime ;

                                            if ( ind_wave >= 0 && ind_wave < Nwave && ind_t >= 0 && ind_t < Ntime ) {
                                                
                                                Iwave_sc[ind_spec] += I_virt * prob * factor;
                                                Qwave_sc[ind_spec] += Q_virt * prob * factor;
                                                Uwave_sc[ind_spec] += U_virt * prob * factor;

                                                if (lf_contr == 1 && Ye_cell > Ye_crit ) {

                                                        Iwave_lf_sc[ind_spec] += I * prob * factor;
                                                        Qwave_lf_sc[ind_spec] += Q * prob * factor;
                                                        Uwave_lf_sc[ind_spec] += U * prob * factor;
                                                }
                                            }
                                        }

                                    }
                                }
                                
                                // ----------------------------------- Real-packets -------------------------------------
                                
                                // Transform Stokes Parameters from the RF to the CMF
                                
                                frame_transform(n_rf,&Q,&U,vel,n_cmf);

                                // New packet direction in the CMF (rejection technique)
                                rejection(&I,&Q,&U,&phisc,&M);
                                tsc=acos(M);
                                newdir(&tsc,&phisc,n_cmf,n_out_cmf);
                                
                                // Need to rotate Stokes Parameters in the scattering plane
                                
                                meridian(n_cmf,ref1,ref2);
                                
                                /* This is the i1 angle of Bulla+2015, obtained by computing the angle between the
                                   reference axes ref1 and ref2 in the meridian frame and the corresponding axes
                                   ref1_sc and ref2_sc in the scattering plane. It is the supplementary angle of the
                                   scatt angle phisc chosen in the rejection technique */
                                i1 = rot_angle(n_cmf,n_out_cmf,ref1,ref2);
                                cos2i1 = cos(2 * i1) ;
                                sin2i1 = sin(2 * i1) ;
                                
                                stokes_rotation_counterclock(&I,&Q,&U,&cos2i1,&sin2i1);

                                // Scattering
                            
                                stokes_scattered_norm(&I,&Q,&U,&tsc);
                                
                                // Need to rotate Stokes Parameters out of the scattering plane to the meridian frame
                                
                                meridian(n_out_cmf,ref1,ref2);
                                
                                /* This is the i2 angle of Bulla+2015, obtained from the angle THETA between the
                                   reference axes ref1_sc and ref2_sc in the scattering plane and ref1 and ref2 in the
                                   meridian frame. NB: we need to add PI to transform THETA to i2 */
                                i2 = PI+rot_angle(n_out_cmf,n_cmf,ref1,ref2);
                                cos2i2 = cos(2 * i2) ;
                                sin2i2 = sin(2 * i2) ;
                                
                                stokes_rotation_clock(&I,&Q,&U,&cos2i2,&sin2i2);
                                
                              
                                // Transform Stokes Parameters from the CMF to the RF
                                
                                vel_rev[0] = - vel[0] ;
                                vel_rev[1] = - vel[1] ;
                                vel_rev[2] = - vel[2] ;
                                
                                frame_transform(n_out_cmf,&Q,&U,vel_rev,n_rf);
                                
                                // Compute new RF frequency
                                nu_rf = nu_cmf / doppler(n_rf,vel);
                                
                            
                            }
                            
                            // ---- Update packet status ----
                            
                            for (ind=0; ind<NLINE; ind++) {
                                
                                s_l[ind] = s_l[ind] - s ;   // Update sobolev point
                                
                                if (s_l[ind]<0) color[ind]=-1 ;   // Electron scattering was after the line. We will never interact with the line again
                            }
                            
                                                        
                        }
                        
                        
                        // ----------------------------------------- Line event --------------------------------------------
                        else if (flag_line == 1) {
                            
                            if (color[which_line]==-1) fprintf(finfo,"Warning! Double-counting the line-absorbed packets \n");
                           
                            // ---- Move packet to the new RF position, check it's between inner and outer regions and update velocity ---
                            
                            x[0] = x[0] + s * n_rf[0] ;
                            x[1] = x[1] + s * n_rf[1] ;
                            x[2] = x[2] + s * n_rf[2] ;
                            
                            inside_inner = check(x,ax_photosph,&tcurrent,&t); // Photosphere is created at time t
                            if (photosph_flag==1 && inside_inner == 1) {
                             
                                N_real_inner += 1 ;
                                break;
                            }

                            inside = check_grid(x,&nx,&dgrid,&tcurrent,&tgrid);
                            if (inside == 0 ) break ;                    

                            get_vel(x,&tcurrent,vel);
                            
                            if (lf_contr == 1) read_cell(x,&tcurrent,&tgrid,&nx,&dgrid,dens,Ye,&dens_cell,&Ye_cell);

                            nu_cmf = nu_cmf - fabs( s * nu_rf / ( CLIGHT * tcurrent ) ) ;      // Update comoving frequency at interaction point
                            
                            // -------- Virtual packets -----
                            
                            for (l=0;l<Nobs;l++) {

                                N_line[l] += 1 ;

                                I_virt=1;  Q_virt=0;  U_virt=0;

                                which_obs(tobs,phiobs,&l,n_obs) ;                               

                                // Remove packets arriving to late at the observer
                                tarrive_vpkt = tcurrent - (dot(x,n_obs)/CLIGHT);
                                
                                if (fabs(tarrive_vpkt) < tf) {

                                    pn=1./(4.*PI);
                                    
                                    tesc = tau_boundary(x,&photosph_flag,ax_photosph,&nx,&dgrid,&ksc_free,&ksc_rich,kabs_par_free,
                                        kabs_par_rich,&tdep,&t,&tcurrent,&tgrid,dens,Ye,&Ye_crit,n_obs,&nmax_eps,&tau_vpkt_max,&dbound,&nu_rf);
                                    
                                    if (tesc == -1) { // If vpkt has entered the inner boundary (tesc==-1), kill it!
                                        
                                        N_inner[l] += 1 ;
                                        break ;
                                    }
                                    
                                    // -- Virtual packets can also come into resonance with REDDER line
                                    
                                    tausob = 0 ;
                                    
                                    // Calculate tau_sob (Sum over REDDER lines)
                                    for (ind=which_line+1;ind<NLINE;ind++) {
                                        
                                        if(color[ind] == 1 && dbound > (s_l[ind]-s) ) {
                                            
                                            xl[0] = x[0] + (s_l[ind]-s) * n_obs[0] ;
                                            xl[1] = x[1] + (s_l[ind]-s) * n_obs[1] ;
                                            xl[2] = x[2] + (s_l[ind]-s) * n_obs[2] ;
                                            tausob += line_opac(&tau_sob0,xl) ;
                                        }
                                    }
                                    
                                    tesc = tesc + tausob ;
                                    
                                    prob=pn*exp(-tesc);
                                    
                                    if (I_virt==I_virt && Q_virt==Q_virt && U_virt==U_virt) { // to avoid NaN                                  
                                        
                                        // Compute new frequency in the RF (packet direction is the observer one) and store I,Q,U values in freqeuncy bins
                                        nu_rf = nu_cmf / doppler(n_obs,vel);

                                        ind_obs = l ;
                                        ind_wave = (CLIGHT * 1e8 / nu_rf - wavei) / step_wave ;
                                        ind_t = (tarrive_vpkt - ti) / step_t ;
                                        ind_spec = ind_wave + ind_t * Nwave + ind_obs * Nwave * Ntime ;

                                        if ( ind_wave >= 0 && ind_wave < Nwave && ind_t >= 0 && ind_t < Ntime ) {
                                            
                                            Iwave_line[ind_spec] += I_virt * prob * factor;
                                            Qwave_line[ind_spec] += Q_virt * prob * factor;
                                            Uwave_line[ind_spec] += U_virt * prob * factor;

                                            if (lf_contr == 1 && Ye_cell > Ye_crit ) {

                                                    Iwave_lf_line[ind_spec] += I * prob * factor;
                                                    Qwave_lf_line[ind_spec] += Q * prob * factor;
                                                    Uwave_lf_line[ind_spec] += U * prob * factor;
                                            }                                           
                                        }
                                        
                                    }
                                }                      
                            }
                            
                            // ------- Real-packets -----
                            
                            
                            // New Stokes Parameter (unpolarized!) and direction in the CMF (istropic!)
                            I=1; Q=0; U=0;
                            u = 2 * gsl_rng_uniform(rng) - 1;
                            v = 2 * PI * gsl_rng_uniform(rng);
                            n_cmf[0] = sqrt(1-pow(u,2.)) * cos(v);
                            n_cmf[1] = sqrt(1-pow(u,2.)) * sin(v);
                            n_cmf[2] = u;
                            
                            // Compute the new packet direction in the RF
                            vel_rev[0] = - vel[0] ;
                            vel_rev[1] = - vel[1] ;
                            vel_rev[2] = - vel[2] ;
                            aberr(n_cmf,vel_rev,n_rf);
                            
                            // Compute new RF frequency
                            nu_rf = nu_cmf / doppler(n_rf,vel);
                            
                            // ---- Update packet status ----
                            
                            color[which_line] = -1 ;   // We will never interact with this line again
                            
                            for (ind=0; ind<NLINE; ind++) s_l[ind] = s_l[ind] - s ;   // Update sobolev point
                            
                        }

                    } // end while loop --> new interaction

                } // end for loop --> new packet
                
                fq=fopen(filespec,"w");
                fqdirect = fopen(filedirec, "w");
                fqother = fopen(fileother, "w");
                fprintf(fq,"%d \n%d\n%d %g %g \n",Nobs,Nwave,Ntime,ti/DAY,tf/DAY);
                fprintf(fqdirect,"%d \n%d\n%d %g %g \n",Nobs,Nwave,Ntime,ti/DAY,tf/DAY);
                fprintf(fqother,"%d \n%d\n%d %g %g \n",Nobs,Nwave,Ntime,ti/DAY,tf/DAY);
                
                for (l = 0; l < Nobs; l++) {

                    for (ind_wave = 0; ind_wave < Nwave ;ind_wave++) {

                        fprintf(fq,"%3.3e \t",wavei+(ind_wave+0.5)*step_wave);
                        fprintf(fqdirect,"%3.3e \t",wavei+(ind_wave+0.5)*step_wave);
                        fprintf(fqother,"%3.3e \t",wavei+(ind_wave+0.5)*step_wave);

                        for (ind_t = 0; ind_t < Ntime ;ind_t++) {
                            
                            ind_spec = ind_wave + ind_t * Nwave + l * Ntime * Nwave ;

                            Iwave[ind_spec] = Iwave_direct[ind_spec] + Iwave_sc[ind_spec] + Iwave_abs[ind_spec] + Iwave_line[ind_spec] ;
                            Qwave[ind_spec] = Qwave_direct[ind_spec] + Qwave_sc[ind_spec] + Qwave_abs[ind_spec] + Qwave_line[ind_spec] ;
                            Uwave[ind_spec] = Uwave_direct[ind_spec] + Uwave_sc[ind_spec] + Uwave_abs[ind_spec] + Uwave_line[ind_spec] ;

                            if (lf_contr == 1) {

                                Iwave_lf[ind_spec] = Iwave_lf_direct[ind_spec] + Iwave_lf_sc[ind_spec] + Iwave_lf_abs[ind_spec] + Iwave_lf_line[ind_spec] ;
                                Qwave_lf[ind_spec] = Qwave_lf_direct[ind_spec] + Qwave_lf_sc[ind_spec] + Qwave_lf_abs[ind_spec] + Qwave_lf_line[ind_spec] ;
                                Uwave_lf[ind_spec] = Uwave_lf_direct[ind_spec] + Uwave_lf_sc[ind_spec] + Uwave_lf_abs[ind_spec] + Uwave_lf_line[ind_spec] ;

                                fprintf(fq,"%3.5g \t %3.5g \t %3.5g \t %3.5g \t %3.5g \t %3.5g \t",
                                    Iwave[ind_spec],Iwave_lf[ind_spec],Qwave[ind_spec],Qwave_lf[ind_spec],
                                    Uwave[ind_spec],Uwave_lf[ind_spec]);
                                fprintf(fqdirect,"%3.5g \t %3.5g \t %3.5g \t %3.5g \t %3.5g \t %3.5g \t",
                                    Iwave_direct[ind_spec],Iwave_lf_direct[ind_spec],Qwave_direct[ind_spec],Qwave_lf_direct[ind_spec],
                                    Uwave_direct[ind_spec],Uwave_lf_direct[ind_spec]);
                                fprintf(fqother,"%3.5g \t %3.5g \t %3.5g \t %3.5g \t %3.5g \t %3.5g \t",
                                    Iwave[ind_spec] - Iwave_direct[ind_spec],Iwave_lf[ind_spec] - Iwave_lf_direct[ind_spec],Qwave[ind_spec] - Qwave_direct[ind_spec],Qwave_lf[ind_spec] - Qwave_lf_direct[ind_spec],
                                    Uwave[ind_spec] - Uwave_direct[ind_spec],Uwave_lf[ind_spec] - Uwave_lf_direct[ind_spec]);

                            }

                            else {
                                fprintf(fq,"%3.5g \t %3.5g \t %3.5g \t",Iwave[ind_spec],Qwave[ind_spec],Uwave[ind_spec]);
                                fprintf(fqdirect,"%3.5g \t %3.5g \t %3.5g \t",Iwave_direct[ind_spec],Qwave_direct[ind_spec],Uwave_direct[ind_spec]);
                                fprintf(fqother,"%3.5g \t %3.5g \t %3.5g \t",Iwave[ind_spec] - Iwave_direct[ind_spec],Qwave[ind_spec] - Qwave_direct[ind_spec],Uwave[ind_spec] - Uwave_direct[ind_spec]);
                            }


                        }      

                        fprintf(fq,"\n") ;
                        fprintf(fqdirect,"\n") ;
                        fprintf(fqother,"\n") ;
                    }

                }

                end = clock();
                time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            
                fprintf(finfo,"Timestep %d out of %d finished (%.2f days). Time: %3.3f \n",itime+1,Ntime,(t+step_t/2.)/DAY,time_spent);

                fclose(fq);
                fclose(fqdirect);
                fclose(fqother);
                fclose(finfo);

            } // end time loop -> new timestep
            
            finfo=fopen(fileinfo,"a+");

            fprintf(finfo,"\n%.3g real packets thermalised in the inner region",N_real_inner);
            fprintf(finfo,"\n%.3g real packets were killed as they travelled in a high-opacity cell \n",N_real_highop);

            for (l = 0; l < Nobs; l++) { 

                Narrived = N_direct[l] + N_abs[l] + N_sc[l] + N_line[l] ;
                
                fprintf(finfo,"\n ---------- OBSERVER %d (costh = %.2f, phi = %.1f) -----------",l+1,cos(tobs[l]),phiobs[l]);
                fprintf(finfo,"\n%.3g virtual packets arrived at the observer\n",Narrived);
                fprintf(finfo,"- %.3g (%.2f%%) arrived directly after creation\n",N_direct[l],100*N_direct[l]/Narrived);
                fprintf(finfo,"- %.3g (%.2f%%) were scattered by a continuum interaction\n",N_abs[l],100*N_abs[l]/Narrived);
                fprintf(finfo,"- %.3g (%.2f%%) were scattered by an electron\n",N_sc[l],100*N_sc[l]/Narrived);
                fprintf(finfo,"- %.3g (%.2f%%) were scattered by a line\n",N_line[l],100*N_line[l]/Narrived);
                fprintf(finfo,"%.3g virtual packets thermalised in the inner region\n",N_inner[l]);
            }

            fclose(fichier);

            end = clock();
            time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
            fprintf(finfo,"\nTotal time: %3.3f \n\n",time_spent);
            
            fclose(finfo);

            free(Iwave_direct); free(Qwave_direct); free(Uwave_direct); 
            free(Iwave_sc); free(Qwave_sc); free(Uwave_sc);
            free(Iwave_abs); free(Qwave_abs); free(Uwave_abs);
            free(Iwave_line); free(Qwave_line); free(Uwave_line);
            free(Iwave); free(Qwave); free(Uwave);

            if (lf_contr == 1){

                free(Iwave_lf_direct); free(Qwave_lf_direct); free(Uwave_lf_direct);
                free(Iwave_lf_sc); free(Qwave_lf_sc); free(Uwave_lf_sc);
                free(Iwave_lf_abs); free(Qwave_lf_abs); free(Uwave_lf_abs);
                free(Iwave_lf_line); free(Qwave_lf_line); free(Uwave_lf_line);
                free(Iwave_lf); free(Qwave_lf); free(Uwave_lf);
            }

        }

    } // end of simulation -> new simulation

return 0;

}

