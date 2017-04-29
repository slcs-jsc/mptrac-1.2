/*! 
  \file
  Lagrangian particle dispersion model.
*/

#include "libtrac.h"

#ifdef MPI
#include "mpi.h"
#endif

/* ------------------------------------------------------------
   Global variables...
   ------------------------------------------------------------ */

/*! Timer for total runtime. */
int timer_total;

/*! Timer for physics calculations. */
int timer_phys;

/*! Timer for file input. */
int timer_input;

/*! Timer for file output. */
int timer_output;

int timer_advection[128];
int timer_diffusion_meso[128];
int timer_diffusion_basic[128];
int timer_variation;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Set simulation time interval. */
void init_simtime(
  ctl_t * ctl,
  atm_t * atm);

/*! Calculate advection of air parcels. */
void module_advection(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt);

void module_advection_petterssen1( ctl_t * ctl,  met_t * met0,  met_t * met1,  atm_t * atm,  int ip,  double dt);
int module_advection_petterssen( ctl_t * ctl,  met_t * met0,  met_t * met1,  atm_t * atm,  int ip,  double dt);
void module_advection_rungekutta( ctl_t * ctl,  met_t * met0,  met_t * met1,  atm_t * atm,  int ip,  double dt);
void module_advection_midpoint( ctl_t * ctl,  met_t * met0,  met_t * met1,  atm_t * atm,  int ip,  double dt);
void module_advection_euler( ctl_t * ctl,  met_t * met0,  met_t * met1,  atm_t * atm,  int ip,  double dt);
void module_advection_heun( ctl_t * ctl,  met_t * met0,  met_t * met1,  atm_t * atm,  int ip,  double dt);
void module_advection_heun_split( ctl_t * ctl,  met_t * met0,  met_t * met1,  atm_t * atm,  int ip,  double dt);
  

double bilinear(met_t *met0, double lon, double lat);

void module_diffusion_basic(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt,
  int rngid);

void module_diffusion_meso_grid(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt,
  int rngid);

void module_diffusion_meso_step(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt,
  int rngid);

/*! Calculate exponential decay of particle mass. */
void module_decay(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt);

/*! Check position of air parcels. */
void module_position(
  met_t * met,
  atm_t * atm,
  int ip);

/*! Calculate sedimentation of air parcels. */
void module_sedi(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt);

/*! Write simulation output. */
void write_output(
  const char *dirname,
  ctl_t * ctl,
  atm_t * atm,
  double t,
  int force);


/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

  gsl_rng *rng[NTHREADS*20];
#define RNG_UP_OFFSET NTHREADS
#define RNG_VP_OFFSET NTHREADS*2
#define RNG_WP_OFFSET NTHREADS*3
#define RNG_LON_OFFSET NTHREADS*4
#define RNG_LAT_OFFSET NTHREADS*5
#define RNG_P_OFFSET NTHREADS*6
  
int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  atm_t *atm;

  met_t *met0, *met1;


  FILE *dirlist;
  
  char buf[LEN];

  char dirname[LEN];

  double dt, t, t0;

  int gridmeso, iterations, i, ip=0, ntask = 0, rank = 0, size = 1;
  
  size_t petterssen_iterations=0, petterssen_calls=0;

#ifdef MPI
  /* Initialize MPI... */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <dirlist> <ctl> <atm_in> <metbase>");

  /* Open directory list... */
  if (!(dirlist = fopen(argv[1], "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(dirlist, "%s", dirname) != EOF) {

    /* MPI parallelization... */
    if ((++ntask) % size != rank)
      continue;

    /* ------------------------------------------------------------
       Initialize model run...
       ------------------------------------------------------------ */

    /* Create timers... */
    timer_total = CREATE_TIMER("total");
    timer_phys = CREATE_TIMER("physics");
    timer_input = CREATE_TIMER("input");
    timer_variation = CREATE_TIMER("variation");
    timer_output = CREATE_TIMER("output");

    /* Start timer for total runtime... */
    START_TIMER(timer_total);

    /* Allocate... */
    ALLOC(atm, atm_t, 1);
    ALLOC(met0, met_t, 1);
    ALLOC(met1, met_t, 1);


    /* Read control parameters... */
    read_ctl(dirname, argv[2], argc, argv, &ctl);

  
#ifdef RK
  printf("Advection-Model: Runge-Kutta %i\n", RK);
#endif
#ifdef HEUN2
  printf("Advection-Model: RK-Heun2\n");
#endif
#ifdef HEUN3
  printf("Advection-Model: RK-Heun3\n");
#endif
#ifdef HEUN_EXPLICIT
  printf("Advection-Model: Heun\n");
#endif
#ifdef HEUN_SPLIT
  printf("Advection-Model: Heun (split)\n");
#endif
#ifdef PETTERSSEN
  printf("Advection-Model: PETTERSSEN\n");
#endif
#ifdef PETTERSSEN1
  printf("Advection-Model: PETTERSSEN1\n");
#endif
#ifdef MIDPOINT
  printf("Advection-Model: MIDPOINT\n");
#endif
#ifdef EULER
  printf("Advection-Model: EULER\n");
#endif

    /* Initialize random number generators...*/ 
    gsl_rng_env_setup();
    for (i = 0; i < NTHREADS*20; i++) {
      rng[i] = gsl_rng_alloc(gsl_rng_taus);
    }
    
    if(ctl.rng_seed != 0) 
      for(ip = 0;ip<NTHREADS*20;ip++)
	gsl_rng_set(rng[ip], (unsigned int)(ctl.rng_seed+ip));
    
    

    /* Read atmospheric data... */
    START_TIMER(timer_input);
    read_atm(dirname, argv[3], atm, &ctl);
    STOP_TIMER(timer_input);
    
    /* Setup diffuison */
#ifdef OLD_DIFFUSION
    gridmeso = 1==2;
    printf("Using the old Diffusion model for meso and basic\n");
#else
    gridmeso = atm->np * ctl.dt_met / ctl.dt_mod > met0->nx*met0->ny*met0->np;
    if(gridmeso) printf("Using the new Diffusion model for meso and basic %g - %u\n",atm->np * ctl.dt_met/ ctl.dt_mod,met0->nx*met0->ny*met0->np);
    else printf("Using the new Diffusion model for basic only\n");
#endif
  
    /* Get simulation time interval... */
    init_simtime(&ctl, atm);

    /* Write initial output... */
    START_TIMER(timer_output);
    write_output(dirname, &ctl, atm, ctl.t_start, 1);
    STOP_TIMER(timer_output);


    /* ------------------------------------------------------------
       Loop over timesteps...
       ------------------------------------------------------------ */

    /* Get rounded start time... */
    if (ctl.direction == 1)
      t0 = floor(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;
    else
      t0 = ceil(ctl.t_start / ctl.dt_mod) * ctl.dt_mod;

    /* Loop over timesteps... */
    for (t = t0; ctl.direction * (t - ctl.t_stop) < ctl.dt_mod;
	 t += ctl.direction * ctl.dt_mod) {

      /* Adjust length of final time step... */
      if (ctl.direction * (t - ctl.t_stop) > 0)
	t = ctl.t_stop;

      /* Get meteorological data... */
      START_TIMER(timer_input);
      get_met(t, ctl.direction, argv[4], ctl.dt_met, ctl.red_met, met0,
	      met1);
#ifndef OLD_DIFFUSION      
    START_TIMER(timer_variation);
      if(!met1->meso && ctl.turb_meso > 0)
	updateVariation(met0, met1);
    STOP_TIMER(timer_variation);
#endif
      STOP_TIMER(timer_input);

      /* Loop over air parcels... */
      START_TIMER(timer_phys);
#pragma omp parallel for default(shared) private(dt,ip) reduction(+:petterssen_iterations, petterssen_calls)
      for (ip = 0; ip < atm->np; ip++)
	if ((ctl.direction * (atm->time[ip] - ctl.t_start) >= 0
	     && ctl.direction * (atm->time[ip] - ctl.t_stop) <= 0
	     && ctl.direction * (atm->time[ip] - t) < 0)) {

	  /* Set time step... */
	  dt = t - atm->time[ip];
	  /* Calculate advection... */
	  
 /* printf("\npre  %8.2f %g,%g\n", atm->time[ip], atm->lon[ip], atm->lat[ip]);
          iterations=1;
          if(fabs(atm->lat[ip]) > 86) iterations*=2;          
	  if(fabs(atm->lat[ip]) > 88) iterations*=2;     
	  if(fabs(atm->lat[ip]) > 89) iterations*=2;
 
          dt=dt/(double)iterations;
	  while((--iterations) >= 0) {*/

	
#ifdef ADVECTION_OFF
	  atm->time[ip]+=dt;
#else
	    
#ifndef PETTERSSEN
	  module_advection(&ctl, met0, met1, atm, ip, dt);
#else
	petterssen_iterations += (size_t)module_advection_petterssen(&ctl, met0, met1, atm, ip, dt);
	petterssen_calls++;  
#endif
#endif

#ifdef DIFFUSION_ON
	  /* Calculate diffusion... */
	  module_diffusion_basic(&ctl, met0, met1, atm, ip, dt,
			   omp_get_thread_num());
	  
	  if(gridmeso) {
	    module_diffusion_meso_grid(&ctl, met0, met1, atm, ip, dt,
			   omp_get_thread_num());
	  } else {
	    module_diffusion_meso_step(&ctl, met0, met1, atm, ip, dt,
			   omp_get_thread_num());
	  }
#endif

	  /* Calculate sedimentation... */
	  module_sedi(&ctl, atm, ip, dt);

	  /* Check position... */
	  module_position(met0, atm, ip);

	  /* Calculate decay of mass... */
	  module_decay(&ctl, atm, ip, dt);
	  }
	/*}*/
      STOP_TIMER(timer_phys);

      
      
      /* Write output... */
      START_TIMER(timer_output);
      write_output(dirname, &ctl, atm, t, t == ctl.t_stop);
      STOP_TIMER(timer_output);
    }

    /* ------------------------------------------------------------
       Finalize model run...
       ------------------------------------------------------------ */

#ifdef PETTERSSEN
    printf("# n=%lu petterssen_iterations = %lu\n", (long int)petterssen_calls , (long int)petterssen_iterations);
    petterssen_calls=0;
    petterssen_iterations=0;
#endif
    
    /* Free random number generators... */
    for (i = 0; i < NTHREADS; i++)
      gsl_rng_free(rng[i]);

    /* Free... */
    free(atm);
    free(met0);
    free(met1);

    /* Stop timer for total runtime... */
    STOP_TIMER(timer_total);

    /* Report timers... */
    printf("# cores = %i\n", omp_get_max_threads());
    PRINT_TIMER(timer_total);
    PRINT_TIMER(timer_phys);
    PRINT_TIMER(timer_input);
    PRINT_TIMER(timer_variation);
    PRINT_TIMER(timer_output);
  }

#ifdef MPI
  /* Finalize MPI... */
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void init_simtime(
  ctl_t * ctl,
  atm_t * atm) {

  /* Set inital and final time... */
  if (ctl->direction == 1) {
    if (ctl->t_start < -1e99)
      ctl->t_start = gsl_stats_min(atm->time, 1, (size_t) atm->np);
    if (ctl->t_stop < -1e99)
      ctl->t_stop = gsl_stats_max(atm->time, 1, (size_t) atm->np);
  } else if (ctl->direction == -1) {
    if (ctl->t_stop < -1e99)
      ctl->t_stop = gsl_stats_min(atm->time, 1, (size_t) atm->np);
    if (ctl->t_start < -1e99)
      ctl->t_start = gsl_stats_max(atm->time, 1, (size_t) atm->np);
  }

  /* Check time... */
  if (ctl->direction * (ctl->t_stop - ctl->t_start) <= 0)
    ERRMSG("Nothing to do!");
}




/*****************************************************************************/

void module_advection(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {
  
#ifdef RK
  module_advection_rungekutta(ctl, met0, met1, atm, ip, dt);
#else
#ifdef PETTERSSEN
  module_advection_petterssen(ctl, met0, met1, atm, ip, dt);
#else
#ifdef PETTERSSEN1
  module_advection_petterssen1(ctl, met0, met1, atm, ip, dt);
#else
#ifdef MIDPOINT
  module_advection_midpoint(ctl, met0, met1, atm, ip, dt);
#else
#ifdef HEUN_EXPLICIT
  module_advection_heun(ctl, met0, met1, atm, ip, dt);
#else
#ifdef HEUN_SPLIT
  module_advection_heun_split(ctl, met0, met1, atm, ip, dt);
#else
  module_advection_euler(ctl, met0, met1, atm, ip, dt);
#endif
#endif
#endif
#endif
#endif
#endif
}




void module_advection_petterssen1(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {
  double dx, dy, dz, t0, t1, v[3], vn[3], x[3], xout[3];
  
    /* Copy air parcel data... */
    x[0]=atm->lon[ip];
    x[1]=atm->lat[ip];
    x[2]=atm->p[ip];
    
    /* Interpolate meteorological data... */
    intpol_met_time(met0, met1, atm->time[ip], x[2], x[0], x[1], &t0,
		    &v[0], &v[1], &v[2]);
    
    /* Advection... */
    dx=dt * v[0];
    dy=dt * v[1];
    dz=dt * v[2];
    
    /* Get new position... */
    xout[0]=x[0]+dx2deg(dx/1000., x[1]);
    xout[1]=x[1]+dy2deg(dy/1000.);
    xout[2]=x[2]+dz;

    /* Petterssen scheme... */
	  
    /* Interpolate meteorological data... */
    intpol_met_time(met0, met1, atm->time[ip] + dt, xout[2], xout[0], xout[1], 
		    &t1, &vn[0], &vn[1], &vn[2] );
      
    /* Advection... */
    dx=0.5*dt*(v[0]+vn[0]);
    dy=0.5*dt*(v[1]+vn[1]);
    dz=0.5*dt*(v[2]+vn[2]);
      
    /* Save new position... */
    atm->time[ip] += dt;
    atm->lon[ip] = x[0]+dx2deg(dx/1000., (xout[1]+x[1])/2.0 );
    atm->lat[ip] = x[1]+dy2deg(dy/1000.);
    atm->p[ip] = x[2]+dz;
  
    /* Extrapolate temperature... */
    if (ctl->qnt_temp >= 0)
      atm->q[ctl->qnt_temp][ip] = (t0+t1)/2;
}

#define PETTERSSEN_MAX_ITER 6
int module_advection_petterssen(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {
  double dx, dy, dz, t1, v[3], vn[3], x[3], xout[3], oldd = 0;
  
  int n = 0;
    /* Copy air parcel data... */
    x[0]=atm->lon[ip];
    x[1]=atm->lat[ip];
    x[2]=atm->p[ip];
    
    /* Interpolate meteorological data... */
    intpol_met_time(met0, met1, atm->time[ip], x[2], x[0], x[1], &t1,
		    &v[0], &v[1], &v[2]);
    
    /* Advection... */
    dx=dt * v[0];
    dy=dt * v[1];
    dz=dt * v[2];
    
    /* Get new position... */
    xout[0]=x[0]+dx2deg(dx/1000., x[1]);
    xout[1]=x[1]+dy2deg(dy/1000.);
    xout[2]=x[2]+dz;

    /* Petterssen scheme... */
    while(fabs(oldd-(fabs(dx)+fabs(dy)+fabs(dz))) > 1e-7
	  && n < PETTERSSEN_MAX_ITER) {
      oldd=fabs(dx)+fabs(dy)+fabs(dz);
	  
      /* Interpolate meteorological data... */
      intpol_met_time(met0, met1, atm->time[ip] + dt, xout[2], xout[0], xout[1], 
		      &t1, &vn[0], &vn[1], &vn[2] );
      
      /* Advection... */
      dx=0.5*dt*(v[0]+vn[0]);
      dy=0.5*dt*(v[1]+vn[1]);
      dz=0.5*dt*(v[2]+vn[2]);
      
      /* Get new position... */
      xout[0]=x[0]+dx2deg(dx/1000., (xout[1]+x[1])/2);
      xout[1]=x[1]+dy2deg(dy/1000.);
      xout[2]=x[2]+dz;
      
      /* Increase iteration counter... */
      n++;     
    }
    
  /* Save new position... */
  atm->time[ip] += dt;
  atm->lon[ip] = xout[0];
  atm->lat[ip] = xout[1];
  atm->p[ip] = xout[2];
  
  /* Extrapolate temperature... */
  if (ctl->qnt_temp >= 0)
    atm->q[ctl->qnt_temp][ip] = t1;
  
  return n;
}


/*****************************************************************************/


void module_advection_midpoint(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {

  double v[3], xm[3], t0, t1;

  /* Interpolate meteorological data... */
  intpol_met_time(met0, met1, atm->time[ip], atm->p[ip],
		  atm->lon[ip], atm->lat[ip], 
		  &t0, &v[0], &v[1], &v[2]);

  /* Get position of the mid point... */
  xm[0] = atm->lon[ip] + dx2deg(0.5 * dt * v[0] / 1000., atm->lat[ip]);
  xm[1] = atm->lat[ip] + dy2deg(0.5 * dt * v[1] / 1000.);
  xm[2] = atm->p[ip] + 0.5 * dt * v[2];

  /* Interpolate meteorological data for mid point... */
  intpol_met_time(met0, met1, atm->time[ip] + 0.5 * dt,
		  xm[2], xm[0], xm[1], 
		  &t1, &v[0], &v[1], &v[2]);

  /* Save new position... */
  atm->time[ip] += dt;
  atm->lon[ip] += dx2deg(dt * v[0] / 1000., xm[1]);
  atm->lat[ip] += dy2deg(dt * v[1] / 1000.);
  atm->p[ip] += dt * v[2];

  /* Extrapolate temperature... */
  if (ctl->qnt_temp >= 0)
    atm->q[ctl->qnt_temp][ip] = 2. * t1 - t0;
}


/*****************************************************************************/

void module_advection_euler(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {

  double t1, v[3], x[3];

  /* Copy air parcel data... */
  x[0] = atm->lon[ip];
  x[1] = atm->lat[ip];
  x[2] = atm->p[ip];

  /* Interpolate meteorological data... */
  intpol_met_time(met0, met1, atm->time[ip],
		  x[2], x[0], x[1], &t1,
		  &v[0], &v[1], &v[2]);

  /* Save new position... */
  atm->time[ip] += dt;
  atm->lon[ip] = x[0] + dx2deg(dt * v[0] / 1000., x[1]);
  atm->lat[ip] = x[1] + dy2deg(dt * v[1] / 1000.);
  atm->p[ip] = x[2] + dt * v[2];

  
  
  /* Extrapolate temperature... */
  if (ctl->qnt_temp >= 0)
    atm->q[ctl->qnt_temp][ip] = t1;
}


/*****************************************************************************/

#ifdef RK
void module_advection_rungekutta(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {
  
#if RK == 4  
  double b[] = {1/6.0, 1/3.0, 1/3.0, 1/6.0};
  double c[] = {0, 0.5, 0.5, 1};
  double a[4][4] = {{0, 0, 0, 0}, {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}};
#endif
    
#if RK == 3 
#ifdef HEUN3
  double b[] = {0.25, 0, 0.75};
  double c[] = {0, 1/3.0, 2/3.0};
  double a[3][3] = {{0, 0, 0}, {1/3.0, 0, 0}, {0, 2/3.0, 0}};
#else
  double b[] = {1/6.0, 4/6.0, 1/6.0};
  double c[] = {0, 0.5, 1};
  double a[3][3] = {{0, 0, 0}, {0.5, 0, 0}, {-1, 2, 0}};
#endif
#endif
  
#if RK == 2 
#ifdef HEUN2
  double b[] = {0.5, 0.5};
  double c[] = {0, 1};
  double a[2][2] = {{0, 0}, {1, 0}};
#else
  double b[] = {0, 1};
  double c[] = {0, 0.5};
  double a[2][2] = {{0, 0}, {0.5, 0}};
#endif
#endif
  
#if RK == 1  
  double b[] = {1};
  double c[] = {0};
#else
  int j;
#endif
  
  
  double sumak[3], sumbk[3];
  double k[3][RK];
  int i;
  double meanlat=0, xout[3], startx[3],  t1;
  
  
  /* Copy air parcel data... */
  startx[0] = atm->lon[ip];
  xout[1] = startx[1] = atm->lat[ip];
  startx[2] = atm->p[ip];
  
  
  sumbk[0] = sumbk[1] = sumbk[2] = 0;
  
  for(i=0; i<RK; i++) {
    sumak[0] = sumak[1] = sumak[2] = 0;
#if RK != 1  
    for(j=0; j<i; j++) {
      sumak[0] += a[i][j] * k[0][j];
      sumak[1] += a[i][j] * k[1][j];
      sumak[2] += a[i][j] * k[2][j];
    }
#endif

    xout[0] = startx[0] + dx2deg(dt * sumak[0] / 1000., xout[1]);
    xout[1] = startx[1] + dy2deg(dt * sumak[1] / 1000.);
    xout[2] = startx[2] +  dt * sumak[2];
    meanlat += xout[1] * b[i];

    intpol_met_time(met0, met1,
		    atm->time[ip] + c[i] * dt,
		    xout[2],
		    xout[0],
		    xout[1],
		    &t1, &k[0][i], &k[1][i], &k[2][i]);
      
    sumbk[0] += b[i] * k[0][i];
    sumbk[1] += b[i] * k[1][i];
    sumbk[2] += b[i] * k[2][i];
  }
  /* Save new position... */
  atm->time[ip] += dt;
  atm->lon[ip] = startx[0] + dx2deg(dt * sumbk[0] / 1000., meanlat);
  atm->lat[ip] = startx[1] + dy2deg(dt * sumbk[1] / 1000.);
  atm->p[ip] = startx[2] + dt * sumbk[2];

  /* Extrapolate temperature... */
  if (ctl->qnt_temp >= 0)
    atm->q[ctl->qnt_temp][ip] = t1 ;
}
#endif


/*****************************************************************************/
void module_advection_heun(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {

  
  double xout[3], startx[3], t1, t2;
  
  double wroot[3];
  double wstep[3];
  
  
  /* Copy air parcel data... */
  startx[0] = atm->lon[ip];
  startx[1] = atm->lat[ip];
  startx[2] = atm->p[ip];
  
  /* wind at start point c=0 */
  intpol_met_time(met0, met1,
		    atm->time[ip],
		    startx[2],
		    startx[0],
		    startx[1],
		    &t1, &wroot[0], &wroot[1], &wroot[2]);
  
  /* help point full euler step */
  xout[0] = startx[0] + dx2deg(dt * wroot[0] / 1000., startx[1]);
  xout[1] = startx[1] + dy2deg(dt * wroot[1] / 1000.);
  xout[2] = startx[2] + dt * wroot[2];
  
  /* wind at step point c=1 */
  intpol_met_time(met0, met1,
		    atm->time[ip] + dt,
		    xout[2],
		    xout[0],
		    xout[1],
		    &t2, &wstep[0], &wstep[1], &wstep[2]);
  
  /* Save new position... b[0]=b[1]=0.5 */
  atm->time[ip] += dt;
  atm->lon[ip] = startx[0] + dx2deg(dt * (wroot[0]+wstep[0]) / 2000., (startx[1]+xout[1])/2);
  atm->lat[ip] = startx[1] + dy2deg(dt * (wroot[1]+wstep[1]) / 2000.);
  atm->p[ip]   = startx[2] + dt * (wroot[2]+wstep[2]) / 2.0;

  /* Extrapolate temperature... */
  if (ctl->qnt_temp >= 0)
    atm->q[ctl->qnt_temp][ip] = (t1+t2)/2.0 ;
}


void module_advection_heun_split(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt) {

  double xout[3], startx[3], t1, t2;
  
  double wroot[3];
  double wstep[3];
  
  
  /* Copy air parcel data... */
  startx[0] = atm->lon[ip];
  startx[1] = atm->lat[ip];
  startx[2] = atm->p[ip];
  
  /* wind at start point c=0 */
  intpol_met_time(met0, met1,
		    atm->time[ip],
		    startx[2],
		    startx[0],
		    startx[1],
		    &t1, &wroot[0], &wroot[1], &wroot[2]);
  
  /* help point full euler step */
  xout[0] = startx[0] + dx2deg(dt * wroot[0] / 1000., startx[1]);
  xout[1] = startx[1] + dy2deg(dt * wroot[1] / 1000.);
  xout[2] = startx[2] + dt * wroot[2];
  
  /* wind at step point c=1 */
  intpol_met_time(met0, met1,
		    atm->time[ip] + dt,
		    xout[2],
		    xout[0],
		    xout[1],
		    &t2, &wstep[0], &wstep[1], &wstep[2]);
  
  /* Save new position... b[0]=b[1]=0.5 */
  atm->time[ip] += dt;
  atm->lon[ip] = startx[0] + dx2deg(dt * wroot[0] / 2000., startx[1]) + dx2deg(dt * wstep[0] / 2000., xout[1]);
  atm->lat[ip] = startx[1] + dy2deg(dt * (wroot[1]+wstep[1]) / 2000.);
  atm->p[ip]   = startx[2] + dt * (wroot[2]+wstep[2]) / 2.0;

  /* Extrapolate temperature... */
  if (ctl->qnt_temp >= 0)
    atm->q[ctl->qnt_temp][ip] = (t1+t2)/2.0 ;
}


/*****************************************************************************/

void module_decay(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt) {

  /* Calculate exponential decay... */
  if (ctl->t12 > 0 && ctl->qnt_mass >= 0)
    atm->q[ctl->qnt_mass][ip] *= exp(-dt / (ctl->t12 * 86400. / log(2.)));
}

/*****************************************************************************/


double bilinear(met_t *met0, double lon, double lat) {
  int ix = locate(met0->lon, met0->nx, lon);
  int iy = locate(met0->lat, met0->ny, lat);  
  double iplat0, iplat1, wlat, wlon;
  /* latitude interpolation */
  wlat=(met0->lat[iy+1]-lat) / (met0->lat[iy+1]-met0->lat[iy]);
  iplat0= wlat * met0->tp[ix][iy]   + (1-wlat) * met0->tp[ix][iy+1];
  iplat1= wlat * met0->tp[ix+1][iy] + (1-wlat) * met0->tp[ix+1][iy+1];
  
  /* longitude interpolation */
  wlon=(met0->lon[ix+1]-lon) / (met0->lon[ix+1]-met0->lon[ix]);
  return wlon*iplat0 + (1-wlon) * iplat1;
}


void module_diffusion_basic (
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt,
  int rngid) {

#ifndef OLD_DIFFUSION
  double  wmet0; 
  double wt, ws;
  double tp;
  
  /* interpolate in time */
#ifdef ADVECTION_OFF
  tp=P(12);
#else
  wmet0 = (met1->time - atm->time[ip]) / (met1->time - met0->time);
  tp = wmet0 * bilinear(met0, atm->lon[ip], atm->lat[ip]) 
    + (1-wmet0) * bilinear(met1, atm->lon[ip], atm->lat[ip]);
#endif
  if(atm->p[ip] < tp) {
    if(Z(atm->p[ip]) >= Z(tp)+1) {
    /* Vertical turbulent diffusion in the stratosphere */
      if (ctl->turb_dz > 0)
	atm->p[ip] += dz2dp(gsl_ran_gaussian_ziggurat(rng[rngid+RNG_P_OFFSET], sqrt(ctl->turb_dz*2 * fabs(dt))) / 1000., atm->p[ip]);
    } else {
      /* Combined turbulent diffusion above tropopause */
      ws=(Z(atm->p[ip])-Z(tp));
      wt=1-ws;
      if (ctl->turb_dx > 0) {
	atm->lon[ip] += dx2deg(gsl_ran_gaussian_ziggurat(rng[rngid+RNG_LON_OFFSET], sqrt(wt * ctl->turb_dx*2 * fabs(dt))) / 1000., atm->lat[ip]);
	atm->lat[ip] += dy2deg(gsl_ran_gaussian_ziggurat(rng[rngid+RNG_LAT_OFFSET], sqrt(wt * ctl->turb_dx*2 * fabs(dt))) / 1000.);
      }
      if (ctl->turb_dz > 0)
	atm->p[ip] += dz2dp(gsl_ran_gaussian_ziggurat(rng[rngid+RNG_P_OFFSET], sqrt(ws * ctl->turb_dz*2 * fabs(dt))) / 1000., atm->p[ip]);
    }      
  } else {
    /* Horizontal turbulent diffusion in the troposphere */
    if (ctl->turb_dx > 0) {
      atm->lon[ip] += dx2deg(gsl_ran_gaussian_ziggurat(rng[rngid+RNG_LON_OFFSET], sqrt(ctl->turb_dx*2 * fabs(dt))) / 1000., atm->lat[ip]);
      atm->lat[ip] += dy2deg(gsl_ran_gaussian_ziggurat(rng[rngid+RNG_LAT_OFFSET], sqrt(ctl->turb_dx*2 * fabs(dt))) / 1000.);
    }
  }
#else
  /* Horizontal turbulent diffusion... */
  if (ctl->turb_dx > 0) {
    atm->lon[ip]
      += dx2deg(gsl_ran_gaussian_ziggurat(rng[rngid], sqrt(ctl->turb_dx*2 / fabs(dt)))
		/ 1000., atm->lat[ip]);
    atm->lat[ip]
      += dy2deg(gsl_ran_gaussian_ziggurat(rng[rngid], sqrt(ctl->turb_dx*2 / fabs(dt)))
		/ 1000.);
  }

  /* Vertical turbulent diffusion... */
  if (ctl->turb_dz > 0)
    atm->p[ip]
      += dz2dp(gsl_ran_gaussian_ziggurat(rng[rngid], sqrt(ctl->turb_dz*2 / fabs(dt)))
	       / 1000., atm->p[ip]);
#endif
}

void module_diffusion_meso_grid(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt,
  int rngid) {
  
  int ix,iy,iz;
  double r,rs;
   if (ctl->turb_meso > 0) {
     /* Get indices... */
     ix = locate(met0->lon, met0->nx, atm->lon[ip]);
     iy = locate(met0->lat, met0->ny, atm->lat[ip]);
     iz = locate(met0->p, met0->np, atm->p[ip]);
     /*r=exp(-2.*abs(dt)/ctl->dt_met)*/
     r = 1 - 2 * fabs(dt) / ctl->dt_met;
     rs = sqrt(1 - r * r);

    /* Calculate mesoscale wind fluctuations... */
    atm->up[ip] =
      r * atm->up[ip] + rs * gsl_ran_gaussian_ziggurat(rng[rngid+RNG_UP_OFFSET],
						       ctl->turb_meso * met1->usig[ix][iy][iz]);
    atm->vp[ip] =
      r * atm->vp[ip] + rs * gsl_ran_gaussian_ziggurat(rng[rngid+RNG_VP_OFFSET],
						       ctl->turb_meso * met1->vsig[ix][iy][iz]);
    atm->wp[ip] =
      r * atm->wp[ip] + rs * gsl_ran_gaussian_ziggurat(rng[rngid+RNG_WP_OFFSET],
						       ctl->turb_meso * met1->wsig[ix][iy][iz]);

    /* Calculate air parcel displacement... */
    atm->lon[ip] += dx2deg(atm->up[ip] * dt / 1000., atm->lat[ip]);
    atm->lat[ip] += dy2deg(atm->vp[ip] * dt / 1000.);
    atm->p[ip] += atm->wp[ip] * dt;
   }
}

void module_diffusion_meso_step(
  ctl_t * ctl,
  met_t * met0,
  met_t * met1,
  atm_t * atm,
  int ip,
  double dt,
  int rngid) {

double r, rs, u[16], v[16], w[16], usig, vsig, wsig;

  int ix, iy, iz, ixh, iyh, izh;

  /* Calculate mesoscale velocity fluctuations... */
  if (ctl->turb_meso > 0) {

    /* Get indices... */
    ix = locate(met0->lon, met0->nx, atm->lon[ip]);
    iy = locate(met0->lat, met0->ny, atm->lat[ip]);
    iz = locate(met0->p, met0->np, atm->p[ip]);

    /* Collect local wind data... */
	if(ix<met0->nx-1) ixh=ix+1;
	else ixh=ix-1;
	if(iy<met0->ny-1) iyh=iy+1;
	else iyh=iy-1;
	if(iz<met0->np-1) izh=iz+1;
	else izh=iz-1;
	
	/* Collect local wind data... */
	u[0] = met0->u[ix][iy][iz];
	u[1] = met0->u[ixh][iy][iz];
	u[2] = met0->u[ix][iyh][iz];
	u[3] = met0->u[ixh][iyh][iz];
	u[4] = met0->u[ix][iy][izh];
	u[5] = met0->u[ixh][iy][izh];
	u[6] = met0->u[ix][iyh][izh];
	u[7] = met0->u[ixh][iyh][izh];

	v[0] = met0->v[ix][iy][iz];
	v[1] = met0->v[ixh][iy][iz];
	v[2] = met0->v[ix][iyh][iz];
	v[3] = met0->v[ixh][iyh][iz];
	v[4] = met0->v[ix][iy][izh];
	v[5] = met0->v[ixh][iy][izh];
	v[6] = met0->v[ix][iyh][izh];
	v[7] = met0->v[ixh][iyh][izh];

	w[0] = met0->w[ix][iy][iz];
	w[1] = met0->w[ixh][iy][iz];
	w[2] = met0->w[ix][iyh][iz];
	w[3] = met0->w[ixh][iyh][iz];
	w[4] = met0->w[ix][iy][izh];
	w[5] = met0->w[ixh][iy][izh];
	w[6] = met0->w[ix][iyh][izh];
	w[7] = met0->w[ixh][iyh][izh];
    

    /* Get indices... */
    ix = locate(met1->lon, met1->nx, atm->lon[ip]);
    iy = locate(met1->lat, met1->ny, atm->lat[ip]);
    iz = locate(met1->p, met1->np, atm->p[ip]);
	if(ix<met0->nx-1) ixh=ix+1;
	else ixh=ix-1;
	if(iy<met0->ny-1) iyh=iy+1;
	else iyh=iy-1;
	if(iz<met0->np-1) izh=iz+1;
	else izh=iz-1;

	/* Collect local wind data... */
	u[8] = met1->u[ix][iy][iz];
	u[9] = met1->u[ixh][iy][iz];
	u[10] = met1->u[ix][iyh][iz];
	u[11] = met1->u[ixh][iyh][iz];
	u[12] = met1->u[ix][iy][izh];
	u[13] = met1->u[ixh][iy][izh];
	u[14] = met1->u[ix][iyh][izh];
	u[15] = met1->u[ixh][iyh][izh];

	v[8] = met1->v[ix][iy][iz];
	v[9] = met1->v[ixh][iy][iz];
	v[10] = met1->v[ix][iyh][iz];
	v[11] = met1->v[ixh][iyh][iz];
	v[12] = met1->v[ix][iy][izh];
	v[13] = met1->v[ixh][iy][izh];
	v[14] = met1->v[ix][iyh][izh];
	v[15] = met1->v[ixh][iyh][izh];
  
	w[8] = met1->w[ix][iy][iz];
	w[9] = met1->w[ixh][iy][iz];
	w[10] = met1->w[ix][iyh][iz];
	w[11] = met1->w[ixh][iyh][iz];
	w[12] = met1->w[ix][iy][izh];
	w[13] = met1->w[ixh][iy][izh];
	w[14] = met1->w[ix][iyh][izh];
	w[15] = met1->w[ixh][iyh][izh];
	
    /* Get standard deviations of local wind data... */
    usig = gsl_stats_sd(u, 1, 16);
    vsig = gsl_stats_sd(v, 1, 16);
    wsig = gsl_stats_sd(w, 1, 16);

    /* Set temporal correlations for mesoscale fluctuations... */
    r = 1 - 2 * fabs(dt) / ctl->dt_met;
    rs = sqrt(1 - r * r);

    /* Calculate mesoscale wind fluctuations... */
    atm->up[ip] =
      r * atm->up[ip] + rs * gsl_ran_gaussian_ziggurat(rng[rngid+RNG_UP_OFFSET],
						       ctl->turb_meso * met1->usig[ix][iy][iz]);
    atm->vp[ip] =
      r * atm->vp[ip] + rs * gsl_ran_gaussian_ziggurat(rng[rngid+RNG_VP_OFFSET],
						       ctl->turb_meso * met1->vsig[ix][iy][iz]);
    atm->wp[ip] =
      r * atm->wp[ip] + rs * gsl_ran_gaussian_ziggurat(rng[rngid+RNG_WP_OFFSET],
						       ctl->turb_meso * met1->wsig[ix][iy][iz]);
    /* Calculate air parcel displacement... */
    atm->lon[ip] += dx2deg(atm->up[ip] * dt / 1000., atm->lat[ip]);
    atm->lat[ip] += dy2deg(atm->vp[ip] * dt / 1000.);
    atm->p[ip] += atm->wp[ip] * dt;
  }
}

/*****************************************************************************/

void module_position(
  met_t * met,
  atm_t * atm,
  int ip) {

  if(fabs(atm->lat[ip]) > 720 || fabs(atm->lon[ip]) > 360) {
    /*printf("Warning! ip=%u lon=%g lat=%g\n", ip, atm->lon[ip], atm->lat[ip]);*/
    
    atm->lon[ip] = fmod(atm->lon[ip], 360);
    atm->lat[ip] = fmod(atm->lat[ip], 360);
  }
  
  
  /* Check latitude... */
  while (atm->lat[ip] < -90 || atm->lat[ip] > 90) {
    if (atm->lat[ip] > 90) {
      atm->lat[ip] = 180 - atm->lat[ip];
      atm->lon[ip] += 180;
    }
    if (atm->lat[ip] < -90) {
      atm->lat[ip] = -180 - atm->lat[ip];
      atm->lon[ip] += 180;
    }
  }
  
  /* Check longitude... */
  while (atm->lon[ip] < -180)
    atm->lon[ip] += 360;
  while (atm->lon[ip] >= 180)
    atm->lon[ip] -= 360;

  /* Check pressure... */
  if (atm->p[ip] > met->p[0])
    atm->p[ip] = met->p[0];
  else if (atm->p[ip] < met->p[met->np - 1])
    atm->p[ip] = met->p[met->np - 1];
}

/*****************************************************************************/

void module_sedi(
  ctl_t * ctl,
  atm_t * atm,
  int ip,
  double dt) {

  /* Coefficients for Cunningham slip-flow correction (Kasten, 1968): */
  const double A = 1.249, B = 0.42, C = 0.87;

  /* Specific gas constant for dry air [J/(kg K)]: */
  const double R = 287.058;

  /* Average mass of an air molecule [kg/molec]: */
  const double m = 4.8096e-26;

  double G, K, eta, lambda, p, r_p, rho, rho_p, T, v, v_p;

  /* Check if parameters are available... */
  if (ctl->qnt_r_p < 0 || ctl->qnt_rho_p < 0 || ctl->qnt_temp < 0)
    return;

  /* Convert units... */
  p = 100 * atm->p[ip];
  T = atm->q[ctl->qnt_temp][ip];
  r_p = 1e-6 * atm->q[ctl->qnt_r_p][ip];
  rho_p = atm->q[ctl->qnt_rho_p][ip];

  /* Density of dry air... */
  rho = p / (R * T);

  /* Dynamic viscosity of air... */
  eta = 1.8325e-5 * (416.16 / (T + 120.)) * pow(T / 296.16, 1.5);

  /* Thermal velocity of an air molecule... */
  v = sqrt(8 * GSL_CONST_MKSA_BOLTZMANN * T / (M_PI * m));

  /* Mean free path of an air molecule... */
  lambda = 2 * eta / (rho * v);

  /* Knudsen number for air... */
  K = lambda / r_p;

  /* Cunningham slip-flow correction... */
  G = 1 + K * (A + B * exp(-C / K));

  /* Sedimentation (fall) velocity... */
  v_p =
    2. * gsl_pow_2(r_p) * (rho_p -
			   rho) * GSL_CONST_MKSA_GRAV_ACCEL / (9. * eta) * G;

  /* Calculate pressure change... */
  atm->p[ip] += dz2dp(v_p * dt / 1000., atm->p[ip]);
}

/*****************************************************************************/

void write_output(
  const char *dirname,
  ctl_t * ctl,
  atm_t * atm,
  double t,
  int force) {

  char filename[LEN];

  double r;

  int year, mon, day, hour, min, sec;

  /* Get time... */
  jsec2time(t, &year, &mon, &day, &hour, &min, &sec, &r);

  /* Write atmospheric data... */
  if (ctl->atm_dt_out > 0 && (fmod(t, ctl->atm_dt_out) == 0 || force)) {
    sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.%s",
	    ctl->atm_basename, year, mon, day, hour, min,
	    (ctl->atm_oformat == 0 ? "tab" : "nc"));
    write_atm(dirname, filename, atm, ctl);
  }

  /* Write CSI data... */
  if (ctl->csi_dt_update > 0 && (fmod(t, ctl->csi_dt_update) == 0 || force)) {
    sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.tab",
	    ctl->csi_basename, year, mon, day, hour, min);
    write_csi(dirname, filename, atm, ctl, t, ctl->csi_dt_update,
	      fmod(t, ctl->csi_dt_out) == 0 || force);
  }

  /* Write gridded data... */
  if (ctl->grid_dt_out > 0 && (fmod(t, ctl->grid_dt_out) == 0 || force)) {
    sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.tab",
	    ctl->grid_basename, year, mon, day, hour, min);
    write_grid(dirname, filename, atm, ctl, t, ctl->grid_dt_out);
  }

  /* Write station data... */
  if (ctl->stat_dt_out > 0 && (fmod(t, ctl->stat_dt_out) == 0 || force)) {
    sprintf(filename, "%s_%04d_%02d_%02d_%02d_%02d.tab",
	    ctl->stat_basename, year, mon, day, hour, min);
    write_station(dirname, filename, atm, ctl, t, ctl->stat_dt_out);
  }
}
