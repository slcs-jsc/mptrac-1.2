/*! 
  \file
  Calculate accuracy of trajectories based on reference calculations.
*/

#include "libtrac.h"
#include <gsl/gsl_sort.h>


int main(
  int argc,
  char *argv[]) {

  atm_t *atm1, *atm2;

  char *name, atmStartName[LEN];
  char *year, *mon, *day, *hour, *min;

  ctl_t ctl;

  FILE *out;

  double maxd=0, d, x0[3], x1[3], *lon1, *lat1, *p1, *dh1, *dv1, *lon2, *lat2, *p2,
    *dh2, *dv2, ahtd, avtd, ahtd2, avtd2, rhtd, rvtd, t;

  static double distv[NP];
  int distiv=0;
  
  static double disth[NP];
  int distih=0;
    
  int ip, f, i=0;

  /* Allocate... */
  ALLOC(atm1, atm_t, 1);
  ALLOC(atm2, atm_t, 1);
  ALLOC(lon1, double, NP);
  ALLOC(lat1, double, NP);
  ALLOC(p1, double, NP);
  ALLOC(dh1, double, NP);
  ALLOC(dv1, double, NP);
  ALLOC(lon2, double, NP);
  ALLOC(lat2, double, NP);
  ALLOC(p2, double, NP);
  ALLOC(dh2, double, NP);
  ALLOC(dv2, double, NP);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG
      ("Give parameters: <outfile> <atm1a> <atm1b> [<atm2a> <atm2b> ...]");

  strcpy(atmStartName, argv[2] );
  if (argv[2][strlen(argv[2]) - 1] == 'c') {
    ctl.atm_iformat = 1;
  }

  /* Write info... */
  printf("Write trajectory analysis data: %s\n", argv[1]);

  /* Create output file... */
  if (!(out = fopen(argv[1], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time [s]\n"
	  "# $2 = AHTD(t) [km]\n"
	  "# $3 = RHTD(t) [%%]\n"
	  "# $4 = AVTD(t) [km]\n"
	  "# $5 = RVTD(t) [%%]\n"
	  "# $6 = h standard deviation (t) [km]\n"
	  "# $7 = v standard deviation (t) [km]\n"
	  "# $8 = h travel distance (t) [km]\n"
	  "# $9 = v travel distance (t) [km]\n");
  fprintf(out,
	  "# $10 = 0%% AHTD(t) [km]\n"
	  "# $11 = 10%% AHTD(t) [km]\n"
	  "# $12 = 25%% AHTD(t) [km]\n"
	  "# $13 = 50%% AHTD(t) [km]\n"
	  "# $14 = 75%% AHTD(t) [km]\n"
	  "# $15 = 90%% AHTD(t) [km]\n"
	  "# $16 = 100%% AHTD(t) [km]\n"
	  "# $17 = 0%% AVTD(t) [km]\n"
	  "# $18 = 10%% AVTD(t) [km]\n"
	  "# $19 = 25%% AVTD(t) [km]\n"
	  "# $20 = 50%% AVTD(t) [km]\n"
	  "# $21 = 75%% AVTD(t) [km]\n"
	  "# $22 = 90%% AVTD(t) [km]\n"
	  "# $23 = 100%% AVTD(t) [km]\n");

  /* Loop over file pairs... */
  for (f = 2; f < argc; f += 2) {

    distiv=0;
    distih=0;
    maxd=0;
    /* Read atmopheric data... */
    read_atm(NULL, argv[f], atm1, &ctl);
    read_atm(NULL, argv[f + 1], atm2, &ctl);

    /* Check if structs match... */
    if (atm1->np != atm2->np)
      ERRMSG("Different numbers of parcels!");
    for (ip = 0; ip < atm1->np; ip++)
      if (fabs(atm1->time[ip] - atm2->time[ip]) > 1)
	printf("WARNING! Times do not match! dt=%f\n",
	       fabs(atm1->time[ip] - atm2->time[ip]));

    /* Init... */
    ahtd = avtd = rhtd = rvtd = 0;
    ahtd2 = avtd2 = 0;

    /* Loop over air parcels... */
    for (ip = 0; ip < atm1->np; ip++) {

      /* Calculate total length of trajectories... */
      if (f > 2) {
	dv1[ip] += fabs(Z(p1[ip]) - Z(atm1->p[ip]));
	geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
	geo2cart(0, lon1[ip], lat1[ip], x1);
	dh1[ip] += DIST(x0, x1);

	dv2[ip] += fabs(Z(p2[ip]) - Z(atm2->p[ip]));
	geo2cart(0, atm2->lon[ip], atm2->lat[ip], x0);
	geo2cart(0, lon2[ip], lat2[ip], x1);
	dh2[ip] += DIST(x0, x1);
      }

      /* Save last locations of air parcels... */
      lon1[ip] = atm1->lon[ip];
      lat1[ip] = atm1->lat[ip];
      p1[ip] = atm1->p[ip];
      lon2[ip] = atm2->lon[ip];
      lat2[ip] = atm2->lat[ip];
      p2[ip] = atm2->p[ip];

      /* Sum up the vertical error... */
      d=fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]));
      avtd +=  d;
      avtd2 +=  d*d;
      distv[distiv++] = d;

      /* Sum up the horizontal error... */
      geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
      geo2cart(0, atm2->lon[ip], atm2->lat[ip], x1);
      d = DIST(x0, x1);
      ahtd += d;
      ahtd2 += d*d;
      disth[distih++] = d;
      if(d>maxd){
	maxd=d;
	i=ip;
      }
      
      /* Sum up the relative transport devation... */
      if (f > 2) {
	t = ((200. * DIST(x0, x1)) / (dh1[ip] + dh2[ip]) / atm1->np);
	if (gsl_finite(t))
	  rhtd += t;
	t = ((200. * fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]))) /
	     (dv1[ip] + dv2[ip]) / atm1->np);
	if (gsl_finite(t))
	  rvtd += t;
      }
    }
    
    d=0;
/*    printf("%8.2f Max Fehler: %d  %g\n", atm1->time[0], i, maxd);*/
    for(ip=0; ip<distih; ip++) {
      if(disth[ip] > d) {
	d=disth[ip];
        i=ip;
      }
    }
/*    printf("%8.2f Max Abweichung: %d  %g\n", atm1->time[0], i, d);*/
    
    gsl_sort(distv, (size_t)1, (size_t)distiv);
    gsl_sort(disth, (size_t)1, (size_t)distih);
    
    

    /* Get date from filename... */
    for (ip = (int) strlen(argv[f]) - 1; argv[f][ip] != '/' || ip == 0; ip--);
    name = strtok(&(argv[f][ip]), "_");
    year = strtok(NULL, "_");
    mon = strtok(NULL, "_");
    day = strtok(NULL, "_");
    hour = strtok(NULL, "_");
    name = strtok(NULL, "_");
    min = strtok(name, ".");
    time2jsec(atoi(year), atoi(mon), atoi(day), atoi(hour), atoi(min), 0, 0,
	      &t);

    /* Write output... */
    d=(double) atm1->np;
    
    /*printf("%.2f %g %g %g %g %g %g %g %g\n", t, ahtd/d, rhtd, avtd/d, rvtd, sqrt(ahtd2/d - gsl_pow_2(ahtd/d)), sqrt(avtd2/d - gsl_pow_2(avtd/d)),  (dh1[ip] + dh2[ip])/2,  (dv1[ip] + dv2[ip])/2);*/
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
	    t, ahtd/d, rhtd, avtd/d, rvtd, sqrt(ahtd2/d - gsl_pow_2(ahtd/d)), sqrt(avtd2/d - gsl_pow_2(avtd/d)),  
	    (dh1[ip] + dh2[ip])/2,  (dv1[ip] + dv2[ip])/2,
	    disth[0], disth[distih/10], disth[distih/4], disth[distih/2], disth[distih-distih/4], disth[distih-distih/10], disth[distih-1],
	    distv[0], distv[distiv/10], distv[distiv/4], distv[distiv/2], distv[distiv-distiv/4], distv[distiv-distiv/10], distv[distiv-1]
 	  );
  }
  
  /* read start atm file */
  /*printf("\n\n#Biggest start-end devations\n");
  read_atm(NULL, atmStartName, atm1, &ctl);
  for (ip = 0; ip < atm1->np; ip++) {
    geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
    geo2cart(0, atm2->lon[ip], atm2->lat[ip], x1);
    disth[ip]=DIST(x0, x1);
  }
  for(distih=0; distih<10; distih++) {
    d=0;
    i=0;
    for(ip=0; ip<atm1->np; ip++) {
      if(disth[ip] > d) {
        d=disth[ip];
        i=ip;
      }
    }
    printf("%8.2f %g %g %g %u %g %8.2f %g %g %g SEDEV%s\n", atm1->time[i], Z(atm1->p[i]), atm1->lon[i],
	   atm1->lat[i], i, d, atm2->time[i], Z(atm2->p[i]), atm2->lon[i], atm2->lat[i], argv[f - 1]);
    disth[i]=0;
  }*/
  
  
  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm1);
  free(atm2);
  free(lon1);
  free(lat1);
  free(p1);
  free(dh1);
  free(dv1);
  free(lon2);
  free(lat2);
  free(p2);
  free(dh2);
  free(dv2);

  return EXIT_SUCCESS;
}
