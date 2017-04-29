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
    *dh2, *dv2, ahtd, axtd, aytd, ytd, xtd, vtd, avtd, ahtd2, axtd2, aytd2, avtd2, rhtd, rvtd, t;

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
	  "# $3 = AVTD(t) [km]\n"
	  "# $4 = AXTD(t) [km]\n"
	  "# $5 = AYTD(t) [km]\n"
	  "# $6 = XTD(t) [km]\n"
	  "# $7 = YTD(t) [km]\n"
	  "# $8 = x standard deviation (t) [km]\n"
	  "# $9 = y standard deviation (t) [km]\n"
	  "# $10 = v standard deviation (t) [km]\n"
	  "# $11 = h standard deviation (t) [km]\n");

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
    xtd = ytd = vtd = ahtd = aytd = axtd = avtd = rhtd = rvtd = 0;
    ahtd2 = aytd2 = axtd2 = avtd2 = 0;

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
      d=Z(atm1->p[ip]) - Z(atm2->p[ip]);
      vtd +=  d;
      avtd +=  fabs(d);
      avtd2 +=  d*d;
      distv[distiv++] = d;
      
      /* Sum up the horizontal-x error... */
      geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
      geo2cart(0, atm2->lon[ip], atm2->lat[ip], x1);
      d = DIST(x0, x1);
      ahtd+= d;
      ahtd2 += d*d;

      /* Sum up the horizontal-x error... */
      geo2cart(0, atm1->lon[ip], (atm1->lat[ip]+atm2->lat[ip])/2, x0);
      geo2cart(0, atm2->lon[ip], (atm1->lat[ip]+atm2->lat[ip])/2, x1);
      d = DIST(x0, x1);
      xtd += deg2dx(atm1->lon[ip]-atm2->lon[ip],  (atm1->lat[ip]+atm2->lat[ip])/2 );
      axtd += d;
      axtd2 += d*d;
            
      /* Sum up the horizontal-y error... */
      geo2cart(0, 0, atm1->lat[ip], x0);
      geo2cart(0, 0, atm2->lat[ip], x1);
      d = DIST(x0, x1);
      ytd += deg2dy(atm1->lat[ip]-atm2->lat[ip] );
      aytd += d;
      aytd2 += d*d;
      
      disth[distih++] = d;
      if(d>maxd){
	maxd=d;
	i=ip;
      }
    }
    
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
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g\n", 
	    t, ahtd/d, avtd/d, axtd/d, aytd/d,
	    xtd/d, ytd/d, sqrt(axtd2/d - gsl_pow_2(xtd/d)), sqrt(aytd2/d - gsl_pow_2(ytd/d)), sqrt(avtd2/d - gsl_pow_2(vtd/d)), sqrt(ahtd2/d - gsl_pow_2(ahtd/d)));
  }
  
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
