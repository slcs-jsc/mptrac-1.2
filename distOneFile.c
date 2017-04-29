/*! 
  \file
  Calculate accuracy of trajectories based on reference calculations.
*/

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  atm_t *atm1, *atm2;

  ctl_t ctl;

  FILE *out;

  double d, x0[3], x1[3], lon1=0, lat1=0, p1=0, dh1=0, dv1=0, lon2=0, lat2=0, p2=0,
    dh2=0, dv2=0, ahtd=0, avtd=0, ahtd2=0, avtd2=0, rhtd=0, rvtd=0, t=0;

  int ip;

  /* Allocate... */
  ALLOC(atm1, atm_t, 1);
  ALLOC(atm2, atm_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG
      ("Give parameters: <outfile> <atm1a> <atm1b> [<atm2a> <atm2b> ...]");

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
	  "# $2 = lon0 [deg]\n"
	  "# $3 = lat0 [deg]\n"
	  "# $4 = lon1 [deg]\n"
	  "# $5 = lat1 [deg]\n"
	  "# $6 = AHTD(t) [km]\n"
	  "# $7 = RHTD(t) [%%]\n"
	  "# $8 = AVTD(t) [km]\n"
	  "# $9 = RVTD(t) [%%]\n"
	  "# $10 = h standard deviation (t) [km]\n"
	  "# $11 = v standard deviation (t) [km]\n"
	  "# $12 = h travel distance (t) [km]\n"
	  "# $13 = v travel distance (t) [km]\n\n");
  read_atm(NULL, argv[2], atm1, &ctl);
  read_atm(NULL, argv[3], atm2, &ctl);

    /* Check if structs match... */
    if (atm1->np != atm2->np)
      ERRMSG("Different numbers of parcels!");
    for (ip = 0; ip < atm1->np; ip++)

    /* Loop over timesteps... */
    for (ip = 0; ip < atm1->np; ip++) {
      ahtd = avtd = rhtd = rvtd = 0;
      ahtd2 = avtd2 = 0;

      /* Calculate total length of trajectories... */
      if (ip > 0) {
	dv1 += fabs(Z(p1) - Z(atm1->p[ip]));
	geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
	geo2cart(0, lon1, lat1, x1);
	dh1 += DIST(x0, x1);

	dv2 += fabs(Z(p2) - Z(atm2->p[ip]));
	geo2cart(0, atm2->lon[ip], atm2->lat[ip], x0);
	geo2cart(0, lon2, lat2, x1);
	dh2 += DIST(x0, x1);
      }

      /* Save last locations of air parcels... */
      lon1 = atm1->lon[ip];
      lat1 = atm1->lat[ip];
      p1 = atm1->p[ip];
      lon2 = atm2->lon[ip];
      lat2 = atm2->lat[ip];
      p2 = atm2->p[ip];

      /* Sum up the vertical error... */
      d=fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]));
      avtd +=  d;
      avtd2 +=  d*d;

      /* Sum up the horizontal error... */
      geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
      geo2cart(0, atm2->lon[ip], atm2->lat[ip], x1);
      d = DIST(x0, x1);
      ahtd += d;
      ahtd2 += d*d;

      /* Sum up the relative transport devation... */
      if (ip > 0) {
	t = ((200. * DIST(x0, x1)) / (dh1 + dh2) / atm1->np);
	if (gsl_finite(t))
	  rhtd += t;
	t = ((200. * fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]))) /
	     (dv1 + dv2) / atm1->np);
	if (gsl_finite(t))
	  rvtd += t;
      }
      
   printf("%.2f (%g,%g)-(%g,%g) %g %g %g %g %g %g %g %g\n", atm1->time[ip],atm1->lon[ip],atm1->lat[ip],  atm2->lon[ip], atm2->lat[ip], ahtd, rhtd, avtd, rvtd,
      sqrt(ahtd2 - gsl_pow_2(ahtd)), sqrt(avtd2 - gsl_pow_2(avtd)),  (dh1 + dh2)/2,  (dv1 + dv2)/2);
   fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %g %g %g\n", atm1->time[ip],atm1->lon[ip],atm1->lat[ip],  atm2->lon[ip], atm2->lat[ip], ahtd, rhtd, avtd, rvtd,
      sqrt(ahtd2 - gsl_pow_2(ahtd)), sqrt(avtd2 - gsl_pow_2(avtd)),  (dh1 + dh2)/2,  (dv1+ dv2)/2);

    }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm1);
  free(atm2);

  return EXIT_SUCCESS;
}
