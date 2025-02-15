/*! 
 *   \file
 *     Estimate horizontal diffusivity based on Smagorinsky theory.
 *     */

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  met_t *met;

  FILE *out;

  static double dz, dzmin = 1e10, z, t, s, ls2, k[EX][EY], c = 0.1;

  static int ip, ip2, ix, iy;

  /* Allocate... */
  ALLOC(met, met_t, 1);

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <map.tab> <met>");

  /* Read control parameters... */
  read_ctl(NULL, argv[1], argc, argv, &ctl);
  z = scan_ctl(NULL, argv[1], argc, argv, "Z", -1, "", NULL);

  /* Read meteorological data... */
  read_met(argv[3], met);

  /* Find nearest pressure level... */
  for (ip2 = 0; ip2 < met->np; ip2++) {
    dz = fabs(Z(met->p[ip2]) - z);
    if (dz < dzmin) {
      dzmin = dz;
      ip = ip2;
    }
  }

  /* Calculate horizontal diffusion coefficients... */
  for (ix = 1; ix < met->nx - 1; ix++)
    for (iy = 1; iy < met->ny - 1; iy++) {
      t = (met->u[ix + 1][iy][ip] - met->u[ix - 1][iy][ip])
	/ (1000. * deg2dx(met->lon[ix + 1] - met->lon[ix - 1], met->lat[iy]))
	- (met->v[ix][iy + 1][ip] - met->v[ix][iy - 1][ip])
	/ (1000. * deg2dy(met->lat[iy + 1] - met->lat[iy - 1]));
      s = (met->u[ix][iy + 1][ip] - met->u[ix][iy - 1][ip])
	/ (1000. * deg2dy(met->lat[iy + 1] - met->lat[iy - 1]))
	+ (met->v[ix + 1][iy][ip] - met->v[ix - 1][iy][ip])
	/ (1000. * deg2dx(met->lon[ix + 1] - met->lon[ix - 1], met->lat[iy]));
      ls2 =
	c / (1. /
	     gsl_pow_2(1000. *
		       deg2dx(met->lon[ix + 1] - met->lon[ix - 1],
			      met->lat[iy]))
	     +
	     1. / gsl_pow_2(1000. *
			    deg2dy(met->lat[iy + 1] - met->lat[iy - 1])));
      k[ix][iy] = ls2 * sqrt(gsl_pow_2(t) + gsl_pow_2(s));
    }

  /* Create output file... */
  printf("Write data file: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = longitude [deg]\n"
	  "# $2 = latitude [deg]\n"
	  "# $3 = zonal wind [m/s]\n"
	  "# $4 = meridional wind [m/s]\n"
	  "# $5 = horizontal diffusivity [m^2/s]\n");

  /* Write data... */
  for (iy = 0; iy < met->ny; iy++) {
    fprintf(out, "\n");
    for (ix = 0; ix < met->nx; ix++)
      if (met->lon[ix] >= 180)
	fprintf(out, "%g %g %g %g %g\n",
		met->lon[ix] - 360.0, met->lat[iy],
		met->u[ix][iy][ip], met->v[ix][iy][ip], k[ix][iy]);
    for (ix = 0; ix < met->nx; ix++)
      if (met->lon[ix] <= 180)
	fprintf(out, "%g %g %g %g %g\n",
		met->lon[ix], met->lat[iy],
		met->u[ix][iy][ip], met->v[ix][iy][ip], k[ix][iy]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(met);

  return EXIT_SUCCESS;
}

