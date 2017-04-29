


/*! 
 *   \file
 *     Convert atmospheric data files between ASCII and netCDF.
 *     */

#include "libtrac.h"

int main(
  int argc,
  char *argv[]) {
  int i, ip;

  static atm_t atm;
  static atm_t full;
  static ctl_t ctl;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> -f {<atm_in>} <atm_out> <grid_out>");

  /* Read control parameter */
  read_ctl(NULL, argv[1], argc, argv, &ctl);

  i=6;

  for(; i<argc-2; i++) {
    /* Read atmospheric data... */
    read_atm(NULL, argv[i], &atm, &ctl);
    
    for(ip=0; ip<atm.np; ip++) {
      full.time[full.np+ip] = atm.time[ip];
      full.lon[full.np+ip] = atm.lon[ip];
      full.lat[full.np+ip] = atm.lat[ip];
      full.p[full.np+ip] = atm.p[ip];
    }    

    full.np+=atm.np;
    
  }
  /* Write atmospheric data... */
  if (argv[3][0] != '-')
    write_atm(NULL, argv[i], &full, &ctl);

  /* Write gridded data... */
  if (argv[4][0] != '-')
    write_grid(NULL, argv[i+1], &full, &ctl, 0, 1e100);

  return EXIT_SUCCESS;
}

