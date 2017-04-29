/*! 
  \file
  Calculate mean longitude, latitude, height and standard deviation.
*/

#include "libtrac.h"
#define GROUPS 22

int main(
  int argc,
  char *argv[]) {

  atm_t *atm1;

  char filename[LEN];
  char *name;
  char *year, *mon, *day, *hour, *min;

  ctl_t ctl;

  FILE *out;

  double tmp, t0, dt=60, t, h[GROUPS], h2[GROUPS], pos[3], meanpos[3];
  double lon[GROUPS], lat[GROUPS], np;
  double dist[GROUPS],dist2[GROUPS], hdist[GROUPS];

  
  double mean_lon=0;
  double mean_lat=0;
  double mean_h=0;
  double mean_h2=0; 
  double mean_hdist=0; 
  double mean_dist=0;
  double mean_dist2=0;
  
  int ip, f, g;
  size_t snp = 0;

  /* Allocate... */
  ALLOC(atm1, atm_t, 1);


  /* Check arguments... */
  if (argc < 3)
    ERRMSG
      ("Give parameters: <outfile> <atm1> {atm2} ");

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
	  "# $2 = mean height [km]\n"
	  "# $3 = mean vertical distance [km]\n"
	  "# $4 = mean vertical standard deviation [km]\n"
          "# $5 = mean lon [deg]\n"
          "# $6 = mean lat [deg]\n"
          "# $7 = mean horizontal distance [km]\n"
          "# $8 = mean horizontal standard deviation [km]\n"
          "# $9 = number of parcels\n\n");
   
  /* Loop over file pairs... */
  for (f = 2; f < argc; f += 1) {

    /* Get date from filename... */
    strcpy(filename, argv[f]);

    for (ip = (int) strlen(argv[f]) - 1; argv[f][ip] != '/' || ip == 0; ip--);
    name = strtok(&(argv[f][ip]), "_");
    year = strtok(NULL, "_");
    mon = strtok(NULL, "_");
    day = strtok(NULL, "_");
    hour = strtok(NULL, "_");
    name = strtok(NULL, "_");
    min = strtok(name, ".");
    time2jsec(atoi(year), atoi(mon), atoi(day), atoi(hour), atoi(min), 0, 0, &t);
    t0=t;

    /* Read atmopheric data... */
    read_atm(NULL, filename, atm1, &ctl);
    
    /* Init... */
    for(g=0;g<GROUPS; g++) {
      lon[g] = lat[g] = h[g] = h2[g] = 0;
      dist[g]=0;
      dist2[g]=0;
      hdist[g]=0;
    }
   mean_lon=0;
   mean_lat=0;
   mean_h=0;
   mean_h2=0;  
   mean_hdist=0;
   mean_dist=0;
   mean_dist2=0;
    /* Loop over air parcels... */
    for (ip = 0; ip < atm1->np; ip++) {
       g=ip/2000;
      if(t0 <= 0 || fabs(atm1->time[ip]-t0) <= dt) {
        geo2cart(Z(atm1->p[ip]), atm1->lon[ip], atm1->lat[ip], pos);
        lon[g] += atm1->lon[ip];
        lat[g] += atm1->lat[ip];
        h[g] += Z(atm1->p[ip]);
        h2[g] += gsl_pow_2(Z(atm1->p[ip]));
      }
    }
    snp = (size_t)atm1->np;
    np=(double)snp/GROUPS;
    
    for (ip = 0; ip < atm1->np; ip++) {
       g=ip/2000;
       geo2cart(h[g]/np, lon[g]/np, lat[g]/np, meanpos);
       geo2cart(h[g]/np, atm1->lon[ip], atm1->lat[ip], pos);
       
       hdist[g]+=fabs(h[g]/np - Z(atm1->p[ip]));
       dist[g]+=DIST(meanpos, pos);
       dist2[g]+=gsl_pow_2(DIST(meanpos, pos));
    }
 
    for(g=0;g<GROUPS; g++) {
      mean_lon+=lon[g];
      mean_lat+=lat[g];
      mean_h+=h[g];
      tmp=sqrt(h2[g]/np - gsl_pow_2(h[g]/np));
      if(gsl_isnan(tmp) || !gsl_finite(tmp)) tmp=0;
      mean_h2+=tmp;
      mean_hdist+=hdist[g];
      mean_dist+=dist[g];
      tmp=sqrt(dist2[g]/np - gsl_pow_2(dist[g]/np));
      if(gsl_isnan(tmp) || !gsl_finite(tmp)) tmp=0;
      mean_dist2+=tmp;
      
    fprintf(out,"# %u %.2f %g %g %g %g %g %g %g %ld\n", g, t, h[g]/np, hdist[g]/np, sqrt(h2[g]/np - gsl_pow_2(h[g]/np)), 
	   lon[g]/np, lat[g]/np, dist[g]/np, sqrt(dist2[g]/np - gsl_pow_2(dist[g]/np)), (long int)snp);
    }
    mean_lon=mean_lon/np/GROUPS;
    mean_lat=mean_lat/np/GROUPS;
    mean_h=mean_h/np/GROUPS;
    mean_h2=mean_h2/GROUPS;
    mean_hdist=mean_hdist/np/GROUPS;
    mean_dist=mean_dist/np/GROUPS;
    mean_dist2=mean_dist2/GROUPS;
    /* Write output... */
    fprintf(out, "%.2f %g %.3f %.3f %g %g %.1f %.1f %ld\n", t, mean_h, mean_hdist, mean_h2, mean_lon, mean_lat, mean_dist, mean_dist2, (long int)snp);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(atm1);

  return EXIT_SUCCESS;
}
