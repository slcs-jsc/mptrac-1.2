/*! 
  \file
  Calculate accuracy of trajectories based on reference calculations.
*/

#include "libtrac.h"
#include <gsl/gsl_sort.h>


#define MAX_CLASSES 1000


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

  double minlon=-190, maxlon=190, minlat=-95, maxlat=95, minz=0, maxz=120;
  
  double sumDistH, sumDistV;
  static double distv[NP];
  int distiv=0;

  double dh, dv, pdfh[MAX_CLASSES], pdfv[MAX_CLASSES], cdfv[MAX_CLASSES], cdfh[MAX_CLASSES]; 
  int class, classes;
  size_t * hrank, *vrank;  
  static int distid[NP];
  static double disth[NP];
  int distih=0;
  int *outOfBox;
  int noutOfBox=0;
  
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
  ALLOC(outOfBox, int, NP);

  /* Check arguments... */
  if (argc < 10)
    ERRMSG
      ("Give parameters: <min lon> <max lon> <min lat> <max lat> <min z> <max z> <outfile> <atm1a> <atm1b> [<atm2a> <atm2b> ...]");

  strcpy(atmStartName, argv[8] );
  if (argv[8][strlen(argv[8]) - 1] == 'c') {
    ctl.atm_iformat = 1;
  }
  minlon=atof(argv[1]);
  maxlon=atof(argv[2]);
  minlat=atof(argv[3]);
  maxlat=atof(argv[4]);
  minz=atof(argv[5]);
  maxz=atof(argv[6]);

  /* Write info... */
  printf("Transport deviations for parcels in lon [%4.2f:%4.2f] lat [%4.2f:%4.2f] z [%3.2f:%3.2f]\n", minlon, maxlon, minlat, maxlat, minz, maxz);
  printf("Write trajectory analysis data: %s\n", argv[7]);

  /* Create output file... */
  if (!(out = fopen(argv[7], "w")))
    ERRMSG("Cannot create file!");


  /* Write header... */
  fprintf(out, "# $1 = parcel number\n"
               "# $2 = cumulative probability\n"
               "# $3 = horizontal deviation\n"
               "# $4 = horizontal parcel index\n"
               "# $5 = vertical deviation\n"
               "# $6 = vertical parcel index\n");

  /* Loop over file pairs... */
  for (f = 8; f < argc; f += 2) {

    distiv=0;
    distih=0;
    maxd=0;
    sumDistH=0;
    sumDistV=0;

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
      
      /* Skip unvalid parcels */
      if(outOfBox[ip] > 0) {
	continue;
      }
      
      if((minlon>atm1->lon[ip]) || (maxlon<atm1->lon[ip]) || 
         (minlat>atm1->lat[ip]) || (maxlat<atm1->lat[ip]) || 
         (minz>Z(atm1->p[ip]))  || (maxz<Z(atm1->p[ip]))) {
	   
	  noutOfBox++;
	  outOfBox[ip]=1;
	}

      /* Calculate total length of trajectories... */
      if (f > 8) {
	dv1[ip] += fabs(Z(p1[ip]) - Z(atm1->p[ip]));
	geo2cart(0, atm1->lon[ip], atm1->lat[ip], x0);
	geo2cart(0, lon1[ip], lat1[ip], x1);
	dh1[ip] += DIST(x0, x1);
        
	dv2[ip] += fabs(Z(p2[ip]) - Z(atm2->p[ip]));
	geo2cart(0, atm2->lon[ip], atm2->lat[ip], x0);
	geo2cart(0, lon2[ip], lat2[ip], x1);
	dh2[ip] += DIST(x0, x1);
      }
      sumDistH += (dh1[ip]+dh2[ip])/2;
      sumDistV += (dv1[ip]+dv2[ip])/2;

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

      distid[distih-1]=ip;

      if(d>maxd){
	maxd=d;
	i=ip;
      }
      
      /* Sum up the relative transport devation... */
      if (f > 8) {
	t = ((200. * DIST(x0, x1)) / (dh1[ip] + dh2[ip]));
	if (gsl_finite(t))
	  rhtd += t;
	t = ((200. * fabs(Z(atm1->p[ip]) - Z(atm2->p[ip]))) /
	     (dv1[ip] + dv2[ip]));
	if (gsl_finite(t))
	  rvtd += t;
      }
    }  
  }
  d = (double) distih;
  fprintf(out, "# parcels=%i\n", distih);
  fprintf(out, "# AHTD=%g\n", ahtd/d);
  fprintf(out, "# SDAHTD=%g\n", sqrt(ahtd2/d - gsl_pow_2(ahtd/d)));
  fprintf(out, "# AVTD=%g\n", avtd/d);
  fprintf(out, "# SDAVTD=%g\n", sqrt(avtd2/d - gsl_pow_2(avtd/d)));

  ALLOC(hrank, size_t, distih);
  ALLOC(vrank, size_t, distiv);
  gsl_sort_index(hrank, disth, (size_t)1, (size_t)distih);
  gsl_sort_index(vrank, distv, (size_t)1, (size_t)distiv);
  gsl_sort(disth, (size_t)1, (size_t)distih);
  gsl_sort(distv, (size_t)1, (size_t)distiv);

  fprintf(out, "# hmedian=%g\n", disth[distih/2]);
  fprintf(out, "# vmedian=%g\n\n", distv[distiv/2]);

  fprintf(out, "%d %g %g %d %g %d\n", 0, 0.0, 0.0, -1, 0.0, -1);

  d=(double)distih;
  for(ip=0; ip<distih; ip++) {
    fprintf(out, "%d %g %g %d %g %d\n",
         (1+ip), (1+ip)/d, disth[ip], distid[hrank[ip]], distv[ip], distid[vrank[ip]] );
  }
  
 /*
  for (ip=0; ip<MAX_CLASSES; ip++) {
    pdfh[ip]=pdfv[ip]=cdfv[ip]=cdfh[ip]=0;
  }
  
  classes = (distih);
  if(classes>MAX_CLASSES) classes=MAX_CLASSES;

  dh = disth[distih-1] / classes;
  dv = distv[distiv-1] / classes;

  for (ip = 0; ip < distih; ip++) {
    class = (int)(disth[ip] / dh);
    pdfh[class]++;        

    class = (int)(distv[ip] / dv);
    pdfv[class]++;
  }

  cdfh[0]=pdfh[0];
  cdfv[0]=pdfv[0];
  for (ip=0; ip<classes-1; ip++) {
    cdfh[ip+1]=cdfh[ip]+pdfh[ip+1];
    cdfv[ip+1]=cdfv[ip]+pdfv[ip+1];
  }
  
  fprintf(out, "# classes=%i\n", classes);
  fprintf(out, "# dh=%g\n", dh);
  fprintf(out, "# dv=%g\n\n", dv);

  d = (double) distih;
  fprintf(out, "%g %g %g %g %g %g\n", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  for (ip=0; ip<classes-1; ip++) {
     fprintf(out, "%g %g %g %g %g %g\n",
        dh * (ip + 0.5), pdfh[ip] / d, cdfh[ip] / d, 
        dv * (ip + 0.5), pdfv[ip] / d, cdfv[ip] / d);
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
