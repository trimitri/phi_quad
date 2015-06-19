#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <limits.h>

/*define the parameters globally */
double kappa;
int main() {
  int l, ii, i0, i1, imax;
  int iniran;
  double coeff[30];
  /* for timing */
  double qq, chi;
  FILE *fp1;

  fp1 = fopen("series_chi", "r");

  for (i0 = 0; i0 < 21; i0++) {
    fscanf(fp1, "%d %lf", &ii, &coeff[i0]);
  }
  printf("give kappa \n");
  scanf("%lf", &kappa);

  qq = 1.;
  for (i0 = 0; i0 < 21; i0++) {
    chi += qq * coeff[i0];
    if (i0 == 20) printf("%lf \n", qq * coeff[i0]);
    qq = qq * (2. * kappa);
  }
  printf("%lf \n", chi);
}
