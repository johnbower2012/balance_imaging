//interpolation.h
// Interpolation of Neville method.
// Improement of Lagrange method.
// Starting with the most vicinity points and taking into consideration the correction from farther points.
// Building up Neville recurrence tablet. Of order O(N^2).
// The last time correction can be regarded as an upper bound of the error. This is not always true.

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <cmath>

//Vector
double * Vector(short nl, short nh)
{
  double *v;
  v = (double *)malloc((size_t) ((nh - nl + 1 + 1) * sizeof(double)));
  if(!v)
    {
      printf("allocation failure in Vector()\n");
    }
  return v-nl+1;
}

void free_Vector(double * v, short nl, short nh)
{
  free((char *)(v + nl - 1));
}

// Return the value of the function by *y.
// *dy returns skeptically the error of the value of the function.
void polint(double xa[], double ya[], short n,
 double x, double * y, double * dy)
{
  short i, m;
  short ns = 1;
  double den, dif, dift, ho, hp, w;
  double * c;
  double * d;

  double tempy = 0.0;
  double tempdy = 0.0;

  c = Vector(1, n);
  d = Vector(1, n);

  dif = fabs(x - xa[1]);
  for(i=1;i<=n;i++)
    {
      if((dift = fabs(x - xa[i])) < dif)
	{
	  ns = i; //ns keeps the nearest point
	  dif = dift;
	} // if this point is nearer to x, calculate it's c[] and d[].
      c[i] = ya[i];
      d[i] = ya[i];
    } //for(i=1;i<=n;i++) //the nearest point as entrance of construction of c[] and d[].
  tempy = ya[ns--];//The initial value y, the nearest left point of y.

  //Actually this doesn't find the NEAREST point. However this won't affect the accuracy of the whole algorithm.
  for(m=1;m<n;m++) //Column of Neville tablet
    {
      for(i=1;i<=n-m;i++)
	{
	  ho = xa[i] - x;
	  hp = xa[i + m] - x;
	  w = c[i + 1] - d[i];
	  if(fabs(den = ho - hp) <= 1e-12)
	    {
	      printf("Error in routine polint!\n");
	    } // Guarantee single value function.
	  den = w / den;
	  d[i] = hp * den;
	  c[i] = ho * den;
	} //for(i=1;i<=n-m;i++)
      tempdy = (2 * ns < (n - m)) ? c[ns +1 ] : d[ns--];
      tempy += tempdy;
    }
  * dy = tempdy;
  * y = tempy;
  free_Vector(d, 1, n);
  free_Vector(c, 1, n);
}
#endif //INTERPOLATION_H
