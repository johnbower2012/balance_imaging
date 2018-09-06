//lagrange.h
//Lagrange Interpolation
#ifndef	LAGRANGE_H
#define	LAGRANGE_H

#include <cmath>

//Vector
double	*	vector(short	nl,
				   short	nh)
{
	double	*v;
	v = (double	*)malloc((size_t)	((nh - nl + 1 + 1) * sizeof(double)));
	if(!v)
	{
		printf("allocation failure in vector()\n");
	}
	return	v-nl+1;
}

void	free_vector(double	*	v,
					short	nl,
					short	nh)
{
	free((char	*)(v + nl - 1));
}


// Return the value of the function by *y.
// *dy returns skeptically the error of the value of the function.
void	polint(double	xa[],
			   double	ya[],
			   short	n,
			   double	x,
			   double	*	y,
			   double	*	dy)
{
	double	tempy = 0.0;
	double	tempdy = 0.0;
	double	*	l;
	l = vector(1,n);

	double p;
	short k, j;
	p = 0.0;
	// Guarantee single value function.
	for(k=1;k<=n;k++)
	{
		for(j=1;j<=n&&j!=k;j++)
		{
			if(xa[k] == xa[j])
			{
				printf(" The function is not a single value one!\n");
				exit(0);
			}
		}
	}
	for(k=1;k<=n;k++)
	{
		l[k]=1.0;
		for(j=1;j<=n;j++)
		{
			if(j != k)
			{
				l[k] *= (x - xa[j]) / (xa[k] - xa[j]);
			}
		}
		p += l[k] * ya[k];
	}
	tempy = p;
	* dy = tempdy;
	* y = tempy;
	free_vector(l, 1, n);
}
#endif