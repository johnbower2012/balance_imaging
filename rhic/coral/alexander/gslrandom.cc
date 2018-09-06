using namespace std;

CRandom::CRandom(int seed){
  //Choose the algorithm, see gsldoc.ps (page 178)
  //randy=gsl_rng_alloc(gsl_rng_taus); // crappy but fast
  //randy=gsl_rng_alloc(gsl_rng_ranlxd2); // for anal-retentive types
  //randy=gsl_rng_alloc(gsl_rng_ran2); // Not really double precision numbers
  //randy=gsl_rng_alloc(gsl_rng_knuthran2); // sort of slow
  randy=gsl_rng_alloc(gsl_rng_ranlxd1); // Just right

  gsl_rng_set(randy,seed);
}

void CRandom::reset(int seed){
  gsl_rng_set(randy,seed);
}

double CRandom::ran(void){
  return gsl_rng_uniform(randy);
}

long unsigned int CRandom::iran(unsigned long int imax){
  return gsl_rng_uniform_int(randy,imax);
}

double CRandom::gauss(void){
  return gsl_ran_ugaussian(randy);
}

void CRandom::gauss2(double *randy1,double *randy2){
  double x,y,r2,r,c,s;
TRY_AGAIN:
  x=1.0-2.0*gsl_rng_uniform(randy);
  y=1.0-2.0*gsl_rng_uniform(randy);
  r2=x*x+y*y;
  if(r2>1.0) goto TRY_AGAIN;
  r=sqrt(r2);
  c=x/r;
  s=y/r;
  *randy1=c*sqrt(-2.0*log(r2));
  *randy2=(s/c)**randy1;
}

void CRandom::thermal(double mass,double T,double *p){
  double pmag,cthet,sthet,phi;
  double e,rootmt,g1,g2;
  const double PI=4.0*atan(1.0);
  if(mass/T<10.0){
  RETRY1:
    pmag=-T*log(ran()*ran()*ran());
    e=sqrt(mass*mass+pmag*pmag);
    if(ran()>exp(-(e-pmag)/T)) goto RETRY1;
    phi=2.0*PI*ran();
    cthet=1.0-2.0*ran();
    sthet=sqrt(1.0-cthet*cthet);
    p[0]=e;
    p[3]=pmag*cthet;
    p[1]=pmag*sthet*cos(phi);
    p[2]=pmag*sthet*sin(phi);
  }
  else{
    rootmt=sqrt(mass*T);
    gauss2(&g1,&g2);
    p[1]=rootmt*g1;p[2]=rootmt*g2;
    gauss2(&g1,&g2);
    p[3]=rootmt*g1;
    p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+mass*mass);
  }
}

