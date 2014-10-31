#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define q 1.602176565E-19   //Carga elemental
#define m 1.67262178E-27    //Masa del Proton
#define c 3E8               //Velocidad de la Luz
#define Bo 3.07E-5          //Campo Magnetico base
#define Rt 6378100        //Radio de la Tierra
#define h  0.00001

void rungeKutta(double t_old, double *t_new, double x_old, double y_old, double z_old, double *xnew, double *ynew, double *znew, double vx_old, double vy_old, double vz_old, double *vxnew, double *vynew, double *vznew, double *gamma);

int main(int argc, char **par){

  double alpha;       //Pitch Angle
  double Eko;   //Energia Cinetica Inicial
  Eko = atoi(par[1]);
  alpha = atoi(par[2]);
  Eko=Eko*(1.60210E-13)/1;
  alpha = alpha*(M_PI/180);

  double T = 10.0;         //Tiempo total de la simulacion
  int npoints = (int)(T/h);     //numero de pasos

  //Vectores
  double *t  = malloc(npoints*sizeof(double));
  double *x  = malloc(npoints*sizeof(double));
  double *y = malloc(npoints*sizeof(double));
  double *z = malloc(npoints*sizeof(double));

  double *vx = malloc(npoints*sizeof(double));
  double *vy = malloc(npoints*sizeof(double));
  double *vz = malloc(npoints*sizeof(double));

  double *tnew = malloc(1*sizeof(double));
  double *xnew = malloc(1*sizeof(double));
  double *ynew = malloc(1*sizeof(double));
  double *znew = malloc(1*sizeof(double));

  double *vxnew = malloc(1*sizeof(double));
  double *vynew = malloc(1*sizeof(double));
  double *vznew = malloc(1*sizeof(double));
  //Condicionles Iniciales
  x[0] = 2*Rt;
  y[0] = 0;
  z[0] = 0;

  //Calcular Velocidad
  double *v = malloc(1*sizeof(double));
  *v = c*sqrt(1-(1/(pow((Eko/(m*c*c)+1),2))));
  vx[0]=0;
  vy[0]= (*v)*(sin(alpha));
  vz[0]= (*v)*cos(alpha);
  //Gamma
  double *gamma = malloc(1*sizeof(double));
  *gamma = 1/(sqrt(1-(pow(*v,2))/(pow(c,2))));
  //Resolver Ecuacion

   FILE *fileout;
   fileout = fopen("trayectoria_E_alpha.dat", "w");
   int i =0;
   fprintf(fileout, "%f\t", vx[0]/Rt);
   fprintf(fileout, "%f\t", vy[0]/Rt);
   fprintf(fileout, "%f\n", vz[0]/Rt);
   for(i=1;i<npoints;i++){
     rungeKutta(t[i-1],tnew, x[i-1],y[i-1],z[i-1],xnew,ynew,znew, vx[i-1],vy[i-1],vz[i-1],vxnew,vynew,vznew,gamma);
     x[i]=*xnew;
     y[i]=*ynew;
     z[i]=*znew;
     vx[i]=*vxnew;
     vy[i]=*vynew;
     vz[i]=*vznew;
     fprintf(fileout, "%f\t", vx[i]/Rt);
     fprintf(fileout, "%f\t", vy[i]/Rt);
     fprintf(fileout, "%f\n", vz[i]/Rt);
   }
   fclose(fileout);

  return 0;

}

double primex(double *gamma, double x, double y, double z, double vx, double vy, double vz){
  double r = sqrt(x*x+y*y+z*z);
  double r5 = pow(r,5);
  double by = (-Bo*pow(Rt,3)/r5)*(3*y*z);
  double bz = (-Bo*pow(Rt,3)/r5)*(2*z*z-x*x-y*y);
  double xp = (q/((*gamma)*m))*(vy*bz-by*vz);
  return xp;
}
double primey(double *gamma, double x, double y, double z, double vx, double vy, double vz){
  double r = sqrt(x*x+y*y+z*z);
  double r5 = pow(r,5);
  double bx = (-Bo*pow(Rt,3)/r5)*(3*x*z);
  double bz = (-Bo*pow(Rt,3)/r5)*(2*z*z-x*x-y*y);
  double yp = -(q/((*gamma)*m))*(vx*bz-bx*vz);
  return yp;
}
double primez(double *gamma, double x, double y, double z, double vx, double vy, double vz){
  double r = sqrt(x*x+y*y+z*z);
  double r5 = pow(r,5);
  double bx = (-Bo*pow(Rt,3)/r5)*(3*x*z);
  double by = (-Bo*pow(Rt,3)/r5)*(3*y*z);
  double zp = (q/((*gamma)*m))*(vx*by-bx*vy);
  return zp;
}


void rungeKutta(double t_old, double *t_new, double x_old, double y_old, double z_old, double *xnew, double *ynew, double *znew, double vx_old, double vy_old, double vz_old, double *vxnew, double *vynew, double *vznew, double *gamma){

  double k_1_primex = primex(gamma,x_old,y_old,z_old,vx_old,vy_old,vz_old);
  double k_1_primey = primey(gamma,x_old,y_old,z_old,vx_old,vy_old,vz_old);
  double k_1_primez = primez(gamma,x_old,y_old,z_old,vx_old,vy_old,vz_old);

  //first step
  double t1 = t_old+ (h/2.0);

  double vx1 = vx_old + (h/2.0) * k_1_primex;
  double vy1 = vy_old + (h/2.0) * k_1_primey;
  double vz1 = vz_old + (h/2.0) * k_1_primez;
  double x1 = x_old + (h/2.0) * vx1;
  double y1 = y_old + (h/2.0) * vy1;
  double z1 = z_old + (h/2.0) * vz1;

  double k_2_primex = primex(gamma,x1,y1,z1,vx1,vy1,vz1);
  double k_2_primey = primey(gamma,x1,y1,z1,vx1,vy1,vz1);
  double k_2_primez = primez(gamma,x1,y1,z1,vx1,vy1,vz1);
  //second step

  double t2 = t_old + (h/2.0);

  double vx2 = vx_old + (h/2.0) * k_2_primex;
  double vy2 = vy_old + (h/2.0) * k_2_primey;
  double vz2 = vz_old + (h/2.0) * k_2_primez;
  double x2 = x_old + (h/2.0) * vx2;
  double y2 = y_old + (h/2.0) * vy2;
  double z2 = z_old + (h/2.0) * vz2;

  double k_3_primex = primex(gamma,x2,y2,z2,vx2,vy2,vz2);
  double k_3_primey = primey(gamma,x2,y2,z2,vx2,vy2,vz2);
  double k_3_primez = primez(gamma,x2,y2,z2,vx2,vy2,vz2);

   //third step
  double t3 = t_old + h;
  double vx3 = vx_old + (h/2.0) * k_3_primex;
  double vy3 = vy_old + (h/2.0) * k_3_primey;
  double vz3 = vz_old + (h/2.0) * k_3_primez;
  double x3 = x_old + (h/2.0) * vx3;
  double y3 = y_old + (h/2.0) * vy3;
  double z3 = z_old + (h/2.0) * vz3;

  double k_4_primex = primex(gamma,x3,y3,z3,vx3,vy3,vz3);
  double k_4_primey = primey(gamma,x3,y3,z3,vx3,vy3,vz3);
  double k_4_primez = primez(gamma,x3,y3,z3,vx3,vy3,vz3);


  //fourth step
 double average_k_x = (1.0/6.0)*(k_1_primex + 2.0*k_2_primex + 2.0*k_3_primex + k_4_primex);
 double average_k_y = (1.0/6.0)*(k_1_primey + 2.0*k_2_primey + 2.0*k_3_primey + k_4_primey);
 double average_k_z = (1.0/6.0)*(k_1_primez + 2.0*k_2_primez + 2.0*k_3_primez + k_4_primez);

 *t_new = t_old + h;
 *vxnew = vx_old + (h * average_k_x);
 *vynew = vy_old + (h * average_k_y);
 *vznew = vz_old + (h * average_k_z);
 *xnew = x_old + (h * (*vxnew));
 *ynew = y_old + (h * (*vynew));
 *znew = z_old + (h * (*vznew));
}
