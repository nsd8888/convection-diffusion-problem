#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define grid 100

using namespace std;


int main()
{

  double T[grid][grid], Tn[grid-1][grid-1];

  int u,v,length,i,j,step;
  float dx,dy,F,D;
  length= 1;
  u=1;
  v=1;
  dx= length/(grid-1);
  dy= length/(grid-1);
  double error=0.00;
  float rho=1.2;


    for(j=0; j<=grid; j++)
    {
        T[0][j]= 0.0;
        T[grid][j]= 20.0;
    }

    for(i=0; i<=grid; i++)
    {
        T[i][0]= 20.0;
        T[i][grid]= 0.0;
    }

   F= (rho*u);
   D= (1.4/(1*dy));
   step=0;

while (error<=0.000001)
{
 for(i=1; i<=(grid-1); i++)
 {
   for(j=1; j<=(grid-1); j++)
   {  
   Tn[i][j]= (((D-(F/2))*dy)*T[i+1][j] + ((D+(F/2))*dy)*T[i-1][j] + ((D+(F/2))*dx)*T[i][j-1] + ((D-(F/2)*dx)*T[i][j+1]))/(((D-(F/2))*dy)+((D+(F/2))*dy)+((D+(F/2))*dx)+(D-(F/2)*dx));                                                                                                          
   }
   }

error=0;
for(i=1; i<=(grid-1); i++)
 {
   for(j=1; j<=(grid-1); j++)
   {  
      error+= abs(T[i][j]-Tn[i][j]);
   } 
 }

cout<<"the error in the step "<<step<<" is "<<error<<endl;

for(i=1; i<=(grid-1); i++)
 {
   for(j=1; j<=(grid-1); j++)
   { 
     T[i][j]= Tn[i][j];
   }
 } 

step++;

}
return 0;
}