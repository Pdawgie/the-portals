/* Consider Case 1 of the Dickey-Fuller (DF) test. Assume that the data generating */
/* process (d.g.p.) is a zero-mean AR(1) process, */
/* 𝑍𝑡 = 𝜌𝑍𝑡−1 + 𝑎𝑡 */

/* with 𝑍0 = 0. Let {𝑎𝑡 */
/* } be a Gaussian white noise process with mean 0 and variance */
/* 𝜎𝑎 */
/* 2 = 1. Recall that the null and alternative hypotheses in terms of 𝜌 are */
/* 𝐻0: 𝜌 = 1 */
/* 𝐻1: 𝜌 < 1 */
/* In this assignment the size, critical values and power of the regression coefficientbased test and the studentized test will be analyzed using Monte Carlo simulation. */
/* This is done by simulating a large number of AR(1) processes under certain */
/* assumptions and calculating 𝜌̂, 𝜈 and 𝜏 for each simulated process. The size, critical */
/* values and power of the tests are then obtained by analyzing the distributions of 𝜈 and */
/* 𝜏. */

/* (a) Use proc iml to simulate 10 000 AR(1) processes under 𝐻0, each with */
/* 𝑛 = 75. For each simulated AR(1) process, calculate 𝜌̂, 𝜈 and 𝜏. Create */
/* a SAS data set containing the 10 000 values of 𝜈 and 𝜏. */
/*  */
/* 𝐻0: 𝜌 = 1 */
/* 𝐻1: 𝜌 < 1 */

proc iml;
rep=10000; 
n=75; 
seed=0; 
matrix=J(rep,2,0);
do i=1 to rep; 
   x=J(n,1,0); 
   y=J(n,1,0); 
   do j=1 to n; 
      at=rannor(J(n,1,0)); 
      if j=1 then x[j,1]=at[j,1]; 
      if j=1 then y[j,1]=0; 
      if j>1 then x[j,1]=x[j-1,1]+at[j,1]; 
      if j>1 then y[j,1]=x[j-1,1]; 
  end;
  b=inv(t(x)*x)*t(x)*y; 
  yhat=x*b; 
  resid=y-yhat; 
  sse=t(resid)*resid; 
  df=n-1;
  mse=sse/df;
  rmse=sqrt(mse);
  stdb=sqrt(inv(t(x)*x)#mse);
  nu=(b-1)*n;
  tau=(b-1)/stdb;
  matrix[i,1]=nu;
  matrix[i,2]=tau;
end;
create matrix from matrix [colname={'nu','tau'}]; 
append from matrix; 
quit;

data Monte_Carlo_1;
set matrix;
if nu<-6.5 then size_v=1;
else size_v=0;
if tau<-1.5 then size_tau=1;
else size_tau=0;
run; 

proc sgplot data=Monte_Carlo_1;
  histogram tau;
  density tau;
run;

proc sgplot data=Monte_Carlo_1;
  histogram nu;
  density nu;
run;

proc means data=Monte_Carlo_1;
var size_v size_tau;
run;

proc univariate data= Monte_Carlo_1 alpha=0.01 alpha=0.05 alpha=0.1;
var nu tau;
run;

/* second simulation */

proc iml;
rep=10000; 
n=75; 
seed=0; 
matrix2=J(rep,2,0);
do i=1 to rep; 
   x=J(n,1,0); 
   y=J(n,1,0); 
   do j=1 to n; 
      at=rannor(J(n,1,0)); 
      if j=1 then x[j,1]=at[j,1]; 
      if j=1 then y[j,1]=0; 
      if j>1 then x[j,1]=0.9*x[j-1,1]+at[j,1]; 
      if j>1 then y[j,1]=x[j-1,1]; 
  end;
  b=inv(t(x)*x)*t(x)*y; 
  yhat=x*b; 
  resid=y-yhat; 
  sse=t(resid)*resid; 
  df=n-1;
  mse=sse/df;
  rmse=sqrt(mse);
  stdb=sqrt(inv(t(x)*x)#mse);
  nu=(b-1)*n;
  tau=(b-1)/stdb;
  matrix2[i,1]=nu;
  matrix2[i,2]=tau;
end;
create matrix2 from matrix2 [colname={'nu','tau'}]; 
append from matrix2; 
quit;

data Monte_Carlo_2;
set matrix2;
if nu<-6.5 then size_v=1;
else size_v=0;
if tau<-2.5 then size_tau=1;
else size_tau=0;
run; 

proc means data=Monte_Carlo_2;
var size_v size_tau;
run;

proc univariate data= Monte_Carlo_2 alpha=0.01 alpha=0.05 alpha=0.1;
var nu tau;
run;