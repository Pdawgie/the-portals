/*a*/
proc iml;
rep=10000;
n=75;
seed=0;
rho=1;
v2=j(rep,1,.);
t2=j(rep,1,.);
do sim=1 to rep;
	at=rannor(j(n,1,seed));
	zt=j(n,1,.);
	zt[1,1]=at[1,1];
	do t=2 to n;
		zt[t,1]=rho#zt[t-1,1]+at[t,1];
	end;
	y=zt[2:n,1];
	x=zt[1:n-1,1];
	b=inv(t(x)*x)*t(x)*y;
	yhat=x*b;
	resid=y-yhat;
	sse=t(resid)*resid;
	df=n-1;
	mse=sse/df;
	rmse=sqrt(mse);
	stdb=sqrt(inv(t(x)*x)#mse);
	rhohat=b;
	s_rhohat=stdb;
	v2[sim,1]=n#(rhohat-1);
	t2[sim,1]=(rhohat-1)/s_rhohat;
end;

v2_t2=v2||t2;
varname='v2'||'t2';
create H0true from v2_t2[colname=varname];
append from v2_t2;
/* print H0true; */
quit;

/*b*/
data size_testb;
set H0true;
size_v2=0;
size_t2=0;
if v2<-6.5 then size_v2=1;
if t2<-1.5 then size_t2=1;
run;

proc sgplot data=size_testb;
  histogram t2;
  density t2;
run;

proc sgplot data=size_testb;
  histogram v2;
  density v2;
run;

proc means data=size_testb;
	var size_v2 size_t2;
run;

/*c*/
proc univariate data=H0true;
	var v2 t2;
run;
/*d*/
proc iml;
rep=10000;
n=75;
seed=0;
rho=0.9;
v2=j(rep,1,.);
t2=j(rep,1,.);
do sim=1 to rep;
	at=rannor(j(n,1,seed));
	zt=j(n,1,.);
	zt[1,1]=at[1,1];
	do t=2 to n;
		zt[t,1]=rho#zt[t-1,1]+at[t,1];
	end;
	y=zt[2:n,1];
	x=zt[1:n-1,1];
	b=inv(t(x)*x)*t(x)*y;
	yhat=x*b;
	resid=y-yhat;
	sse=t(resid)*resid;
	df=n-1;
	mse=sse/df;
	rmse=sqrt(mse);
	stdb=sqrt(inv(t(x)*x)#mse);
	rhohat=b;
	s_rhohat=stdb;
	v2[sim,1]=n#(rhohat-1);
	t2[sim,1]=(rhohat-1)/s_rhohat;
end;

v2_t2=v2||t2;
varname='v2'||'t2';
create H1true from v2_t2[colname=varname];
append from v2_t2;
/*print nu_tau;*/
quit;

/*e*/
data powercv;
set H1true;
power_v=0;
power_t=0;
if v2<-6.5 then power_v=1;
if t2<-2.5 then power_t=1;
run;


goptions reset=all;
title1 'power of all the tests';

proc means data=powercv;
	var power_v power_t;
run;

/*f*/
proc iml;
rep=10000;
n=75;
seed=0;
rho=0.8;
v2=j(rep,1,.);
t2=j(rep,1,.);
do sim=1 to rep;
	at=rannor(j(n,1,seed));
	zt=j(n,1,.);
	zt[1,1]=at[1,1];
	do t=2 to n;
		zt[t,1]=rho#zt[t-1,1]+at[t,1];
	end;
	y=zt[2:n,1];
	x=zt[1:n-1,1];
	b=inv(t(x)*x)*t(x)*y;
	yhat=x*b;
	resid=y-yhat;
	sse=t(resid)*resid;
	df=n-1;
	mse=sse/df;
	rmse=sqrt(mse);
	stdb=sqrt(inv(t(x)*x)#mse);
	rhohat=b;
	s_rhohat=stdb;
	v2[sim,1]=n#(rhohat-1);
	t2[sim,1]=(rhohat-1)/s_rhohat;
end;

v2_t2=v2||t2;
varname='v2'||'t2';
create H1true from v2_t2[colname=varname];
append from v2_t2;
/*print nu_tau;*/
quit;

/*g*/
data powercv2;
set H1true;
power_v2=0;
power_t2=0;
if v2<-6.5 then power_v2=1;
if t2<-2.5 then power_t2=1;
run;

goptions reset=all;
title1 'power of the all the tests';

proc means data=powercv;
	var power_v2 power_t2;
run;
