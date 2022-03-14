proc iml; 
use work.query;
read all into xy;
n=nrow(xy);
s2=sample(1:n,2,"NoReplace");
/* print s2; */
/* cen=xy[s2, ]; */
x1=xy[,1];
x2=xy[,2];

nx1={10,20,49};
nx2={30,40,40};
x1cen=x1[nx1,];
x2cen=x2[nx2,];
cen=x1cen||x2cen;
icen=cen;
print icen;
con=1010101;

do tt=1 to n while (con>0);
print "iteration:" tt;

v1=cen[1,]//xy;
d1=(distance(v1,"L2"))[,1];
v2=cen[2,]//xy;
d2=(distance(v2,"L2"))[,1];
v3=cen[3,]//xy;
d3=(distance(v3,"L2"))[,1];
d=(d1||d2||d3)[2:nrow(v1),];
md=d[,><];
mdi=d[,>:<];
gr=(d=(md));

mean_cen1=(gr[,1]#xy)[+,]/(gr[,1])[+,];
mean_cen2=(gr[,2]#xy)[+,]/(gr[,2])[+,];
mean_cen3=(gr[,3]#xy)[+,]/(gr[,3])[+,];

cenres=cenres//(i||mean_cen1||mean_cen2);

print cenres;

ocen=cen;
cen=mean_cen1//mean_cen2//mean_cen3;
con=(abs(ocen-cen))[<>];
print con cen;
end;

/*  */
/* ;.................................................... */

pdat=xy||mdi;
call sort(pdat,{2});
/* print pdat; */
xy=pdat;
cnt=J(1,nrow(gr),1)*gr;
print cnt;

cnt=cusum(cnt);
print cnt;
te=(xy[,1:2]-repeat(J(1,n,1/n)*xy[,1:2],n))##2;
tte=sum((xy[,1:2]-repeat(J(1,n,1/n)*xy[,1:2],n))##2);
print tt;

cdat=pdat[,1:2];
mdat=cen[pdat[,3],];

www=(cdat-mdat)##2;
w=sum((cdat-mdat)##2);

r2=1-w/tte;
psF=((tte-w)/nrow(cen)-1)/(w/(n-nrow(cen)));

print w r2 psF;

quit;

ods graphics on;
proc cluster data = work.query method = centroid ccc print=30 outtree=Tree;
var x1 x2;
run;
ods graphics off;

proc tree noprint ncl=3 out=out;
copy x1 x2;
run;

proc candisc out = can;
class cluster;
var x1 x2;
run;
proc sgplot data = can;
title "Cluster Analysis for q1.sas7bdat 3 clusters";
scatter y = can2 x = can1 / group = cluster;
label can2="X2" can1="X1";

run;


