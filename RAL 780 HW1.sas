data q1;
input y x;
datalines;
16.0    1.0
9.0    0.0
17.0    2.0
12.0    0.0
22.0    3.0
13.0    1.0
8.0    0.0
15.0    1.0
19.0    2.0
11.0    0.0
;

run;


proc iml;
use q1;
read all into matrix;
y=matrix[,1];
x=matrix[,2];
Xbar=(1/10)*sum(x);
Ybar=(1/10)*sum(y);
b1num=sum((x-xbar)#(y-ybar));
b1denum=sum((x-xbar)##2);
b1=b1num/b1denum;
print b1;
b0=ybar-b1*xbar;
print b0;
/* therfore Yhat=10.2+4X; */
Yhat=10.2+4#X;
print yhat;

SSE=sum((Y-Yhat)##2);
SSR=sum((Yhat-Ybar)##2);
SSTO=sum((Y-Ybar)##2);
SSTOU=sum((Y)##2);

SSTOv2= sum((Y)##2-(10#Ybar)##2);
/* query */
print SSE SSR SSTO SSTOU SSTOv2;
df=10-2;
MSE= SSE/df;
MSR= SSR/1;

R2=SSR/SSTO;

r=(b1/b1)*Sqrt(r2);

print mse msr r2 r;
/* Very strong, positive linear relationship between X and Y */

/* H0 : b1 = 0 */
/* H1 : b1 != 0 */

F=MSR/MSE;

Print f;

seb1=sqrt(mse)*sqrt((sum((x)##2))/(10*sum((x)##2)));
print seb1;

t=b1/seb1;
print t;

ucl=b1+tinv(0.975,8)*seb1;
lcl=b1-tinv(0.975,8)*seb1;
print ucl lcl;




ods graphics on;
proc reg data=q1;
model y=x / alpha=0.1 CLB;
model y=x / CLM CLI;
id x;
run;
ods graphics off;

















