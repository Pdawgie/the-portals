
ods graphics on;
proc contents data=work.yoyo ;
run;
proc means data=work.yoyo;
run;

PROC FREQ DATA=work.yoyo order=freq ;
    TABLES cut color clarity/ binomial missing plots=freqplot;
RUN;

title "Data manipulation for regression";
data yoyo_m1;
set work.yoyo;

ind_cut_good=0;
if cut= 'Good' then ind_cut_good=1;  
ind_cut_ideal=0;
if cut= 'Ideal' then ind_cut_ideal=1;
ind_cut_Premium=0;
if cut= 'Premium' then ind_cut_Premium=1;
ind_cut_Very_Good=0;
if cut= 'Very Good' then ind_cut_Very_Good=1;


ind_colorD=0;
if color='D' then ind_colorD=1; 
ind_colorE=0;
if color='E' then ind_colorE=1;
ind_colorF=0;
if color='F' then ind_colorF=1;
ind_colorG=0;
if color='G' then ind_colorG=1;
ind_colorH=0;
if color='H' then ind_colorH=1;
ind_colorI=0;
if color='I' then ind_colorI=1;

ind_clarity_IF =0;
if clarity= 'IF' then ind_clarity_IF=1; 
ind_clarity_SI1 =0; 
if clarity= 'SI1' then ind_clarity_SI1=1;
ind_clarity_SI2 =0;
if clarity= 'SI2' then ind_clarity_SI2=1;
ind_clarity_VS1 =0;
if clarity= 'VS1' then ind_clarity_VS1=1;
ind_clarity_VS2 =0;
if clarity= 'VS2' then ind_clarity_VS2=1;
ind_clarity_VVS1 =0;
if clarity= 'VVS1' then ind_clarity_VVS1=1;
ind_clarity_VVS2=0;
if clarity= 'VVS2' then ind_clarity_VVS2=1;

volume = (x*y)*z;

ln_volume=log(volume);

ln_price=log(price);

BC_price=price**0.5;
LN_BC_price=log(BC_price);

ln_carat=log(carat);
colors=color;

keep LN_BC_price ln_carat BC_price cut colors clarity price volume carat coconuts x y z ln_volume ln_price 
ln_carat ind_cut_good ind_cut_ideal ind_cut_Premium ind_cut_Very_Good  
ind_colorD ind_colorE ind_colorF ind_colorG ind_colorH ind_colorI  
ind_clarity_IF ind_clarity_SI1 ind_clarity_SI2 ind_clarity_VS1 ind_clarity_VS2 ind_clarity_VVS1 ind_clarity_VVS2;
run;

proc sgscatter data=yoyo_m1; 
title "Data scatter plot matrix";
matrix price carat volume coconuts  / diagonal=(histogram kernel) ellipse; 
run; 

proc corr data=yoyo_m1;
var price carat volume coconuts;
run;

proc transreg data=yoyo_m1 ss2 details ;
   model BoxCox(price) = identity( carat ind_cut_good ind_cut_ideal ind_cut_Premium ind_cut_Very_Good 
 ind_colorD ind_colorE ind_colorF ind_colorG ind_colorH ind_colorI  ind_clarity_IF ind_clarity_SI1 ind_clarity_SI2 
 ind_clarity_VS1 ind_clarity_VS2 ind_clarity_VVS1 ind_clarity_VVS2) ;
run;

proc transreg data=yoyo_m1 ss2 details ;
   model BoxCox(carat) = identity( LN_BC_price) ;
run;


PROC UNIVARIATE DATA=yoyo_m1 NORMAL PLOT;
VAR ln_BC_price ln_carat;
HISTOGRAM  ln_BC_price ln_carat  / normal kernel;
RUN;

title "Auxillary regression";
proc reg data=yoyo_m1 plots=fitplot ;
model ln_carat = coconuts ind_cut_good ind_cut_ideal ind_cut_Premium ind_cut_Very_Good 
 ind_colorD ind_colorE ind_colorF ind_colorG ind_colorH ind_colorI  ind_clarity_IF ind_clarity_SI1 ind_clarity_SI2 
 ind_clarity_VS1 ind_clarity_VS2 ind_clarity_VVS1 ind_clarity_VVS2 / vif collin clb;
run;


title "Regression first attempt";
proc reg data=yoyo_m1 plots=fitplot ;
Ln_BC_transformed_run_1: model LN_BC_price= coconuts ln_carat  ind_cut_good ind_cut_ideal ind_cut_Premium ind_cut_Very_Good 
 ind_colorD ind_colorE ind_colorF ind_colorG ind_colorH ind_colorI  ind_clarity_IF ind_clarity_SI1 ind_clarity_SI2 
 ind_clarity_VS1 ind_clarity_VS2 ind_clarity_VVS1 ind_clarity_VVS2/ vif collin clb ;
 
 title "Final Regression var coconuts removed";
proc reg data=yoyo_m1 plots=fitplot ;
 Ln_BC_Transformed: model LN_BC_price= ln_carat  ind_cut_good ind_cut_ideal ind_cut_Premium ind_cut_Very_Good 
 ind_colorD ind_colorE ind_colorF ind_colorG ind_colorH ind_colorI  ind_clarity_IF ind_clarity_SI1 ind_clarity_SI2 
 ind_clarity_VS1 ind_clarity_VS2 ind_clarity_VVS1 ind_clarity_VVS2/ vif collin clb ;
 output out=outs residual=ers p=ln_ythat;
run;

proc univariate data=outs normal  alpha=0.05;
var ers ln_ythat;
histogram ers /normal kernel;
histogram ln_ythat /normal kernel; 

run;

proc sgplot data=outs;
/* scatter x=ln_ythat y=ln_price;  */
scatter x=ln_ythat y=ln_BC_price;
ellipse x=ln_ythat y=ln_BC_price; 
run;

data test;
set outs;
pre=exp(ln_ythat);
pred=pre**2;
keep price pred;
run;

proc sgplot data=test;
scatter x=pred y=price;
ellipse x=pred y=price; 
run;













