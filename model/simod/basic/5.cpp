$PROB CEFTOZOLAN/TAZOBAKTAM

$CMT DEP1 PER1 DEP2 PER2

$PARAM
 X1 = 0.5;
 X2 = 0.5;
 X3 = 0.5;
 X4 = 0.5;
 X6 = 0.5;
 X7 = 0.5;
 X8 = 0.5;
 X9 = 0.5;


$MAIN
double CLC=X1;     
double QC=X2;
double V1C=X3;
double V2C=X4;
double CLT=X6;     
double QT=X7;
double V1T=X8;
double V2T=X9;
  
double K10 = CLC/V1C;
double K12 = QC/V1C;
double K21 = QC/V2C;
double K30 = CLT/V1T;
double K34 = QT/V1T;
double K43 = QT/V2T;

$DES
dxdt_DEP1= -K12*DEP1+K21*PER1-K10*DEP1;
dxdt_PER1=  K12*DEP1-K21*PER1 ;
dxdt_DEP2= -K34*DEP2+K43*PER2-K30*DEP2;
dxdt_PER2=  K34*DEP2-K43*PER2;
  
$TABLE
double  CEF=DEP1/V1C;
double TAZ=DEP2/V1T;

$CAPTURE @etas 1:LAST
CEF TAZ