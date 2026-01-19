$PROB vdv

$CMT DEP1 PER1 DEP2

$PARAM
THETA1 = 1;
THETA2 = 1;
THETA3 = 1;
THETA4 = 1;
THETA5 = 1;

$MAIN
double V1C=THETA1*exp(ETA(1)) ;      
double CLC=THETA2*exp(ETA(2)) ;
double V1T=THETA3*exp(ETA(3))     ;  
double CLT=THETA4*exp(THETA5*ETA(3)) ;

$DES

double AC1=DEP1/V1C;
double AT1=DEP2/V1T;

dxdt_DEP1=  - CLC * AC1;
dxdt_PER1=  0 ;
dxdt_DEP2=  - CLT * AT1;

$OMEGA @block 
0
0 0
0 0 0

$TABLE

double  CEF=DEP1/V1C;
double TAZ=DEP2/V1T;

$CAPTURE @etas 1:LAST
CEF TAZ