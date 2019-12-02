$Title Gams code of supply chain network design problem using two-stage stochastic programming model



Sets
         i       "factory index"               /1*3/
         j       "Distribution center index"   /1*10/
         k       "customer zone index"         /1*32/
         s       "senario index"               /1*20/ ;
Alias (s, sPrime)
Positive variables
         x(i,j,s)        "amount of products shipped from factory i to DC j under scenario s"
         u(j,k,s)        "amount of products shipped from DC j to customer zone k under scenario s"
         epsilonplus(s)  "its hard to explain!"
         senarioObj(s)   "objective function for each senario"
         H(s)            "A variable to monitor the variance ish..."

        ;


Binary variable
         y(j)            "binary variable representing activation of DC" ;

Variables
         z         "objective function value" ;

Parameters
         d(k,s)          "demand of customer k under scenario s"
         c(i,j)          "transportation cost of unit product from factory i to DC j"
         a(j,k)          "transportation cost of unit product from DC j to customer zone k"
         b(k)            "shortage costs"
         lambda          "how important it is for us to decrease the expected objective"
         gama            "how important it is for us to decrease the variance of objective values"
         omega           "how bad it is not to satisfy constraints";
lambda = 1;
gama = 1;
omega = 10000;

Parameters  g(j) "fixed cost of opening cost"
         /1      813400
          2      600200
          3      853400
          4      935600
          5      933600
          6      777800
          7      793200
          8      980200
          9      891200
          10     806200 /

          cw(i)  "capacity of factory i"
         /1      30000
          2      32000
          3      20000/


         cy(j)           "capacity of DC j"
         /1      4067
          2      3001
          3      4267
          4      4678
          5      4668
          6      3889
          7      3966
          8      4901
          9      4456
          10     4031 /

         pi(s)        "probability of occurence of scenario s"

         /1        0.05
         2        0.05
         3        0.05
         4        0.05
         5        0.05
         6        0.05
         7        0.05
         8        0.05
         9        0.05
         10        0.05
         11        0.05
         12        0.05
         13        0.05
         14        0.05
         15        0.05
         16        0.05
         17        0.05
         18        0.05
         19        0.05
         20        0.05 /

         b(k)    "shortage costs"
              /   1        610
                  2        690
                  3        637
                  4        611
                  5        623
                  6        736
                  7        694
                  8        665
                  9        629
                  10        656
                  11        718
                  12        675
                  13        658
                  14        641
                  15        692
                  16        775
                  17        770
                  18        710
                  19        794
                  20        618
                  21        689
                  22        772
                  23        776
                  24        686
                  25        791
                  26        744
                  27        797
                  28        613
                  29        798
                  30        621
                  31        749
                  32        727    /;

* read some data from excel sheets to gdx file and load

$CALL GDXXRW.EXE stochastic.xlsx   par=d rng=d!B3:V36 Rdim=1 Cdim=1     par=c rng=c!B3:L7 Rdim=1 Cdim=1       par=a rng=a!B3:AH14 Rdim=1 Cdim=1

$GDXIN  stochastic.gdx
$LOAD  d, c, a
$GDXIN

Display d, c, a, b ;


* constraints declaration
Equations
         demand(k,s)             "demand satisfy constraint"
         balance(j,s)            "balance constraint"
         cap1(i,s)               "capacity limitation at factory i"
         cap2(j,s)               "capacity limitation at DC j"
         objec(s)                "say what senarioObj is"
         objective               "objective function"

         cons1(s)
                 ;


*Constraints definition
objective..      z =e= lambda * sum(s, senarioObj(s)* pi(s))+gama * sum(s, pi(s)*(senarioObj(s)-sum(sPrime, pi(sPrime)*senarioObj(sPrime))+2*h(s)))+omega* sum(s,pi(s)*epsilonplus(s))

                         ;

demand(k,s)..    sum(j,u(j,k,s))+epsilonplus(s) =g=  d(k,s) ;

balance(j,s)..   sum(i, x(i,j,s)) =e= sum(k, u(j,k,s)) ;

cap1(i,s)..      sum(j, x(i,j,s))=l= cw(i) ;

cap2(j,s)..      sum(i, x(i,j,s)) =l= cy(j)*y(j) ;

objec(s)..       senarioObj(s) =e= sum(j, g(j) * y(j))+sum((i,j), c(i,j) * x(i,j,s)) + sum((j,k), a(j,k) * u(j,k,s)) ;

cons1(s)..       0 =l= senarioObj(s)-sum(sPrime,pi(sPrime)*senarioObj(sPrime))+h(s);




Model stochastic /all/ ;
stochastic.reslim=10000;
stochastic.optCR= 0;
Solve stochastic using MIP minimizing z ;

Display z.l, y.l, x.l, u.l;
