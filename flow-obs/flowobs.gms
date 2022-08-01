$ontext
Stationary flow of an incompressible fluid in a rectangular area in the
presence of an obstacle.
$offtext

set X /u1*u20/;
set Y /u1*u20/;

* Determination of zone for water movement equation
Set vyside(X,Y);
Set vxside(X,Y);
    vxside(X,Y)                  = yes;
    vxside(X,Y)$(ord(X)=1)       = no;
    vxside(X,Y)$(ord(X)=card(X)) = no;
    vxside(X,Y)$(ord(Y)=1)       = no;
    vxside(X,Y)$(ord(Y)=card(Y)) = no;
    vxside('u10','u10')            = no;
    vxside('u10','u11')            = no;
    vyside(X,Y) = vxside(X,Y);

* Parameters
scalar dx  step space in x direction  ;
scalar dy  step space in y direction  ;
scalar r   density of the fluid        /1000/;

$if not set gdxincname $abort 'no include file name for data file provided'
$gdxin %gdxincname%
$load dx dy
$gdxin

parameter m(X,Y) kinematic viscosity ;
          m(X,Y) := 0.5;

variables
  obj        objective variable
  D(X,Y)     error
  P(X,Y)     pressure
  Vx(X,Y)    x-direction velocity
  Vy(X,Y)    y-direction velocity  ;

* Variable bounds and initialization
  Vx.l(X,Y)     =  0.0;
  Vy.l(X,Y)     =  0.0;
  Vy.l(x,y)$(vxside(x,y)) = 0.0;

  D.lo(X,Y)     = 0.0;
  D.up(X,Y)     = 10.0;

  P.up(X,Y)     =  20000;
  P.lo(X,Y)     = -20000;
  P.l(X,Y)      =  0.0;

*Boundary conditions
  Vx.fx('u1',Y)   = 0.5;
  Vx.fx('u20',Y)  = 0.5;
  Vx.fx(X,'u1')   = 0;
  Vx.fx(X,'u20')  = 0;
  Vy.fx('u1',Y)   = 0;
  Vy.fx('u20',Y)  = 0;
  Vy.fx(X,'u1')   = 0;
  Vy.fx(X,'u20')  = 0;

* Obstacle description
  vx.fx('u10','u10') = 0;
  vx.fx('u10','u11') = 0;

Equations
     For_Vx(X,Y)
     For_Vy(X,Y)
     Div_Vxy(X,Y)
     eobj  objective function ;

For_Vx(X,Y)$(vxside(X,Y))..
   (P(X+1,Y)-P(X,Y))/(r*dx) =e=
     m(X,Y)*((Vx(X+1,Y)-2*Vx(X,Y)+Vx(x-1,Y))/(dx*dy)
    +        (Vx(X,Y+1)-2*Vx(X,Y)+Vx(X,Y-1))/(dx*dy));

For_Vy(X,Y)$(vxside(X,Y))..
   (P(X,Y+1)-P(X,Y))/(r*dy) =e=
     m(X,Y)*((Vy(X+1,Y)-2*Vy(X,Y)+Vy(X-1,Y))/(dx*dy)
    +        (Vy(X,Y+1)-2*Vy(X,Y)+Vy(X,Y-1))/(dy*dx));
*--
*For_Vx(X,Y)$(Vxside(X,Y))..
*    Vx(X,Y)*(Vx(X+1,Y)-Vx(X-1,Y))/(2*dx) +
*    0.25*(Vy(X+1,Y-1)+Vy(X+1,Y)+Vy(X,Y-1)+Vy(X,Y)) *
*         (Vx(X,Y+1)-Vx(X,Y-1))/(2*dy) +
*         (P(X+1,Y)-P(X,Y))/(r*dx)
*    =e=
*    m(X,Y)*((Vx(X+1,Y)-2*Vx(X,Y)+Vx(X-1,Y))/(dx*dx) +
*            (Vx(X,Y+1)-2*Vx(X,Y)+Vx(X,Y-1))/(dy*dy));
*
*For_Vy(X,Y)$(Vyside(X,Y))..
*    0.25*(Vx(X-1,Y+1)+Vx(X-1,Y)+Vx(X,Y+1)+Vx(X,Y)) *
*         (Vy(X+1,Y)-Vy(X-1,Y))/(2*dy) +
*         Vy(X,Y)*(Vy(X,Y+1)-Vy(X,Y-1))/(2*dy) +
*         (P(X,Y+1)-P(X,Y))/(r*dy)
*    =e=
*    m(X,Y)*((Vy(X+1,Y)-2*Vy(X,Y)+Vy(X-1,Y))/(dx*dx) +
*            (Vy(X,Y+1)-2*Vy(X,Y)+Vy(X,Y-1))/(dy*dy));
*--
Div_Vxy(X,Y)$((ord(X) > 1) $ (ord(Y) > 1))..
    (Vx(X,Y)-Vx(X-1,Y))/dx + (Vy(X,Y)-Vy(X,Y-1))/dy =e=
         D(X,Y);

eobj.. obj =e= SUM((X,Y),(D(X,Y)*D(X,Y)));

Model flowobs /all/;

$onecho >bench.opt
  solvers conopt knitro minos snopt
$offecho

flowobs.optfile=1;
option nlp=bench;

$gdxout tran
$unload X Y vxside vyside Vx Vy P D m
$gdxout

Solve flowobs using nlp minimizing obj;
execute_unload "flowobs_result.gdx", X, Y, Vx, Vy, vxside, vyside, P, D, m;
* end flowobs