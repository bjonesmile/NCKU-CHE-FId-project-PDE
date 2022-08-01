$ontext
An insulated bar with an initial temperature distribution at t=0, having ends
that are subsequently maintain at temperature which may be function of time.
PDE:
T = (theta - theta0)/(theta - theta0), tau = alpha*t/L^2, X = x/L
dT/dtau = d^2T/dX^2, for 0<x<1, 0<t<T
I.C.: T(X,0)=f(x), 0<=x<=1
B.C.: T(0,tau)=g0(tau), 0<t<=1 : T = 0, for 0<=X<=1
      T(1,tau)=g1(tau), 0<t<=1 : T = 1, for X=0 and X=1
$offtext

Set T /t1*t350/ ;
Set X /x1*x11/ ;

*Determination of zone for temperature distribute equation
Set region(X,T) network of grid point;
    region(X,T) = yes;
    region(X,T)$(ord(T)=card(T)) = no;
    region(X,T)$(ord(X)=1) = no;
    region(X,T)$(ord(X)=card(X)) = no;
    
Set init_region(X,T) for init setting;
    init_region(X,T) =no;
    init_region(X,T)$(ord(T)=1) = yes;
    init_region(X,T)$(ord(X)=1) = no;
    init_region(X,T)$(ord(X)=card(X)) = no;

Set limit_region(X,T) ;
    limit_region(X,T) = yes ;
    limit_region(X,T)$(ord(X)=1) = no;
    limit_region(X,T)$(ord(X)=card(X)) = no;

* Parameters
scalar dt  step space in x direction  ;
scalar dx  step space in y direction  ;
parameter init_state(X) slab initial temperature state
        /x1 0.4
         x2 0.49
         x3 0.56
         x4 0.61
         x5 0.64
         x6 0.65
         x7 0.64
         x8 0.61
         x9 0.56
         x10 0.49
         x11 0.4 /;

dt = 0.001;
dx = 0.1;

Scalars
          Tout_n  /0.75/
          k0    0-oreder heat source
          k1    1-oreder heat source
          alpha
          c     
          r1    "dt/(dx*dx)"
          r2    "dt/dx     ";
          r1 = dt/(dx*dx) ;
          r2 = dt/dx ;
          k1 = 0.000 ;
          k0 = 0.600 ;
          c = 10;
          alpha = 1.0;

Parameter
          Tout_x234_vertex(t)
          Tout_x567_vertex(t)
          Tout_x8910_vertex(t)
          ;
          Tout_x234_vertex(t) = 1;
          Tout_x567_vertex(t) = 1;
          Tout_x8910_vertex(t) = 1;

variables
          FId           flexibility index
          Tout(X,T)     uncertain parameter Tout
          obj           objective variable
          u(X,T)        temperature distribution
          bc_L(T)       bondary condition of left side
          bc_R(T)       bondary condition of right side ;
          
* Variable bounds and initialization
    u.l(X,T) = 0;
    u.l(X,T) = 0.8;
    u.lo(X,T) = 0.0;
    
    u.fx(X,T)$(init_region(X,T)) = init_state(X) ;
    loop((T),u.l(X,T)=init_state(X));
    FId.up = 1.0;
    FId.lo = 0;
* Boundary conditions
    bc_L.lo(T) = 0.3;
    bc_R.lo(T) = 0.3;
    bc_R.up(T) = 0.6;
    bc_L.up(T) = 0.6;
    
Equations
        Tout_L
        Tout_R
        Tout_uncertain_x2(X,T)        uncertain parameter Tout at x2
        Tout_uncertain_x3(X,T)        uncertain parameter Tout at x3
        Tout_uncertain_x4(X,T)        uncertain parameter Tout at x4
        Tout_uncertain_x5(X,T)        uncertain parameter Tout at x5
        Tout_uncertain_x6(X,T)        uncertain parameter Tout at x6
        Tout_uncertain_x7(X,T)        uncertain parameter Tout at x7
        Tout_uncertain_x8(X,T)        uncertain parameter Tout at x8
        Tout_uncertain_x9(X,T)        uncertain parameter Tout at x9
        Tout_uncertain_x10(X,T)       uncertain parameter Tout at x10
        eobj                          objective function
        eq_IC(X,T)                    initial condition setting
        BC_L1(T)
        BC_R1(T)
        BC_L2(T)
        BC_R2(T)
        BC_L3(T)
        BC_R3(T)
        BC_L4(T)
        BC_R4(T)
        BC_L5(T)
        BC_R5(T)
        For_u(X,T)
        left_bc(T)
        right_bc(T) ;
    eq_IC(X,T)$(init_region(X,T)) .. u(X,T) =e= init_state(X) ;
    Tout_L(X,T) .. Tout('x1',T) =e= Tout('x2',T) ;
    Tout_R(X,T) .. Tout('x10',T) =e= Tout('x11',T) ;
    Tout_uncertain_x2(X,T)  ..  Tout('x2',T) =e= Tout_n+0.5*FId*Tout_x234_vertex(T) ;
    Tout_uncertain_x3(X,T)  ..  Tout('x3',T) =e= Tout_n+0.5*FId*Tout_x234_vertex(T) ;
    Tout_uncertain_x4(X,T)  ..  Tout('x4',T) =e= Tout_n+0.5*FId*Tout_x234_vertex(T) ;
    Tout_uncertain_x5(X,T)  ..  Tout('x5',T) =e= Tout_n+0.5*FId*Tout_x567_vertex(T) ;
    Tout_uncertain_x6(X,T)  ..  Tout('x6',T) =e= Tout_n+0.5*FId*Tout_x567_vertex(T) ;
    Tout_uncertain_x7(X,T)  ..  Tout('x7',T) =e= Tout_n+0.5*FId*Tout_x567_vertex(T) ;
    Tout_uncertain_x8(X,T)  ..  Tout('x8',T) =e= Tout_n+0.5*FId*Tout_x8910_vertex(T) ;
    Tout_uncertain_x9(X,T)  ..  Tout('x9',T) =e= Tout_n+0.5*FId*Tout_x8910_vertex(T) ;
    Tout_uncertain_x10(X,T) ..  Tout('x10',T)=e= Tout_n+0.5*FId*Tout_x8910_vertex(T) ;
    
*   forward
    For_u(X,T)$(region(X,T)) ..  u(X,T+1) =e= u(X,T)+alpha*(r1*u(X-1,T)-2*r1*u(X,T)+r1*u(X+1,T))-c*dt*(u(X,T)-Tout(X,T)) ;
*   central
*   For_u(X,T)$(region(X,T)) ..  u(X,T+1) =e= u(X,T-1)+alpha(T)*r1*(u(X-1,T)-2*u(X,T)+u(X+1,T))-v*r2*(u(X+1,T)-u(X-1,T))-c*2*dt*(u(X,T)-k0) ;
    
    left_bc(T).. u('x1',T) =e= bc_L(T);
    right_bc(T).. u('x11',T) =e= bc_R(T);
    BC_L1(T)$(ord(T) lt 70).. bc_L(T) =e= bc_L(T+1);
    BC_R1(T)$(ord(T) lt 70).. bc_R(T) =e= bc_R(T+1);
    BC_L2(T)$((ord(T) gt 70) and (ord(T) lt 140)).. bc_L(T) =e= bc_L(T+1);
    BC_R2(T)$((ord(T) gt 70) and (ord(T) lt 140)).. bc_R(T) =e= bc_R(T+1);
    BC_L3(T)$((ord(T) gt 140) and (ord(T) lt 210)).. bc_L(T) =e= bc_L(T+1);
    BC_R3(T)$((ord(T) gt 140) and (ord(T) lt 210)).. bc_R(T) =e= bc_R(T+1);
    BC_L4(T)$((ord(T) gt 210) and (ord(T) lt 280)).. bc_L(T) =e= bc_L(T+1);
    BC_R4(T)$((ord(T) gt 210) and (ord(T) lt 280)).. bc_R(T) =e= bc_R(T+1);
    BC_L5(T)$((ord(T) gt 280) and (ord(T) lt 350)).. bc_L(T) =e= bc_L(T+1);
    BC_R5(T)$((ord(T) gt 280) and (ord(T) lt 350)).. bc_R(T) =e= bc_R(T+1);

    eobj .. obj =e= FId;
    
Model CHS /all/;

*$onecho >bench.opt
*  solvers conopt knitro minos snopt
*$offecho

*temp_slab_FId.optfile=1;
*option nlp=bench;
option reslim = 300;
Scalar ms 'model status', ss 'solve status' ;