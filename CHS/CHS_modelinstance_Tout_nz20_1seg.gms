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
          r2    "dt/dx     "
          v     0 Vx velocitgy;
          r1 = dt/(dx*dx) ;
          r2 = dt/dx ;
          k1 = 0.000 ;
          k0 = 0.600 ;
          c = 10;
          alpha = 1.0;
          v = 0;

Parameter
          Tout_vertex(t)
          ;
          Tout_vertex(t) = 1;

variables
          FId           flexibility index
          Tout(T)       uncertain parameter T_out
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
        Tout_uncertain(T)     uncertain parameter
        eobj                  objective function
        eq_IC(X,T)            initial condition setting
        BC_L1(T)
        BC_L2(T)
        BC_L3(T)
        BC_L4(T)
        BC_L5(T)
        BC_L6(T)
        BC_L7(T)
        BC_L8(T)
        BC_L9(T)
        BC_L10(T)
        BC_L11(T)
        BC_L12(T)
        BC_L13(T)
        BC_L14(T)
        BC_L15(T)
        BC_L16(T)
        BC_L17(T)
        BC_L18(T)
        BC_L19(T)
        BC_L20(T)
        BC_R1(T)
        BC_R2(T)
        BC_R3(T)
        BC_R4(T)
        BC_R5(T)
        BC_R6(T)
        BC_R7(T)
        BC_R8(T)
        BC_R9(T)
        BC_R10(T)
        BC_R11(T)
        BC_R12(T)
        BC_R13(T)
        BC_R14(T)
        BC_R15(T)
        BC_R16(T)
        BC_R17(T)
        BC_R18(T)
        BC_R19(T)
        BC_R20(T)
        For_u(X,T)
          
*       sym_bc(T)
        left_bc(T)
        right_bc(T) ;
    eq_IC(X,T)$(init_region(X,T)) .. u(X,T) =e= init_state(X) ;
    Tout_uncertain(T) ..  Tout(T) =e= Tout_n+0.5*FId*Tout_vertex(T) ;
    
*   forward
    For_u(X,T)$(region(X,T)) ..  u(X,T+1) =e= u(X,T)+alpha*(r1*u(X-1,T)-2*r1*u(X,T)+r1*u(X+1,T))-v*r2*(u(X+1,T)-u(X,T))-c*dt*(u(X,T)-Tout(T)) ;
*   central
*   For_u(X,T)$(region(X,T)) ..  u(X,T+1) =e= u(X,T-1)+alpha(T)*r1*(u(X-1,T)-2*u(X,T)+u(X+1,T))-v*r2*(u(X+1,T)-u(X-1,T))-c*2*dt*(u(X,T)-k0) ;
    
*    sym_bc(T).. bc_L(T) =e= bc_R(T);
    left_bc(T).. u('x1',T) =e= bc_L(T);
    right_bc(T).. u('x11',T) =e= bc_R(T);
    BC_L1(T)$(ord(T) lt 17).. BC_L(T) =e= BC_L(T+1);
    BC_L2(T)$((ord(T) gt 17) and (ord(T) lt 35))..                BC_L(T) =e= BC_L(T+1);
    BC_L3(T)$((ord(T) gt 35) and (ord(T) lt 52))..                BC_L(T) =e= BC_L(T+1);
    BC_L4(T)$((ord(T) gt 52) and (ord(T) lt 70))..                BC_L(T) =e= BC_L(T+1);
    BC_L5(T)$((ord(T) gt 70) and (ord(T) lt 87))..                BC_L(T) =e= BC_L(T+1);
    BC_L6(T)$((ord(T) gt 87) and (ord(T) lt 105))..                BC_L(T) =e= BC_L(T+1);
    BC_L7(T)$((ord(T) gt 105) and (ord(T) lt 122))..                BC_L(T) =e= BC_L(T+1);
    BC_L8(T)$((ord(T) gt 122) and (ord(T) lt 140))..                BC_L(T) =e= BC_L(T+1);
    BC_L9(T)$((ord(T) gt 140) and (ord(T) lt 157))..                BC_L(T) =e= BC_L(T+1);
    BC_L10(T)$((ord(T) gt 157) and (ord(T) lt 175))..                BC_L(T) =e= BC_L(T+1);
    BC_L11(T)$((ord(T) gt 175) and (ord(T) lt 192))..                BC_L(T) =e= BC_L(T+1);
    BC_L12(T)$((ord(T) gt 192) and (ord(T) lt 210))..                BC_L(T) =e= BC_L(T+1);
    BC_L13(T)$((ord(T) gt 210) and (ord(T) lt 227))..                BC_L(T) =e= BC_L(T+1);
    BC_L14(T)$((ord(T) gt 227) and (ord(T) lt 245))..                BC_L(T) =e= BC_L(T+1);
    BC_L15(T)$((ord(T) gt 245) and (ord(T) lt 262))..                BC_L(T) =e= BC_L(T+1);
    BC_L16(T)$((ord(T) gt 262) and (ord(T) lt 280))..                BC_L(T) =e= BC_L(T+1);
    BC_L17(T)$((ord(T) gt 280) and (ord(T) lt 297))..                BC_L(T) =e= BC_L(T+1);
    BC_L18(T)$((ord(T) gt 297) and (ord(T) lt 315))..                BC_L(T) =e= BC_L(T+1);
    BC_L19(T)$((ord(T) gt 315) and (ord(T) lt 332))..                BC_L(T) =e= BC_L(T+1);
    BC_L20(T)$((ord(T) gt 332) and (ord(T) lt 350))..                BC_L(T) =e= BC_L(T+1);
    BC_R1(T)$(ord(T) lt 17).. BC_R(T) =e= BC_R(T+1);
    BC_R2(T)$((ord(T) gt 17) and (ord(T) lt 35))..                BC_R(T) =e= BC_R(T+1);
    BC_R3(T)$((ord(T) gt 35) and (ord(T) lt 52))..                BC_R(T) =e= BC_R(T+1);
    BC_R4(T)$((ord(T) gt 52) and (ord(T) lt 70))..                BC_R(T) =e= BC_R(T+1);
    BC_R5(T)$((ord(T) gt 70) and (ord(T) lt 87))..                BC_R(T) =e= BC_R(T+1);
    BC_R6(T)$((ord(T) gt 87) and (ord(T) lt 105))..                BC_R(T) =e= BC_R(T+1);
    BC_R7(T)$((ord(T) gt 105) and (ord(T) lt 122))..                BC_R(T) =e= BC_R(T+1);
    BC_R8(T)$((ord(T) gt 122) and (ord(T) lt 140))..                BC_R(T) =e= BC_R(T+1);
    BC_R9(T)$((ord(T) gt 140) and (ord(T) lt 157))..                BC_R(T) =e= BC_R(T+1);
    BC_R10(T)$((ord(T) gt 157) and (ord(T) lt 175))..                BC_R(T) =e= BC_R(T+1);
    BC_R11(T)$((ord(T) gt 175) and (ord(T) lt 192))..                BC_R(T) =e= BC_R(T+1);
    BC_R12(T)$((ord(T) gt 192) and (ord(T) lt 210))..                BC_R(T) =e= BC_R(T+1);
    BC_R13(T)$((ord(T) gt 210) and (ord(T) lt 227))..                BC_R(T) =e= BC_R(T+1);
    BC_R14(T)$((ord(T) gt 227) and (ord(T) lt 245))..                BC_R(T) =e= BC_R(T+1);
    BC_R15(T)$((ord(T) gt 245) and (ord(T) lt 262))..                BC_R(T) =e= BC_R(T+1);
    BC_R16(T)$((ord(T) gt 262) and (ord(T) lt 280))..                BC_R(T) =e= BC_R(T+1);
    BC_R17(T)$((ord(T) gt 280) and (ord(T) lt 297))..                BC_R(T) =e= BC_R(T+1);
    BC_R18(T)$((ord(T) gt 297) and (ord(T) lt 315))..                BC_R(T) =e= BC_R(T+1);
    BC_R19(T)$((ord(T) gt 315) and (ord(T) lt 332))..                BC_R(T) =e= BC_R(T+1);
    BC_R20(T)$((ord(T) gt 332) and (ord(T) lt 350))..                BC_R(T) =e= BC_R(T+1);
    
    eobj .. obj =e= FId;
    
Model CHS /all/;

*$onecho >bench.opt
*  solvers conopt knitro minos snopt
*$offecho

*temp_slab_FId.optfile=1;
*option nlp=bench;
option reslim = 300;
Scalar ms 'model status', ss 'solve status' ;