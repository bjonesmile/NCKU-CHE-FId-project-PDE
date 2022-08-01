Sets
    t   /1*240/
    x   /1*5/
    j   /1*5/
    ;
    
Alias(j,k);

Scalar  nfe /5/; 


Table A(j,j)
     1        2         3        4         5
1   -13.0     14.7883  -2.6667   1.8784   -1.000
2   -5.3238   3.8730    2.0656  -1.2910    0.6762
3    1.5000  -3.2275    0.0000   3.2275   -1.5000
4   -0.6762   1.2910   -2.0656  -3.8730    5.3238
5    1.0000  -1.8784    2.6667  -14.7883  13.0000;

Table B(j,j)
    1           2           3           4           5
1   84.0000     -122.0632   58.6667     -44.6035    24.0000
2   53.2379     -73.3333    26.6667     -13.3333    6.7621
3   -6.0000     16.6667    -21.3333     16.6667     -6.0000
4   6.7621      -13.3333    26.6667     -73.3333    53.2379
5   24.0000     -44.6035    58.6667     -122.0632   84.0000;

Parameters
    I_n(t)
    I_vertex(t)
    ;

$gdxin Db_vertex.gdx
$load I_vertex
$gdxin

sets
    bc_in_region(t,x,j)
    bc_out_region(t,x,j)
    c_region(t,x,j)
    d_region(t,x,j)
    i_region(t,x,j)
    ;
    bc_in_region(t,x,j) = no;
    bc_in_region(t,x,j)$(ord(x) = 1 and ord(j) = 1) = yes;
    bc_out_region(t,x,j) = no;
    bc_out_region(t,x,j)$(ord(x) = card(x) and ord(x) = card(x)) = yes;
    
    c_region(t,x,j) = no;
    c_region(t,x,j)$(ord(j) = 1 or ord(j) = card(j)) = yes;
    c_region(t,x,j)$(ord(x) = 1 and ord(j) = 1) = no;
    c_region(t,x,j)$(ord(x) = card(x) and ord(j) = card(j)) = no;
    
    d_region(t,x,j) = yes;
    d_region(t,x,j)$(bc_in_region(t,x,j)) = no;
    d_region(t,x,j)$(bc_out_region(t,x,j)) = no;
    d_region(t,x,j)$(c_region(t,x,j)) = no;
    
    i_region(t,x,j) = no;
    i_region(t,x,j)$(ord(t) = 1) = yes;
    
parameters
    dt  /60/
    dx  /100/
    ;
    
Scalars
    eta         /0.4/
*oil
    rho_o       /932/
    cp_o        /2090/
    A_i         /1.9e-3/
    
    p_i         /0.157/
    r_i         /0.025/
*absorber pipe
    rho_A       /7850/
    cp_A        /460/
    A_o         /1.9e-3/
    p_o         /0.188/
    r_o         /0.03/
*glass envelope
    rho_E       /2400/
    cp_E        /840/
    A_E         /2.8e-3/
    
    h_air       /25/
    h_p         /3017/
    xi_A        /0.18/
    xi_E        /0.9/
    SB_const    /5.67e-8/
    W           /5.75/
    c
    
    T_air       /25/
    T_sky       /40/
    
    ;
    c = 1/xi_A*(1-xi_E)/xi_E*r_i/r_o;
    
    I_n(t)$(ord(t) <= 48) = 400;
    I_n(t)$(ord(t) > 48 and ord(t) <= 96) = 500;
    I_n(t)$(ord(t) > 96 and ord(t) <= 144) = 600;
    I_n(t)$(ord(t) > 144 and ord(t) <= 192) = 500;
    I_n(t)$(ord(t) > 192 and ord(t) <= 240) = 400;


Variables
    obj
    FId

    m(t)

    I(t,x,j)
    T_o(t,x,j)
    T_a(t,x,j)
    T_e(t,x,j)
    
    T_o_in(t)
    T_o_out(t)
    ;
    FId.l = 1.0;
    FId.lo = 0.0;
    FId.up = 1.25;

    m.lo(t) = 1.5;
    m.up(t) = 4;

    T_a.lo(t,x,j) = 100;
    T_a.up(t,x,j) = 400;
    T_a.l(t,x,j) = 220;
    T_a.fx('1',x,j) = 220;
    
    T_o.lo(t,x,j) = 200;
    T_o.up(t,x,j) = 280;
    T_o.l(t,x,j) = 220;
    T_o.fx('1',x,j) = 220;
    
    T_e.lo(t,x,j) = 25;
    T_e.l(t,x,j) = 35;
    T_e.fx('1',x,j) = 35;

    T_o_in.l(t) = 220;
    T_o_in.fx('1') = 220;
    
Parameters
    T_o_lb(t,x,j)
    T_o_ub(t,x,j)
    ;
$if not exist T_o_lb.gdx $abort 'no include To lower biund file for data file provided'
$gdxin T_o_lb.gdx
$load T_o_lb = T_o_lb
$gdxin
$if not exist T_o_ub.gdx $abort 'no include To upper biund file for data file provided'
$gdxin T_o_ub.gdx
$load T_o_ub = T_o_ub
$gdxin
    T_o.lo(t,x,j) = T_o_lb(t,x,j);
    T_o.up(t,x,j) = T_o_ub(t,x,j);

Equations
    fobj
    I_uncertain_x1_1(t)
    I_uncertain_x1_2(t)
    I_uncertain_x1_3(t)
    I_uncertain_x1_4(t)
    I_uncertain_x1_5(t)
    I_uncertain_x2_1(t)
    I_uncertain_x2_2(t)
    I_uncertain_x2_3(t)
    I_uncertain_x2_4(t)
    I_uncertain_x2_5(t)
    I_uncertain_x3_1(t)
    I_uncertain_x3_2(t)
    I_uncertain_x3_3(t)
    I_uncertain_x3_4(t)
    I_uncertain_x3_5(t)
    I_uncertain_x4_1(t)
    I_uncertain_x4_2(t)
    I_uncertain_x4_3(t)
    I_uncertain_x4_4(t)
    I_uncertain_x4_5(t)
    I_uncertain_x5_1(t)
    I_uncertain_x5_2(t)
    I_uncertain_x5_3(t)
    I_uncertain_x5_4(t)
    I_uncertain_x5_5(t)
*init    
    eq_init_T_o(t,x,j)
    eq_init_T_a(t,x,j)
    eq_init_T_e(t,x,j)
*C.V    
	m1(t)
	m2(t)
	m3(t)
	m4(t)
*S,V
    eq_PDE_To(t,x,j)
    eq_cont_To_x(t,x,j)
    eq_cont_To_dx(t,x,j)
    eq_BC_To_IN(t,x,j)
    eq_BC_To_OUT(t,x,j)
    
    eq_PDE_Ta(t,x,j)
    
    eq_PDE_Te(t,x,j)
    
    eq_T_o_in(t)
    eq_T_o_out(t)
    ;
    
    eq_init_T_o(t,x,j)$(i_region(t,x,j)).. T_o(t,x,j) =e= 220;
    eq_init_T_a(t,x,j)$(i_region(t,x,j)).. T_a(t,x,j) =e= 220;
    eq_init_T_e(t,x,j)$(i_region(t,x,j)).. T_e(t,x,j) =e= 35;
    
    eq_PDE_To(t,x,j)$(d_region(t,x,j) and ord(t) < card(t))..
    T_o(t+1,x,j) =e= T_o(t,x,j)+dt/2/(rho_o*cp_o*A_i)*(-m(t  )*cp_o*sum(k,A(j,k)*T_o(t  ,x,k))/(dx)+p_i*h_p*(T_a(t  ,x,j)-T_o(t  ,x,j))
                                                       -m(t+1)*cp_o*sum(k,A(j,k)*T_o(t+1,x,k))/(dx)+p_i*h_p*(T_a(t+1,x,j)-T_o(t+1,x,j)));
                                                      
    eq_cont_To_dx(t,x,j)$(ord(x) < card(x))..
    sum(k,A('5',k)*T_o(t,x,k)) =e=  sum(k,A('1',k)*T_o(t,x+1,k));
    eq_cont_To_x(t,x,j)$(ord(x) < card(x))..
    T_o(t,x,'5') =e= T_o(t,x+1,'1');
    eq_BC_To_IN(t,x,j)$(bc_in_region(t,x,j) and ord(t) > 1)..
    T_o(t,x,j)-T_o_in(t) =e= 0;
    eq_BC_To_OUT(t,x,j)$(bc_out_region(t,x,j) and ord(t) > 1)..
    sum(k,A(j,k)*T_o(t,x,k))/(dx) =e= 0;
    
    eq_PDE_Ta(t,x,j)$(ord(t) < card(t))..
    T_a(t+1,x,j) =e= T_a(t,x,j)+dt/2/(rho_a*cp_a*A_o)*(p_i*h_p*(T_o(t  ,x,j)-T_a(t  ,x,j))-SB_const/c*p_i*(T_a(t  ,x,j)**4-T_e(t  ,x,j)**4)+I(t  ,x,j)*eta*W
                                                      +p_i*h_p*(T_o(t+1,x,j)-T_a(t+1,x,j))-SB_const/c*p_i*(T_a(t+1,x,j)**4-T_e(t+1,x,j)**4)+I(t+1,x,j)*eta*W);
                                                      
    eq_PDE_Te(t,x,j)$(ord(t) < card(t))..
    T_e(t+1,x,j) =e= T_e(t,x,j)+dt/2/(rho_e*cp_e*A_e)*(SB_const/c*p_o*(T_a(t  ,x,j)**4-T_e(t  ,x,j)**4)-SB_const*xi_E*p_o*(T_e(t  ,x,j)**4-T_sky**4)-h_air*p_o*(T_e(t  ,x,j)-T_air)
                                                      +SB_const/c*p_o*(T_a(t+1,x,j)**4-T_e(t+1,x,j)**4)-SB_const*xi_E*p_o*(T_e(t+1,x,j)**4-T_sky**4)-h_air*p_o*(T_e(t+1,x,j)-T_air));
                                                      
    eq_T_o_out(t).. T_o_out(t) =e= T_o(t,'5','5');
    eq_T_o_in(t)$(ord(t) > 1).. T_o_in(t) =e=  0.995*T_o_in(t-1) + 0.005*(T_o_out(t-1)-90);

    I_uncertain_x1_1(t)  ..  I(t,'1','1') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x1_2(t)  ..  I(t,'1','2') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x1_3(t)  ..  I(t,'1','3') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x1_4(t)  ..  I(t,'1','4') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x1_5(t)  ..  I(t,'1','5') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x2_1(t)  ..  I(t,'2','1') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x2_2(t)  ..  I(t,'2','2') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x2_3(t)  ..  I(t,'2','3') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x2_4(t)  ..  I(t,'2','4') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x2_5(t)  ..  I(t,'2','5') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x3_1(t)  ..  I(t,'3','1') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x3_2(t)  ..  I(t,'3','2') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x3_3(t)  ..  I(t,'3','3') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x3_4(t)  ..  I(t,'3','4') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x3_5(t)  ..  I(t,'3','5') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x4_1(t)  ..  I(t,'4','1') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x4_2(t)  ..  I(t,'4','2') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x4_3(t)  ..  I(t,'4','3') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x4_4(t)  ..  I(t,'4','4') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x4_5(t)  ..  I(t,'4','5') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x5_1(t)  ..  I(t,'5','1') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x5_2(t)  ..  I(t,'5','2') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x5_3(t)  ..  I(t,'5','3') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x5_4(t)  ..  I(t,'5','4') =e= I_n(t)+300*FId*I_vertex(t) ;
    I_uncertain_x5_5(t)  ..  I(t,'5','5') =e= I_n(t)+300*FId*I_vertex(t) ;

	m1(t)$(ord(t) lt 60).. m(t) =e= m(t+1);
	m2(t)$((ord(t) gt 60) and (ord(t) lt 120)).. m(t) =e= m(t+1);
	m3(t)$((ord(t) gt 120) and (ord(t) lt 180)).. m(t) =e= m(t+1);
	m4(t)$((ord(t) gt 180) and (ord(t) lt 240)).. m(t) =e= m(t+1);

    fobj..  obj =e= FId;
        
    model PTC /all/;
    Scalar ms 'model status', ss 'solve status' ;
    
    option nlp = conopt;
    solve ptc using nlp maximizing obj;