function out1 = model1_dfdp_fn(in1,in2)
%MODEL1_DFDP_FN
%    OUT1 = MODEL1_DFDP_FN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    06-Feb-2014 03:10:09

A = in1(1,:);
B = in1(2,:);
Ka = in2(7,:);
Kb = in2(8,:);
k3 = in2(5,:);
m = in2(9,:);
n = in2(10,:);
nu = in2(1,:);
t2 = 1.0./Ka;
t3 = B.*t2;
t4 = t3.^m;
t5 = t4+1.0;
t6 = 1.0./t5.^2;
t7 = B.^n;
t8 = Kb.^n;
t9 = t7+t8;
t10 = 1.0./t9.^2;
t11 = log(B);
t12 = 1.0./t9;
out1 = reshape([1.0./t5,0.0,1.0,0.0,-A,0.0,0.0,A,0.0,A.*t7.*t12,0.0,-B,B.*1.0./Ka.^2.*m.*nu.*t3.^(m-1.0).*t6,0.0,0.0,-A.*Kb.^(n-1.0).*k3.*n.*t7.*t10,-nu.*t4.*t6.*log(t3),0.0,0.0,-A.*k3.*t7.*t10.*(t7.*t11+t8.*log(Kb))+A.*k3.*t7.*t11.*t12],[2, 10]);