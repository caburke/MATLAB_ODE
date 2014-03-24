function jac = FHN_jac_fn(in1,in2)
%FHN_JAC_FN
%    JAC = FHN_JAC_FN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    04-Feb-2014 15:34:30

V = in1(1,:);
b = in2(2,:);
c = in2(3,:);
t2 = 1.0./c;
jac = reshape([-c.*(V.^2-1.0),-t2,c,-b.*t2],[2, 2]);