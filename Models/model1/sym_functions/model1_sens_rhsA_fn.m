function SAdot = model1_sens_rhsA_fn(in1,in2,in3,in4)
%MODEL1_SENS_RHSA_FN
%    SADOT = MODEL1_SENS_RHSA_FN(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    06-Feb-2014 03:10:10

A = in1(1,:);
B = in1(2,:);
Ka = in4(7,:);
SA1 = in2(:,1);
SA2 = in2(:,2);
SA3 = in2(:,3);
SA4 = in2(:,4);
SA5 = in2(:,5);
SA6 = in2(:,6);
SA7 = in2(:,7);
SA8 = in2(:,8);
SA9 = in2(:,9);
SA10 = in2(:,10);
SB1 = in3(:,1);
SB2 = in3(:,2);
SB3 = in3(:,3);
SB4 = in3(:,4);
SB5 = in3(:,5);
SB6 = in3(:,6);
SB7 = in3(:,7);
SB8 = in3(:,8);
SB9 = in3(:,9);
SB10 = in3(:,10);
k1 = in4(3,:);
m = in4(9,:);
nu = in4(1,:);
t2 = 1.0./Ka;
t3 = B.*t2;
t4 = m-1.0;
t5 = t3.^t4;
t6 = t3.^m;
t7 = t6+1.0;
t8 = 1.0./t7.^2;
SAdot = [-SA1.*k1-SB1.*m.*nu.*t2.*t5.*t8,-SA2.*k1-SB2.*m.*nu.*t2.*t5.*t8,-SA3.*k1+1.0./t7-SB3.*m.*nu.*t2.*t5.*t8,-SA4.*k1-SB4.*m.*nu.*t2.*t5.*t8+1.0,-A-SA5.*k1-SB5.*m.*nu.*t2.*t5.*t8,-SA6.*k1-SB6.*m.*nu.*t2.*t5.*t8,-SA7.*k1-SB7.*m.*nu.*t2.*t5.*t8,-SA8.*k1-SB8.*m.*nu.*t2.*t5.*t8,-SA9.*k1+B.*1.0./Ka.^2.*m.*nu.*t5.*t8-SB9.*m.*nu.*t2.*t5.*t8,-SA10.*k1-SB10.*m.*nu.*t2.*t5.*t8];
