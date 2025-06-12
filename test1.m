%test

B1 = testBt1(1:11,1:172);
B2 = testBCPU2;

P1 = testP;
P2 = testPCPU2;

T1  = testT;
T2 = testTCPU2;

max(B1-B2);

figure()
hold on
plot(abs(B1'-B2'));

figure()
plot(abs(P1'-P2'))
figure()
plot(abs(T1'-T2'))
%%
P1-P2

%%
T1 - T2
