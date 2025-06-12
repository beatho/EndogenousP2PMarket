%% Recup data
clear mpc

name = "TestFeeder"; 
moment = 1; % <1440

% case3 case9 case39 case85 case141 case_ACTIVSg200 case_ACTIVSg500
% case1888rte  case2383wp
mpopt = mpoption('verbose',1, 'pf.alg', 'NR');

%cas = name+ ".m";
%open(cas);
define_constants;
%oldmpc = loadcase(cas);
name1 = strcat('ACGrid/Agent', name, '.txt');
name2 = strcat('ACGrid/Branch', name, '.txt');
name3 = strcat('ACGrid/Case', name, '.txt');
name4 = strcat('ACGrid/AgentConsumption', name, '.txt');
%name4 = strcat('ACGrid/Bus', name, '.txt');
%name5 = strcat('ACGrid/Sol', name, '.txt');

AgentInfo = load(name1); % bus, a, b, P, Pmin, Pmax, Qobj, Qmin, Qma, zone
BranchInfo = load(name2); % from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;
CaseInfo = load(name3); % Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
ConsumptionInfo = load(name4);
%BusInfo = load(name4); % Gs, Bs, min, max, V0, theta0
%SolInfo = load(name5); % VM, VA, P, Q

SizeBus = 13;
SizeGen = 21;
SizeBranch = 13;
SizeGenCost = 7;

nAgent = CaseInfo(1,3);
nCons = CaseInfo(1,4);
nBus = CaseInfo(1,6);
nLine = CaseInfo(1,7);
V0 = CaseInfo(1,8);
theta0 = CaseInfo(1,9); 
Sb1= CaseInfo(1,1);
%% creation du cas

mpc.version = '2';
mpc.baseMVA = Sb1;
mpc.bus = zeros(nBus,SizeBus);
mpc.branch = zeros(nLine,SizeBranch);
mpc.gen = zeros(1,SizeGen);
mpc.gencost = zeros(1,SizeGenCost);

%% Bus %	bus_i	type	Pd	Qd	Gs	Bs	area Vm	Va baseKV zone	Vmax	Vmin

for bus=(1:nBus)
    mpc.bus(bus,1) = bus;
end
mpc.bus(:,2) = PQ;
mpc.bus(1,2) = REF;
mpc.bus(:,7) = 1;
mpc.bus(:,5) = 0;
mpc.bus(:,6) = 0;
mpc.bus(:,8) = 1;
mpc.bus(:,9) = 0;
mpc.bus(:,12) = 1.2;
mpc.bus(:,13) = 0.8;

mpc.bus(1, 12) = 1;
mpc.bus(1, 13) = 1;
%%
for i = 1:nAgent
    Pobj = -ConsumptionInfo(moment, i);
	bus = AgentInfo(i, 1);
	factor = AgentInfo(i, 2);
	PF = AgentInfo(i, 3);
    S0 = factor * Pobj;
	P0 = S0 * PF;
    Q0 = S0 * sqrt(1 - PF * PF);

    mpc.bus(bus,3) = mpc.bus(bus,3) - P0;
    mpc.bus(bus,4) = mpc.bus(bus,4) - Q0;
end

%% gen % bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf

mpc.gen(1,1) = 1;
mpc.gen(1,QMAX) = 100;
mpc.gen(1,QMIN) = -100;
mpc.gen(1,VG) = V0;
mpc.gen(1,GEN_STATUS) = 1;
mpc.gen(:,MBASE) = 100;
mpc.gen(1,PMIN) = -100;
mpc.gen(1,PMAX) = 100;


%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0

mpc.gencost(1, 1) = 2;
mpc.gencost(1, 4) = 3;
mpc.gencost(1, 5) = 0.01;
mpc.gencost(1, 6) = 0.03;


%% branch fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
% from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;
mpc.branch(:,F_BUS:T_BUS)= BranchInfo(:,1:2) + 1;

for l=(1:nLine)
    ys = BranchInfo(l,3) + 1j * BranchInfo(l,4);
    zs = 1/ys;
    r = real(zs);
    x = imag(zs);
    b = 2 * BranchInfo(l,5);
    tau = BranchInfo(l,6);
    theta = BranchInfo(l,7) * 180/pi;
    mpc.branch(l,BR_R) = r;
    mpc.branch(l,BR_X) = x;
    mpc.branch(l,BR_B) = b;
    mpc.branch(l,TAP) = tau;
    mpc.branch(l,SHIFT) = theta;
    mpc.branch(l,BR_STATUS) = 1;
    mpc.branch(l,ANGMIN) = -360;
    mpc.branch(l,ANGMAX) = 360;
end
%% Comparaison impédance


[YbusPerso, ~, ~] = makeYbus(ext2int(mpc));
YbusPerso2 = full(YbusPerso);

%%

savecase('testFeeder.m', mpc)


%%
tic
result = runpf(mpc,mpopt);
toc


SolInfo3 = zeros(nBus, 4);% VM, VA, P, Q

SolInfo3(:,1) = result.bus(:,VM);
SolInfo3(:,2) = result.bus(:,VA) *pi/180;
SolInfo3(:,3) = -result.bus(:,PD)/Sb1;
SolInfo3(:,4) = -result.bus(:,QD)/Sb1;


idGen = result.gen(1,BUS_I);
idBus = find(result.bus(:,1)==idGen);
SolInfo3(idBus,3) = SolInfo3(idBus,3) + result.gen(1,PG)/Sb1;
SolInfo3(idBus,4) = SolInfo3(idBus,4) + result.gen(1,QG)/Sb1;


%result.iterations
a = SolInfo3(1,3)
b = SolInfo3(1,4)

