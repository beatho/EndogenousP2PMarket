%% create data for my PF
name = "case10ba"; 

% dist : case4_dist case10ba case85
% case3 case9 case39 case85 case141 case_ACTIVSg200 case_ACTIVSg500 case1888rte 
% 1951rte 2848rte 2868rte


%  case14 case18 case_ieee30  case57 case69  case118  case145  
%case1888rte case2383wp case300  case_ACTIVSg2000 case136ma case_ACTIVSg10k
% case_ACTIVSg25k case_ACTIVSg70k case3012wp
cas = name + ".m";
mpopt = mpoption('verbose',1, 'pf.alg', 'NR'); % GS NR
open(cas);
name1 = strcat('ACGrid/Agent', name, '.txt');
name2 = strcat('ACGrid/Branch', name, '.txt');
name3 = strcat('ACGrid/Case', name, '.txt');
name4 = strcat('ACGrid/Bus', name, '.txt');
name5 = strcat('ACGrid/Sol', name, '.txt');
name6 = strcat('ACGrid/Bgrid', name, '.txt');
name7 = strcat('ACGrid/Ggrid', name, '.txt');
define_constants;
oldmpc = loadcase(cas);
mpc = loadcase(cas);
oldref = find(mpc.bus(:,BUS_TYPE)==REF,1);
oldidref = mpc.bus(oldref,1);
idref = mpc.bus(1,BUS_I) ;
ref = 1;


if(oldidref ~=1)
    % changement de la numerotation des bus
    mpc.bus(oldref, 1) = idref;
    mpc.bus(ref,1) = oldidref;
    temp = mpc.bus(oldref, :);
    mpc.bus(oldref, :) = mpc.bus(ref, :);
    mpc.bus(ref,:) = temp;
    % changement des generateurs
    i = find(mpc.gen(:,GEN_BUS)==oldidref);
    k = find(mpc.gen(:,GEN_BUS)==idref);
    for indice=i
        mpc.gen(indice,GEN_BUS) = idref;
    end
    for indice=k
        mpc.gen(indice,GEN_BUS) = oldidref;
    end
    % changement des branches
    ifr = find(mpc.branch(:,F_BUS)==oldidref);
    ito = find(mpc.branch(:,T_BUS)==oldidref);
    kfr = find(mpc.branch(:,F_BUS)==idref);
    kto = find(mpc.branch(:,T_BUS)==idref);

    for indice=ifr
        mpc.branch(indice,F_BUS)=idref;
    end
    for indice=kfr
        mpc.branch(indice,F_BUS)=oldidref;
    end
    for indice=ito
       mpc.branch(indice,T_BUS)=idref;
    end
    for indice=kto
       mpc.branch(indice,T_BUS)=oldidref;
    end

end
idGenRef=find(mpc.gen(:,GEN_BUS)==idref);

test = mpc; 
test2 = mpc; 
tic
%profile on
Solmpc = runpf(mpc);
%profile viewer
toc
%%
disp("*********************************************************************************");
 %desequilibre et intro
Pgen = Solmpc.gen(1:end,PG);
Pgen(idGenRef) = 0;
test2.gen(:,PG) = Pgen;
Pbus = Solmpc.bus(1:end,PD);
Pdes = sum(Pgen) - sum(Pbus);

Qgen = Solmpc.gen(:,QG);
Qgen(idGenRef) = 0;
test2.gen(:,QG) = Qgen;
Qbus = Solmpc.bus(:,QD);
Qdes = sum(Qgen) - sum(Qbus);


test2.bus(ref,QD) = test2.bus(ref,QD) + Qdes;
test2.bus(ref,PD) =  test2.bus(ref,PD) + Pdes;

V0 = Solmpc.bus(ref,VM);
theta0 = Solmpc.bus(ref,VA);

% adaptation du cas (cela ne converge plus si on fait les 2...)
% test2.bus(:,VM)= V0;
% test2.bus(:,VA) = 0;
% test2.gen(:,VG)= V0;
test2.bus(:,BUS_TYPE) = PQ;
test2.bus(ref,BUS_TYPE) = REF;

tic
Soltest2 = runpf(test2);
toc
disp("********************************************************************************");

%%
% if oldref~=1
%     mpc.bus(1,BUS_TYPE) = REF;
%     %mpc.bus(1,VM) = 1;
%     %mpc.bus(1,VA) = 0;
%     mpc.bus(oldref,BUS_TYPE) = PQ;
%     oldidGenRef = find(mpc.gen(:,BUS_I)==oldidref,1);
%     
%     %desequilibre et intro
%     Pgen = mpc.gen(:,PG);
%     Pbus = mpc.bus(:,PD);
%    
%     taille = size(mpc.gen);
%     Ngen = taille(1);
%     nParamGen = taille(2);
%     Pdes = sum(Pgen) - sum(Pbus);
%     
%     mpc.bus(1,PD) =  mpc.bus(1,PD) + Pdes;
%     
%     
%     idGenRef = find(mpc.gen(:,BUS_I)==idref,1);
%     if isempty(idGenRef) % il faut ajouter un generateur sur le noeuf de ref
%         %nouveau gen
%         oldmpcGen = mpc.gen;
%         newmpcGen = zeros(Ngen+1,nParamGen);
%         newmpcGen(2:end,:) = oldmpcGen;
%         newmpcGen(1,:) = oldmpcGen(oldidGenRef,:);
%         newmpcGen(1,1) = idref;
%         newmpcGen(1,PG) = 0;
%         newmpcGen(1,QG) = 0;
%         newmpcGen(1,VG) = 1;
% 
%         mpc.gen = newmpcGen;
% 
%         %nouveau cout de gen
%         oldmpcGenCost = mpc.gencost;
%         newmpcGenCost = zeros(Ngen+1,7);
%         newmpcGenCost(2:end,:) = oldmpcGenCost;
%         newmpcGenCost(1,:) = oldmpcGenCost(oldidGenRef,:);
%         mpc.gencost = newmpcGenCost;
% 
%         % nouveau type de fuel et type
%         mpc.gentype={};
%         mpc.genfuel={};
% 
%         mpc = rmfield(mpc,{'gentype','genfuel'});
%     else
%         mpc.gen(idGenRef,VG) = 1;
%         
%     end
% else
%       %desequilibre et intro
%     Pgen = mpc.gen(:,PG);
%     Pbus = mpc.bus(:,PD);
%     Pdes = sum(Pgen) - sum(Pbus);
%     mpc.bus(1,PD) =  mpc.bus(1,PD) + Pdes;
%     Pdes=0;
%     
% 
% end


Ub1 = test2.bus(1,BASE_KV);
if(Ub1==0)
    Ub1=1;
end
Sb1 = test2.baseMVA;
Zb1 = Ub1*Ub1/Sb1;

taille = size(test2.bus);

NBus = taille(1);
taille = size(test2.gen);
Ngen = taille(1);
Pgen = Soltest2.gen(:,PG);
Qgen = Soltest2.gen(:,QG);

Pbus = Soltest2.bus(:,PD);
Qbus = Soltest2.bus(:,QD);
busRed = Soltest2.bus(Pbus>0 | (Pbus==0 & Qbus~=0),:);
busRed2 = Soltest2.bus(Pbus<0,:);
nGenSup = length(Pbus(Pbus<0)); % generateur PQ ?
nCons = length(Pbus(Pbus>0 | (Pbus==0 & Qbus~=0)));

% on produit Ploss de "trop" (car cela va dans les lignes), donc il faut
% enlever Ploss de la production d'un des générateurs, pour que l'agent
% qui s'occupe des pertes force tous les producteurs à participer.
Ploss = sum(Pgen) - sum(Pbus);
coresBus = find(test2 ...
    .bus(:,1)); % coresBus(i) : new id of mpc.bus(i,1) !
% ext2int(mpc) pour numeroter les bus peut marcher aussi ?


Qloss = sum(Qgen) - sum(Qbus);
%idem
%% Agent
AgentInfo = zeros(nCons+Ngen+nGenSup, 10); % bus, a, b, P, Pmin, Pmax, Qobj, Qmin, Qma, zone

% Conso
P = -busRed(:,PD);
Q = -busRed(:,QD);
a = 0.1 * ones(nCons,1);
b = -P.*a; 
Pmin = 1.2*P; % P<0 !
Pmax = 0.8*P;
Qmin = Q .* (1.05 - 0.1*(Q>0));
Qmax = Q .* (0.95 + 0.1*(Q>0));

debut = 1;
fin = nCons;

AgentInfo(debut:fin,1) = coresBus(Pbus>0 | (Pbus==0 & Qbus~=0));
AgentInfo(debut:fin,2) = a;
AgentInfo(debut:fin,3) = b;
AgentInfo(debut:fin,4) = P;
AgentInfo(debut:fin,5) = Pmin;
AgentInfo(debut:fin,6) = Pmax;
AgentInfo(debut:fin,7) = Q;
AgentInfo(debut:fin,8) = Qmin;
AgentInfo(debut:fin,9) = Qmax;
AgentInfo(debut:fin,10) = busRed(:,BUS_AREA);

% Prod peu controlable (se comporte comme un conso)
debut = fin + 1;
fin = debut + nGenSup -1;

P = -busRed2(:,PD);
Q = -busRed2(:,QD);
a = 0.1 * ones(nGenSup,1);
b = -P.*a; 
Pmin = 0.8*P; % P>0 !
Pmax = 1.2*P;
Qmin = Q .* (1.05 - 0.1*(Q>0));
Qmax = Q .* (0.95 + 0.1*(Q>0));

AgentInfo(debut:fin,1) = coresBus(Pbus<0);
AgentInfo(debut:fin,2) = a;
AgentInfo(debut:fin,3) = b;
AgentInfo(debut:fin,4) = P;
AgentInfo(debut:fin,5) = Pmin;
AgentInfo(debut:fin,6) = Pmax;
AgentInfo(debut:fin,7) = Q;
AgentInfo(debut:fin,8) = Qmin;
AgentInfo(debut:fin,9) = Qmax;
AgentInfo(debut:fin,10) = busRed2(:,BUS_AREA);




% Prod

debut = fin + 1;
fin = debut + Ngen -1;

P = Pgen;
P(idGenRef) = P(idGenRef) - Ploss;
if P(idGenRef) <0
    P(idGenRef) = 0;
end
Q = Qgen;
%Q(idGenRef) = Q(idGenRef) - Qloss;
try
    a = mpc.gencost(1:Ngen,5);
    b = mpc.gencost(1:Ngen,6);
catch
    a = 0.1 * ones(Ngen,1);
    b = 1 * ones(Ngen,1);
end
a(a==0) = 0.1;

Pmin = mpc.gen(:,PMIN);
Pmax = mpc.gen(:,PMAX);
Qmin = mpc.gen(:,QMIN);
Qmin(Qmin<-100000) = -100000;
Qmax = mpc.gen(:,QMAX);
Qmax(Qmax>100000) = 100000;
for i=1:Ngen
    bus = coresBus(mpc.bus(:,1) == mpc.gen(i,BUS_I));
    AgentInfo(debut + i-1,1) = bus;
end

AgentInfo(debut:fin,2) = a;
AgentInfo(debut:fin,3) = b;
AgentInfo(debut:fin,4) = P;
AgentInfo(debut:fin,5) = Pmin;
AgentInfo(debut:fin,6) = Pmax;
AgentInfo(debut:fin,7) = Q;
AgentInfo(debut:fin,8) = Qmin;
AgentInfo(debut:fin,9) = Qmax;

%% Branch

nLine = length(mpc.branch(:,1));

lineInfo = zeros(nLine, 10); % from, to, Ys Real, Ys Im, Yp, tau, theta, Limit=0, zs Real, zs Imag;

for l=(1:nLine)
    i = coresBus(mpc.bus(:,1) == mpc.branch(l,F_BUS));
    j = coresBus(mpc.bus(:,1) == mpc.branch(l,T_BUS));
    r = mpc.branch(l,BR_R);
    x = mpc.branch(l,BR_X);
    b = mpc.branch(l,BR_B);

    tau = mpc.branch(l,TAP); 
    theta = mpc.branch(l,SHIFT);
    
    z = r + 1i*x;
    y = 1/z;
    lineInfo(l,1)  = i;
    lineInfo(l,2)  = j;
    lineInfo(l,3)  = real(y);
    lineInfo(l,4)  = imag(y);
    lineInfo(l,5)  = b/2;
    lineInfo(l,6)  = tau;
    lineInfo(l,7)  = theta*pi/180;
    lineInfo(l,8)  = 0;
    lineInfo(l,9)  = real(z);
    lineInfo(l,10) = imag(z);


end


%% Constraint

ConstraintInfo = zeros(NBus,6); % Gs, Bs, min, max, V0, theta0

ConstraintInfo(:,1) = mpc.bus(:,GS);
ConstraintInfo(:,2) = mpc.bus(:,BS);

ConstraintInfo(:,3) = mpc.bus(:,VMIN);
ConstraintInfo(:,4) = mpc.bus(:,VMAX);

for i = (1:NBus)
    if ConstraintInfo(i,4) < Solmpc.bus(i,VM)
        ConstraintInfo(i,4) = floor(Solmpc.bus(i,VM)*10 + 1)/10;
    end
    if ConstraintInfo(i,3) > Solmpc.bus(i,VM)
        ConstraintInfo(i,3) = floor(Solmpc.bus(i,VM)*10)/10;
    end
end


ConstraintInfo(:,5) = mpc.bus(:,VM);
ConstraintInfo(:,6) = mpc.bus(:,VA) * pi/180;

for i=(1:Ngen)
    bus = mpc.gen(i,GEN_BUS);
    idbus = find(mpc.bus(:,BUS_I)==bus);
    ConstraintInfo(idbus,5) = mpc.gen(i,VG);
end

%% Case Info

CaseInfo    = zeros(1,9); % Sbase, Vbase, nAgent, nCons,nGenSup, nBus, nLine, V0, theta0
CaseInfo(1) = Sb1;
CaseInfo(2) = Ub1;
CaseInfo(3) = nCons+Ngen+nGenSup;
CaseInfo(4) = nCons;
CaseInfo(5) = nGenSup;
CaseInfo(6) = NBus;
CaseInfo(7) = nLine;
CaseInfo(8) = V0;
CaseInfo(9) = theta0 *pi/180;

%% PF solution

SolInfo = zeros(NBus, 4);% VM, VA, P, Q

SolInfo(:,1) = Solmpc.bus(:,VM);
SolInfo(:,2) = Solmpc.bus(:,VA) *pi/180;
SolInfo(:,3) = -Solmpc.bus(:,PD)/Sb1;
SolInfo(:,4) = -Solmpc.bus(:,QD)/Sb1;

for i=(1:Ngen)
    idGen = Solmpc.gen(i,BUS_I);
    idBus = find(Solmpc.bus(:,1)==idGen);
    SolInfo(idBus,3) = SolInfo(idBus,3) + Solmpc.gen(i,PG)/Sb1;
    SolInfo(idBus,4) = SolInfo(idBus,4) + Solmpc.gen(i,QG)/Sb1;
end

SolInfo3 = zeros(NBus, 4);% VM, VA, P, Q

SolInfo3(:,1) = Soltest2.bus(:,VM);
SolInfo3(:,2) = Soltest2.bus(:,VA) *pi/180;
SolInfo3(:,3) = -Soltest2.bus(:,PD)/Sb1;
SolInfo3(:,4) = -Soltest2.bus(:,QD)/Sb1;

for i=(1:Ngen)
    idGen = Soltest2.gen(i,BUS_I);
    idBus = find(Soltest2.bus(:,1)==idGen);
    SolInfo3(idBus,3) = SolInfo3(idBus,3) + Soltest2.gen(i,PG)/Sb1;
    SolInfo3(idBus,4) = SolInfo3(idBus,4) + Soltest2.gen(i,QG)/Sb1;
end




%% 

writematrix(AgentInfo, name1, 'Delimiter',' ');
writematrix(lineInfo, name2, 'Delimiter',' ');
writematrix(CaseInfo, name3, 'Delimiter',' ');
writematrix(ConstraintInfo, name4,'Delimiter',' ');
writematrix(SolInfo3, name5,'Delimiter',' ');
%writematrix(Bgrid1,name6,'Delimiter',' ');
%writematrix(Ggrid1,name7,'Delimiter',' ');

%%

test.bus(:,BUS_TYPE) = PQ; % que des noeuds PQ
test.bus(1,BUS_TYPE) = REF;

test.branch(:,F_BUS:T_BUS)= lineInfo(:,1:2);

test.bus(:,BUS_I) = coresBus;

offset = nCons+nGenSup;
for n=(1:offset)
    idBus = AgentInfo(n,1);
    Pd = -AgentInfo(n,4);
    Qd = -AgentInfo(n,7);
    test.bus(idBus,PD) = Pd;
    test.bus(idBus,QD) = Qd;
end
%test.bus(:,VM)= SolInfo(1,1);
%test.bus(:,VA) = SolInfo(1,2)*180/pi;





for n=(1:Ngen)
    idBus = AgentInfo(offset+n,1);
    Pg = AgentInfo(offset+n,4);
    Qg = AgentInfo(offset+n,7);
    test.bus(idBus,PD) = test.bus(idBus,PD) - Pg;
    test.bus(idBus,QD) = test.bus(idBus,QD) - Qg;
end
test.gen=mpc.gen(idGenRef,:);
test.gen(1,PG) = 0;
test.gen(1,QG) = 0;
try
    test.gencost =mpc.gencost(idGenRef,:);
catch
    
end


atest = runpf(test,mpopt);
err1 = (SolInfo(:,3) + atest.bus(:,PD))';
err1(1) = err1(1) - atest.gen(1,PG);
err2 = (SolInfo(:,4) + atest.bus(:,QD))';
err2(1) = err2(1) - atest.gen(1,QG);

SolInfo2 = zeros(NBus, 4);% VM, VA, P, Q

SolInfo2(:,1) = atest.bus(:,VM);
SolInfo2(:,2) = atest.bus(:,VA) *pi/180;
SolInfo2(:,3) = -atest.bus(:,PD)/Sb1;
SolInfo2(:,4) = -atest.bus(:,QD)/Sb1;

SolInfo2(1,3) = SolInfo2(1,3) + atest.gen(1,PG)/Sb1;
SolInfo2(1,4) = SolInfo2(1,4) + atest.gen(1,QG)/Sb1;


sum(SolInfo2-SolInfo, "all")
