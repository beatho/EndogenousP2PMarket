%% compare impedance

mpc = case14();

[Ybus, Yf, Yt] = makeYbus(ext2int(mpc));

G = real(Ybus);
B = imag(Ybus);

B2 = full(B);


Ybus2 = full(Ybus);
%Yperso = zeros(NBus,NBus);
Shunt = diag(ConstraintInfo(:,1)) + 1j * diag(ConstraintInfo(:,2));
Yperso = Shunt/Sb1;

for l=(1:nLine)
    i = coresBus(mpc.bus(:,1) == mpc.branch(l,F_BUS));
    j = coresBus(mpc.bus(:,1) == mpc.branch(l,T_BUS));
    
    r = mpc.branch(l,BR_R);
    x = mpc.branch(l,BR_X);
    b = mpc.branch(l,BR_B);
    tau = mpc.branch(l,TAP);
    theta = mpc.branch(l,SHIFT)*pi/180; 
    
    z = r + 1i*x;
    y = 1/z;
    if(tau>0)
        Yperso(i,i) =  Yperso(i,i) + (y + b/2 *1j)/(tau*tau); % Yff
        Yperso(j,j) =  Yperso(j,j) + (y + b/2 *1j); % Ytt
        Yperso(i,j) =  Yperso(i,j) - y / (exp(-1j*theta) *tau );
        Yperso(j,i) =  Yperso(j,i) - y / (exp(1j*theta) *tau );
        %Yperso(i,j) =  Yperso(i,j) - y*(cos(theta) + 1j*sin(theta)) / (tau ); % Yft
        %Yperso(j,i) =  Yperso(j,i) - y*(cos(theta) - 1j*sin(theta)) / (tau ); % Ytf
    else
        Yperso(i,i) =  Yperso(i,i) + (y+b/2 *1j); % Yff
        Yperso(j,j) =  Yperso(j,j) + (y+b/2 *1j); % Ytt
    
        Yperso(i,j) =  Yperso(i,j)-y; % Yft
        Yperso(j,i) =  Yperso(j,i)-y; % Ytf
    end
end

sum(Yperso-Ybus, 'all')
[i j] = find((Yperso-Ybus2).*conj(Yperso-Ybus2)>0.00001);
n = size(i);
for k=(1:n)
    buses = [mpc.bus(i(k),1) mpc.bus(j(k),1) ];
    buses2 = [mpc.bus(j(k),1) mpc.bus(i(k),1) ]
    branchTest = find(mpc.branch(:,1:2)==buses2);
    [Yperso(i(k),j(k)) Ybus2(i(k),j(k))];
    
end


