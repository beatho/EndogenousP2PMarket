define_constants;
mpc = case4_perso();
Ybus = makeYbus(ext2int(mpc));
Ybus2 = full(YbusPerso2);

YbusPerso = zeros(4,4);
inverse=false;
for l=(1:3)
    i = mpc.branch(l,F_BUS);
    j = mpc.branch(l,T_BUS);
    r = mpc.branch(l,BR_R);
    x = mpc.branch(l,BR_X);
    b = mpc.branch(l,BR_B);

    tau = mpc.branch(l,TAP); 
    theta = mpc.branch(l,SHIFT);
    
    z = r + 1i*x;
    y = 1/z;
    if(tau >0) 
        if(i>j)
            t = i;
            i = j;
            j = t;
            inverse = true;
        end
        Yij = -y*exp(1i*theta)/tau;
        Yji = -y*exp(-1i*theta)/tau;
        if(inverse)
            Yjj = y /(tau*tau);
            Yii = y;
        else
            Yii = y /(tau*tau);
            Yjj = y;
        end
        
    else
        Yij = -y;
        Yji = -y;
        Yii = y + b/2;
        Yjj = y + b/2;
    end
    YbusPerso(i,j) = YbusPerso(i,j) + Yij;
    YbusPerso(j,i) = YbusPerso(j,i) + Yji;
    YbusPerso(i,i) = YbusPerso(i,i) + Yii;
    YbusPerso(j,j) = YbusPerso(j,j) + Yjj;
    inverse=false;
end

dY = YbusPerso -Ybus2
runopf(case4_dist)
