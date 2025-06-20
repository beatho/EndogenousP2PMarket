using PowerModels, ExaModels, NLPModelsIpopt, CUDA
using MadNLPGPU, KernelAbstractions

function termination_code(status::MadNLP.Status)
    if status == MadNLP.SOLVE_SUCCEEDED
        return 1
    elseif status == MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL
        return 2
    elseif status == MadNLP.DIVERGING_ITERATES || status == MadNLP.DIVERGING_ITERATES
        return 3
    else
        return 4
    end
end

function parse_ac_power_data(filename)
    data = PowerModels.parse_file(filename)
    flexible = true;
    
    PowerModels.standardize_cost_terms!(data, order = 2)
    #PowerModels.print_summary(data)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    # arcs : (l,i,j) dans les deux sens, arcs_from et arc_to pour dans un seul sens
    
    #=println(ref[:bus])
    println(ref[:ref_buses])
    println(typeof(ref[:arcs]))
    println(ref[:arcs]);
    println("*******************************************************");
    println(keys(data))
    for key in keys(ref)
        println(key);
        println(ref[key]);
         println("------------------------");
    end
    println("*******************************************************");
    
    for (i, branch) in ref[:branch]
        println(branch);
        println(PowerModels.calc_branch_y(branch))
    end=#

    S0 = data["baseMVA"];
    
    B = length(data["bus"]);
    L = length(data["branch"]);
    if((B-1) != L)
        println("grid not radial");
        exit(-1);
    end

    temp = ref[:arcs];
    ref[:arcs] = Vector{Tuple{Int64, Int64, Int64}}(undef, 0);

    for (k, (l, i, j)) in enumerate(temp)
        if(i<j)
            push!(ref[:arcs], (l,i,j));
        end
    end
    
    #println(ref[:arcs]);

    linePreB = Vector{Int64}(undef, B); 
    linePreLine = Vector{Int64}(undef, L);
    fromb =  Vector{Int64}(undef, L);
    to_b  =  Vector{Int64}(undef, L);

    for (k, (l, i, j)) in enumerate(ref[:arcs])
        fromb[l] = i;
        to = j;
        to_b[l] = j;
        linePreB[to] = i;
    end

    #println(linePreB)

    for l=2:L
        linePreLine[l] = linePreB[fromb[l]];
    end
    linePreLine[1] = 0;
    linePreB[1] = 0;

    ## recuperation de tous les agents
    N_c = length(data["load"]);
    N_g = length(data["gen"]);
    N = N_c + N_g;
    N_gnc = 0;

    a = Vector{Float64}(undef, 2*N);
    b = Vector{Float64}(undef, 2*N);
    pstart = Vector{Float64}(undef, 2*N);
    pmin = Vector{Float64}(undef, 2*N);
    pmax = Vector{Float64}(undef, 2*N);
    bus = Vector{Int64}(undef, N);
    AgentByBus = Vector{Vector{Int64}}(undef, B);
    lineAgent = Vector{Int64}(undef, N);
    oldIndice = Vector{Int64}(undef, N);

    for i=(1:B)
        AgentByBus[i] = Vector{Float64}(undef, 0);
    end

    for load in ref[:load]
        truc = load[2]
        Pobj = - truc["pd"];
        if Pobj > 0 # gen 
            N_gnc += 1;
        end
    end
    offset = 0
    #Nc -= N_gnc;
    indice = N_gnc - 1;

    for k = (1:N_c)
        load = ref[:load][k]
        truc = load;
        i = Int64(truc["index"])
        #println("k : $(k) ", i, " ", truc["pd"])
        Pobj = - truc["pd"];
        if Pobj>0
            oldIndice[i] = N_c - indice;
            indice -= 1;
            offset += 1;
        elseif Pobj<0
            oldIndice[i] = i - offset;
        end
        a[i] = 0.1 * S0*S0;
        b[i] = - Pobj * a[i];
        pstart[i] = Pobj;
        if flexible
            pmin[i] = Pobj < 0 ? 1.2 * Pobj : 0.8 * Pobj;
            pmax[i] = Pobj < 0 ? 0.8 * Pobj : 1.2 * Pobj;
        else
            pmin[i] = Pobj;
            pmax[i] = Pobj;
        end


        bus[i] = truc["load_bus"];
        push!(AgentByBus[bus[i]], i );
        lineAgent[i] = linePreB[bus[i]];

        Qobj = - truc["qd"];
        a[i + N] = 0.1 * S0 * S0;
        b[i + N] = - Qobj * a[i];
        pstart[i + N] = Qobj;
        if flexible
            pmin[i + N] = Qobj < 0 ? 1.05 * Qobj : 0.95 * Qobj;
            pmax[i + N] = Qobj < 0 ? 0.95 * Qobj : 1.05 * Qobj;
        else
            pmin[i + N] = Qobj;
            pmax[i + N] = Qobj;
        end
        
    end

  

    for i = 1:N_g
        a[N_c + i]      = 0.1 * S0 * S0;
        b[N_c + i]      = ref[:gen][i]["cost"][2];
        pstart[N_c + i] = ref[:gen][i]["pg"];
        pmin[N_c + i]   = max(ref[:gen][i]["pmin"], 0);
        pmax[N_c + i]   = ref[:gen][i]["pmax"];

        oldIndice[N_c + i] = N_c + i;

        bus[N_c + i] = ref[:gen][i]["gen_bus"];
        push!(AgentByBus[bus[N_c + i]], N_c + i );
        lineAgent[N_c + i] = linePreB[bus[N_c + i]];

        a[N_c + i + N]      = 0.1 * S0 * S0;
        b[N_c + i + N]      = - ref[:gen][i]["qg"] * a[N_c + i + N] ;
        pstart[N_c + i + N] =   ref[:gen][i]["qg"];
        pmin[N_c + i + N]   =   ref[:gen][i]["qmin"];
        pmax[N_c + i + N]   =   ref[:gen][i]["qmax"];
    end
    ## creation 
    ref[:gen] = Dict{Int64, Any}();

    for n in (1:N)
        #println("indice : ", n, " type ", typeof(n));
        ref[:gen][n] = Dict{String, Any}();
        ref[:gen][n]["cost"]     = [a[n] b[n] 0];
        ref[:gen][n]["gen_bus"]  = bus[n];
        ref[:gen][n]["line_pre"] = lineAgent[n]; 
        ref[:gen][n]["pmin"]     = pmin[n];
        ref[:gen][n]["pmax"]     = pmax[n];
        ref[:gen][n]["pstart"]   = pstart[n];
        ref[:gen][n]["costg"]    = [a[n + N] b[n + N] 0];
        ref[:gen][n]["qmin"]     = pmin[n + N];
        ref[:gen][n]["qmax"]     = pmax[n + N];
        ref[:gen][n]["qstart"]   = pstart[n + N];
        ref[:gen][n]["oldi"]     = oldIndice[n];
    end

    arcdict = Dict(a => k for (k, a) in enumerate(ref[:arcs])) # (id_branche, bus1, bus2) => idCores et (id_branche, bus2, bus1) => idCores
    busdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:bus]))
    gendict = Dict(k => i for (i, (k, v)) in enumerate(ref[:gen]))
    branchdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:branch]))

    #println(arcdict)
    #println(gendict)

        
    return (
        bus = [
            begin
                bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][k]]
                gs = sum(shunt["gs"] for shunt in bus_shunts; init = 0.0)
                bs = sum(shunt["bs"] for shunt in bus_shunts; init = 0.0)
                (i = busdict[k],
                gs = gs,
                bs = bs
                #vmax   = v["vmax"],
                #vmin   = v["vmin"] 
                )
            end for (k, v) in ref[:bus]
        ],
        agent = [
            (
                i = gendict[k], #gendict[k]
                oldi = v["oldi"],
                cost1 = v["cost"][1],
                cost2 = v["cost"][2],
                cost1q = v["costg"][1],
                cost2q = v["costg"][2],
                #=pmax   = v["pmax"] ,
                pmin   = v["pmin"] ,
                qmax   = v["qmax"] ,
                qmin   = v["qmin"] ,
                pstart = v["pstart"],
                qstart = v["qstart"],=#
                bus = busdict[v["gen_bus"]], # busdict[v["gen_bus"]]
                line = v["line_pre"]>0 ? branchdict[v["line_pre"]] : 0  #branchdict[v["line_pre"]]
            ) for (k, v) in ref[:gen]
        ],
        arc = [
            (
                i = k,
                rate_a = ref[:branch][l]["rate_a"],
                bus_f = busdict[i],
                bus_t = busdict[j],
            ) for (k, (l, i, j)) in enumerate(ref[:arcs])
        ],
        branch = [
            begin
                f_idx = arcdict[i, branch["f_bus"], branch["t_bus"]] # identification f_idx, t_idx sont les indices dans arcdict selon ordre des bus
                #t_idx = arcdict[i, branch["t_bus"], branch["f_bus"]]
                (
                    i = branchdict[i], #branchdict
                    #j = 1, # à quoi cela sert ?????????????
                    f_idx = f_idx, 
                    #t_idx = t_idx,
                    bus_f = busdict[branch["f_bus"]], # busdict
                    bus_t = busdict[branch["t_bus"]], # busdict
                    x = branch["br_x"],
                    r = branch["br_r"],
                    rate_a_sq = branch["rate_a"]^2,
                    linePre = linePreLine[i] > 0 ? branchdict[linePreLine[i]] : 0 #branchdict[linePreLine[i]]
                )
            end for (i, branch) in ref[:branch]
        ],
        ref_buses = [busdict[i] for (i, k) in ref[:ref_buses]],
        vmax = [v["vmax"] for (k, v) in ref[:bus]],
        vmin = [v["vmin"] for (k, v) in ref[:bus]],
        pmax = [v["pmax"] for (k, v) in ref[:gen]],
        pmin = [v["pmin"] for (k, v) in ref[:gen]],
        qmax = [v["qmax"] for (k, v) in ref[:gen]],
        qmin = [v["qmin"] for (k, v) in ref[:gen]],
        pstart = [v["pstart"] for (k, v) in ref[:gen]],
        qstart = [v["qstart"] for (k, v) in ref[:gen]],
        rate_a = [ref[:branch][l]["rate_a"] for (k, (l, i, j)) in enumerate(ref[:arcs])],
        angmax = [b["angmax"] for (i, b) in ref[:branch]],
        angmin = [b["angmin"] for (i, b) in ref[:branch]],
    )
end

function lecture_file(filename, size)
    fileMat = Vector{Vector{Float64}}(undef, size)
    k = 1;
    for i=(1:size)
        fileMat[i] = Vector{Float64}(undef, 0);
    end 
    open(filename) do file
        for l in eachline(file)
            test = split(l, " ")
            for t in test
                push!(fileMat[k], parse(Float64, t))
            end
            k+=1;
            if k==size+1
                break;
            end
        end
    end
    return fileMat;
end

function save_file(filename, Mat)
    open(filename, "a") do io
        for a in Mat
            for b in a
                write(io, "$(b); ")
            end
            write(io, "\n")
        end
    end

end

function Data_from_AC(path, name, sizeRandom = nothing)
    fileName1 = path * "Case"   * name * ".txt";
	fileName2 = path * "Agent"  * name * ".txt";
	fileName3 = path * "Branch" * name * ".txt";
	fileName4 = path * "Bus"	* name * ".txt";
    fileName5 = path * "AgentConsumption" * name * ".txt";

    ref = Dict{Symbol, Any}();

    Info = lecture_file(fileName1, 1); # Sbase, Vbase, nAgent, nCons, nGenSup, nBus, nLine, V0, theta0
	_Sbase = Info[1][1];
    _Vbase = Info[1][2];
    _Zbase = _Vbase * _Vbase / _Sbase;
	_V0 = Info[1][8];
	_theta0 = Info[1][9];
	_AC = true;
    step = 1;

	N =  round(Int64, Info[1][3]); # + the loss agent and grid agent for TestFeeder
    B =  round(Int64, Info[1][6]);
    L =  round(Int64, Info[1][7]);
    
    Mat = lecture_file(fileName2, N);#    (_nAgent - 1, 10); # bus, a, b, P, Pmin, Pmax, Qobj, Qmin, Qmax, zone
    MatLine = lecture_file(fileName3, L); # from, to, Ys Real, Ys Im, Yp , tau, theta, Limit=0, zs Real, Zs imag
	MatBus = lecture_file(fileName4, B); # Gs, Bs, min, max, V0, theta0, zone � rajouter ici ?
    MatConso = [[0] for i=(1:N)];
    
    try
        test = lecture_file(fileName5, N);
        MatConso = test;
    catch
        #nothing to do
    end
   
    
    offsetbus = 0;
    #inversionLine = Vector{Int64}(undef, L);
	if (B == L + 1) 
		for l = 1:L 
			if (MatLine[l][1] > MatLine[l][2]) 
				temp = round(Int64, MatLine[l][1]);
				MatLine[l][1] =  MatLine[l][2];
				MatLine[l][2] = temp;
				MatLine[l][7] = -MatLine[l][7];	# ????
            end
            if (MatLine[l][1] == 0 || MatLine[l][2] == 0)
                offsetbus = 1;
            end
        end
		# il faut ordonner pour que Matline(k,1) = k+1;
		radial = true;
	else 
		radial = false;
        for l = 1:L 
            if (MatLine[l][1] == 0 || MatLine[l][2] == 0)
                offsetbus = 1;
            end
        end
    end
	
	ref[:branch] = Dict{Int64, Any}();
    ref[:arcs]   = Vector{Tuple{Int64, Int64, Int64}}(undef, L);
   	
	for l = (1:L) 
        ref[:branch][l] = Dict{String, Any}();
		i = MatLine[l][1] + offsetbus;
		j = MatLine[l][2] + offsetbus;
		
		limit = MatLine[l][8] > 0 ? MatLine[l][8] : 10000;
		ZsRe = MatLine[l][9];
		ZsIm = MatLine[l][10];
		YlsRe = MatLine[l][3]; # re(1/Z)
        YlsIm = MatLine[l][4]; # im(1/Z)
        Ylp = MatLine[l][5]; # b/2
        tau = MatLine[l][6] > 0 ? MatLine[l][6] : 1;
        theta = MatLine[l][7];

        ref[:branch][l]["f_bus"] = i;
        ref[:branch][l]["t_bus"] = j;
        ref[:branch][l]["rate_a"] = limit;
        ref[:branch][l]["br_r"] = ZsRe;
        ref[:branch][l]["br_x"] = ZsIm;
        ref[:branch][l]["br_g"] = YlsRe;
        ref[:branch][l]["br_b"] = YlsIm;
        ref[:branch][l]["br_c"] = Ylp;
        ref[:branch][l]["tap"] = tau;
        ref[:branch][l]["shift"] = theta;

        ref[:arcs][l] = (l, i, j);
       
    end
	
    ref[:bus] = Dict{Int64, Any}();
	for i = 1:B # bound on voltage angle rad
        ref[:bus][i] = Dict{String, Any}();
        ref[:bus][i]["angmax"] = 3;
        ref[:bus][i]["angmin"] = -3;
        ref[:bus][i]["angtart"] = MatBus[i][6];
        ref[:bus][i]["vstart"]  = MatBus[i][5];

		ref[:bus][i]["gs"] = MatBus[i][1] / _Zbase;
        ref[:bus][i]["bs"] = MatBus[i][2] / _Zbase;

        ref[:bus][i]["vmin"] = MatBus[i][3];
        ref[:bus][i]["vmax"] = MatBus[i][4];

    end

    ref[:ref_buses] = Dict{Int64, Any}();
    ref[:ref_buses][1] = ref[:bus][1];
    
    linePreB = Vector{Int64}(undef, B); 
    linePreLine = Vector{Int64}(undef, L);
    fromb =  Vector{Int64}(undef, L);
    to_b  =  Vector{Int64}(undef, L);

    for (k, (l, i, j)) in enumerate(ref[:arcs])
        fromb[l] = i;
        to = j;
        to_b[l] = j;
        linePreB[to] = i;
    end
    for l=2:L
        linePreLine[l] = linePreB[fromb[l]];
    end
    linePreLine[1] = 0;
    linePreB[1] = 0;
    Nc, Ng, Np = (0, 0, 0)

    if isnothing(sizeRandom) 
        if name == "TestFeeder"
            ref, Nc, Ng, Np = Agent_From_TestFeeder(ref, Mat, MatConso, Info, linePreB);
        else
            ref, Nc, Ng, Np = Agent_from_File(ref, Mat, Info, linePreB);
        end
    else
        ref, Nc, Ng, Np = Agent_from_random(ref, sizeRandom, B, linePreB)
    end

    for key in keys(ref)
        ref[key]= sort(ref[key], by =  x->x[1])
    end

    #=ref[:branch] = sort(ref[:branch], by =  x->x[1])
    ref[:bus] = sort(ref[:bus], by =  x->x[1])
    ref[:arcs] = sort(ref[:arcs], by =  x->x[1])
    ref[:gen] = sort(ref[:gen], by =  x->x[1])=#


    arcdict = Dict(a => k for (k, a) in enumerate(ref[:arcs])) # (id_branche, bus1, bus2) => idCores et (id_branche, bus2, bus1) => idCores
    busdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:bus]))
    gendict = Dict(k => i for (i, (k, v)) in enumerate(ref[:gen]))
    branchdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:branch]))

    #for i = (1:N+2)
        #println(ref[:gen][i]["pstart"], " ", ref[:gen][i]["qstart"] )
    #end

    for a in sort(ref[:gen], by = x-> x[1])
        #println(a)
        #println(gendict[a[1]], " ", a[2]["oldi"], " ", a[2]["gen_bus"], " ",  a[2]["pstart"]);
    end


	return (
        bus = [
            (
                i = busdict[k],
                gs = v["gs"],
                bs = v["bs"]
                #vmax   = v["vmax"],
                #vmin   = v["vmin"] 
                )
             for (k, v) in ref[:bus]
        ],
        arc = [
        (
            i = k,
            rate_a = ref[:branch][l]["rate_a"],
            bus_f = busdict[i],
            bus_t = busdict[j],
        ) for (k, (l, i, j)) in enumerate(ref[:arcs])
        ],
        agent = [
            (
                i = gendict[k], 
                oldi = v["oldi"],
                cost1 = v["cost"][1],
                cost2 = v["cost"][2],
                cost1q = v["costg"][1],
                cost2q = v["costg"][2],
                lb = v["lb"],
                ub = v["ub"],
                lbq = v["lbq"],
                ubq = v["ubq"],
                Mn = v["nVoisin"],
                gamma = v["gamma"],
                bus = busdict[v["gen_bus"]], # busdict[v["gen_bus"]]
                line = v["line_pre"]>0 ? branchdict[v["line_pre"]] : 0  #branchdict[v["line_pre"]]
            ) for (k, v) in ref[:gen]
        ],
        branch = [
            begin
                f_idx = arcdict[i, branch["f_bus"], branch["t_bus"]] # identification f_idx, t_idx sont les indices dans arcdict selon ordre des bus
                g = branch["br_g"]
                b = branch["br_b"]
                tr = branch["tap"] .* cos.(branch["shift"])
                ti = branch["tap"] .* sin.(branch["shift"])
                ttm = tr^2 + ti^2
                g_fr = 0 # je ne sais pas ce que c'est mais cela vaut 0 ? Peut être transfo
                b_fr = 0
                g_to = 0
                b_to = 0
                c1 = (-g * tr - b * ti) / ttm
                c2 = (-b * tr + g * ti) / ttm
                c3 = (-g * tr + b * ti) / ttm
                c4 = (-b * tr - g * ti) / ttm
                c5 = (g + g_fr) / ttm
                c6 = (b + b_fr) / ttm
                c7 = (g + g_to)
                c8 = (b + b_to)
                (
                    i = branchdict[i], #branchdict
                    f_idx = f_idx, 
                    bus_f = busdict[branch["f_bus"]], # busdict
                    bus_t = busdict[branch["t_bus"]], # busdict
                    x = branch["br_x"],
                    r = branch["br_r"],
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    c4 = c4,
                    c5 = c5,
                    c6 = c6,
                    c7 = c7,
                    c8 = c8,
                    rate_a_sq = branch["rate_a"]^2,
                    linePre = linePreLine[i] > 0 ? branchdict[linePreLine[i]] : 0 #branchdict[linePreLine[i]]
                )
            end for (i, branch) in ref[:branch]
        ],
        ref_buses = [busdict[i] for (i, k) in ref[:ref_buses]],
        vmax = [v["vmax"] for (k, v) in ref[:bus]],
        vmin = [v["vmin"] for (k, v) in ref[:bus]],
        pmax = [v["pmax"] for (k, v) in ref[:gen]],
        pmin = [v["pmin"] for (k, v) in ref[:gen]],
        qmax = [v["qmax"] for (k, v) in ref[:gen]],
        qmin = [v["qmin"] for (k, v) in ref[:gen]],
        pstart = [v["pstart"] for (k, v) in ref[:gen]],
        qstart = [v["qstart"] for (k, v) in ref[:gen]],
        rate_a = [ref[:branch][l]["rate_a"] for (k, (l, i, j)) in enumerate(ref[:arcs])],
        angmax = [b["angmax"] for (i, b) in ref[:bus]],
        angmin = [b["angmin"] for (i, b) in ref[:bus]],
        Nc = Nc, 
        Ng = Ng,
        Np = Np,
        ConsoAll = MatConso,
        MatAgent = Mat,
        step = 1,
    )
end

function Agent_from_File(ref, Mat, Info, linePreB)
    _Sbase = Info[1][1];
    N =  round(Int64, Info[1][3]) + 1; # + the loss agent
    nCons  =  round(Int64, Info[1][4]) + 1; # + the loss agent
	nGen = N - nCons;
	nPro = 0;

    a = Vector{Float64}(undef, 2*N);
    b = Vector{Float64}(undef, 2*N);
    pstart = Vector{Float64}(undef, 2*N);
    pmin = Vector{Float64}(undef, 2*N);
    pmax = Vector{Float64}(undef, 2*N);
    Lb = Vector{Float64}(undef, 2*N);
    Ub = Vector{Float64}(undef, 2*N);
    nVoisin = Vector{Float64}(undef, 2*N);
    bus = Vector{Int64}(undef, N);
    lineAgent = Vector{Int64}(undef, N);
    oldIndice = Vector{Int64}(undef, N);

  
    offsetbus = 0;
    for i = 1:N-1
        if (Mat[i][1] == 0)
            offsetbus = 1;
            break
        end
    end

    oldIndice[1] = 1;
    bus[1] = 1;
    lineAgent[1] = linePreB[bus[1]];
	Lb[1] =  -10000; # pour ne pas avoir besoin de le modifier
	Ub[1] =  0;
    Ub[N + 1] =   10000; # idem
    Lb[N + 1] =  -10000;
    a[1] = 0;
    b[1] = 0;
    pstart[1] = 0;
    pmax[1] =  0.00001;
    pmin[1] = -0.00001;
    pstart[N + 1] = 0;
    pmax[N + 1] =  0.000001;
    pmin[N + 1] = -0.000001;
    a[N+1] = 0;
    b[1] = 0;
	nVoisin[1] = nGen;
	nVoisin[N + 1] = N - 1;
	for i = 2:N
        oldIndice[i] = i;
        bus[i] =  Mat[i - 1][1] + offsetbus;
        lineAgent[i] = linePreB[bus[i]];
		cost1 = Mat[i - 1][2] * (_Sbase * _Sbase);
		cost2 = Mat[i - 1][3] * _Sbase;
		pstart[i] = Mat[i - 1][4] / _Sbase;
		pLim1 = Mat[i - 1][5] / _Sbase;
		pLim2 = Mat[i - 1][6] / _Sbase;
		pstart[N + i]  = Mat[i - 1][7] / _Sbase;

		qLim1 = Mat[ i - 1 ][8] / _Sbase;
		qLim2 = Mat[ i - 1 ][9] / _Sbase;
		
		costQ1 = 0.1 * (_Sbase * _Sbase);
		costQ2 = -costQ1 * pstart[N + i];
		
		if i <= nCons 
			Mn = nGen;
			Ub[ i ] = 0;
			Lb[ i ] = pLim1;
		else 
			Mn = nCons;
			Ub[ i ] = pLim2;
			Lb[ i ] = 0;
        end
		pmin[ i ] = pLim1;
		pmax[ i ] = pLim2;


		a[ i ] = cost1;
		b[ i ] = cost2;
		a[ i + N ] = costQ1;
		b[ i + N ] = costQ2;

		
		pmin[ i + N ] = qLim1;
		pmax[ i + N ] = qLim2;
		nVoisin[ i ] = Mn;
		nVoisin[ i + N ] = N - 1;
		Ub[ i + N ] = qLim2 * (qLim2 > 0);
		Lb[ i + N ] = qLim1 * (qLim1 < 0);
    end

    ref[:gen] = Dict{Int64, Any}();
    for n in (1:N)
        ref[:gen][n] = Dict{String, Any}();
        ref[:gen][n]["cost"]     = [a[n] b[n] 0];
        ref[:gen][n]["gen_bus"]  = bus[n];
        ref[:gen][n]["line_pre"] = lineAgent[n]; 
        ref[:gen][n]["pmin"]     = pmin[n];
        ref[:gen][n]["pmax"]     = pmax[n];
        ref[:gen][n]["pstart"]   = pstart[n];
        ref[:gen][n]["costg"]    = [a[n + N] b[n + N] 0];
        ref[:gen][n]["qmin"]     = pmin[n + N];
        ref[:gen][n]["qmax"]     = pmax[n + N];
        ref[:gen][n]["qstart"]   = pstart[n + N];
        ref[:gen][n]["oldi"]     = oldIndice[n];
        ref[:gen][n]["gamma"]    = 0;

        ref[:gen][n]["nVoisin"]  = nVoisin[n];
        ref[:gen][n]["lb"]       = Lb[n];
        ref[:gen][n]["ub"]       = Ub[n];
        ref[:gen][n]["lbq"]      = Lb[n + N];
        ref[:gen][n]["ubq"]      = Ub[n + N];
    end

    return ref, nCons, nGen, nPro
end

function Agent_From_TestFeeder(ref, MatAgent, ConsoAll, Info, linePreB)
    _Sbase = Info[1][1];
    N =  round(Int64, Info[1][3]) + 2; # + the loss agent
    nCons  =  round(Int64, Info[1][4]) + 1; # + the loss agent
	nGen = 0;
	nPro = 1;
    dP = 0.1; # 10%

    a = Vector{Float64}(undef, 2*N);
    b = Vector{Float64}(undef, 2*N);
    
    pstart = Vector{Float64}(undef, 2*N);
    pmin = Vector{Float64}(undef, 2*N);
    pmax = Vector{Float64}(undef, 2*N);
    Lb = Vector{Float64}(undef, 2*N);
    Ub = Vector{Float64}(undef, 2*N);
    nVoisin = Vector{Float64}(undef, 2*N);
    
    lineAgent = Vector{Int64}(undef, N);
    oldIndice = Vector{Int64}(undef, N);
   
    
    bus = Vector{Int64}(undef, N);
    offsetbus = 0;
    for i = 1:N-2
        if (MatAgent[i][1] == 0)
            offsetbus = 1;
            break
        end
    end
   
    # loss agent
    oldIndice[1] = 1;
    bus[1] = 1;
    lineAgent[1] = linePreB[bus[1]];
	Lb[1] =  -10000; # pour ne pas avoir besoin de le modifier
    Ub[1] =  0;
	Ub[N + 1] =   10000; # idem
    Lb[N + 1] =  -10000;
    a[1] = 0;
    b[1] = 0;
    pstart[1] = 0;
    pmax[1] = 0;
    pmin[1] = 0;
    pstart[N + 1] = 0;
    pmax[N + 1] = 0;
    pmin[N + 1] = 0;
	nVoisin[1] = nGen + nPro;
	nVoisin[N + 1] = N - 1;
    a[N + 1] = 0;
    b[N + 1] = 0;

    # grid agent
    oldIndice[N] = N;
    bus[N] = 1;
    lineAgent[N] = linePreB[bus[N]];
	Lb[N] =  -10000; # pour ne pas avoir besoin de le modifier
    Ub[N] =  10000;
	Ub[2*N] =   10000; # idem
    Lb[2*N] =  -10000;
    a[N] = 0.01;
    b[N] = 0.03;
    pstart[N] = 0;
    pmax[N] = 10000;
    pmin[N] = -10000;
    pstart[2*N] = 0;
    pmax[2*N] = 10000;
    pmin[2*N] = -10000;
	nVoisin[N] = nCons;
	nVoisin[2*N] = N - 1;
    a[2*N] = 0;
    b[2*N] = 0;

	for i = 2:N-1
        bus[i] =  MatAgent[i - 1][1] + offsetbus;
        lineAgent[i] = linePreB[bus[i]];
		factor = MatAgent[i - 1][2];
		PF =  MatAgent[i - 1][3]
        Pobj = ConsoAll[i - 1][1];
        oldIndice[i] = i;
        
		
        S0 = factor* Pobj;
        P0 = S0 * PF;
        pstart[i] = - P0;

        cost1 = 1;
		cost2 = P0 * cost1;

		pLim1 = -(1 + dP) * P0;
		pLim2 = -(1 - dP) * P0;

        Q0 = S0 * sqrt(1 - PF * PF );
        signe = 2 * ((i-1) % 5) - 1; # un sur 5 ?
		Q0 = signe * Q0; 
        
        pstart[N + i] = Q0;

        qLim1 = Q0 * (1 - dP + 2 * dP * (Q0 < 0));
        qLim2 = Q0 * (1 + dP - 2 * dP * (Q0 < 0));

		costQ1 = 1 * (_Sbase * _Sbase);
		costQ2 = -costQ1 * Q0;
		
		 
        nVoisin[ i ] = nGen + nPro;
        Ub[ i ] = 0;
        Lb[ i ] = pLim1;
		
		pmin[ i ] = pLim1;
		pmax[ i ] = pLim2;

		a[ i ] = cost1;
		b[ i ] = cost2;

		a[ i + N ] = costQ1;
		b[ i + N ] = costQ2;
		
		pmin[ i + N ] = qLim1;
		pmax[ i + N ] = qLim2;
		
		nVoisin[ i + N ] = N - 1;
		Ub[ i + N ] = qLim2 * (qLim2 > 0);
		Lb[ i + N ] = qLim1 * (qLim1 < 0);
    end

    ref[:gen] = Dict{Int64, Any}();
    for n in (1:N)
        ref[:gen][n] = Dict{String, Any}();
        ref[:gen][n]["cost"]     = [a[n] b[n] 0];
        ref[:gen][n]["gen_bus"]  = bus[n];
        ref[:gen][n]["line_pre"] = lineAgent[n]; 
        ref[:gen][n]["pmin"]     = pmin[n];
        ref[:gen][n]["pmax"]     = pmax[n];
        ref[:gen][n]["pstart"]   = pstart[n];
        ref[:gen][n]["costg"]    = [a[n + N] b[n + N] 0];
        ref[:gen][n]["qmin"]     = pmin[n + N];
        ref[:gen][n]["qmax"]     = pmax[n + N];
        ref[:gen][n]["qstart"]   = pstart[n + N];
        ref[:gen][n]["oldi"]     = oldIndice[n];
        ref[:gen][n]["gamma"]    = 0;

        ref[:gen][n]["nVoisin"]  = nVoisin[n];
        ref[:gen][n]["lb"]       = Lb[n];
        ref[:gen][n]["ub"]       = Ub[n];
        ref[:gen][n]["lbq"]      = Lb[n + N];
        ref[:gen][n]["ubq"]      = Ub[n + N];
    end

    return ref, nCons, nGen, nPro
end

function Data_from_random(nAgent, nBus, nDeep, nBranch)
    ref = Dict{Symbol, Any}();
    BTLINE = [0.306, 0.29, 13.2, 0.35]; # 94-AL1/15-ST1A 0.4 : r x b Imax (ohm, ohm, nF) /km kA
    billion = 1e9;

    lengthline = 0.01;
    dlength = 0.005;

    _Sbase = 1;
    _Vbase = 0.4;
    _Zbase = _Vbase * _Vbase / _Sbase;
	_V0 = 1;
	_theta0 = 0;
	_AC = true;
    step = 1;

	N =  nAgent;
    B =  nBus;
    L =  nBus - 1;
    rapContr = nDeep / nBranch;
	proba = rapContr / (1 + rapContr);
    radial = true;
    MatConso = [[0] for i=(1:N)];
    Mat = [[]];
    	
	ref[:branch] = Dict{Int64, Any}();
    ref[:arcs]   = Vector{Tuple{Int64, Int64, Int64}}(undef, L);
   	
	branch = Vector{Int64}(undef, 0);
	distZero = Vector{Int64}(undef, B);
    for i=(1:B)
        distZero[i] = 0;
    end
	sizeVector = 1;
    push!(branch, 1);
	
	for l = 1:L 
		j = l + 1;
		i = 0; # random dans vecteur !!!
		nRandom = rand();
		dist = distZero[j];
		sizeVector = length(branch);
		if (sizeVector > nBranch || nRandom < proba) # on poursuit une branche existante
			indice = 0;
			while true
				if (sizeVector > nBranch) # si trop de branche on évite de partir du premier bus
					indice = floor(Int64, (rand() * (sizeVector - 2) + 1)) + 1;
				else 
					indice = floor(Int64, (rand() * (sizeVector - 1))) + 1;
                end
				i = branch[indice];
				dist = distZero[i] + 1;
			    if (dist<=nDeep)
                    break;
                end
            end
			distZero[j] = dist;
			
			if (indice>1) #pas relie au premier bus
				branch[indice] = j;
			else 
                push!(branch, j);
            end
		else 
            while true
                i = floor(Int64, rand() * l) + 1;
				dist = distZero[i] + 1;
                if (dist <= nDeep)
                    break;
                end
            end
			distZero[j] = dist;
            push!(branch, j);
        end
		
		Longueur = lengthline + 2 * (rand() - 0.5) * dlength;
	    ZsRe = Longueur * BTLINE[1] / _Zbase;
	    ZsIm = Longueur * BTLINE[2] / _Zbase;
	    limitLine = 2 * BTLINE[4] * _Vbase / _Sbase; # limite puissance active
	    limitI = limitLine * limitLine; # limite courrant
	    YlsRe = ZsRe / (ZsRe * ZsRe + ZsIm * ZsIm);
	    YlsIm = -ZsIm / (ZsRe * ZsRe + ZsIm * ZsIm);
	    Ylp = 314.15 * Longueur * BTLINE[3] / (2 * billion * _Zbase); # b/2

	
        ref[:branch][l] = Dict{String, Any}();
        ref[:branch][l]["f_bus"] = i;
        ref[:branch][l]["t_bus"] = j;
        ref[:branch][l]["rate_a"] = limitLine;
        ref[:branch][l]["br_r"] = ZsRe;
        ref[:branch][l]["br_x"] = ZsIm;
        ref[:branch][l]["br_g"] = YlsRe;
        ref[:branch][l]["br_b"] = YlsIm;
        ref[:branch][l]["br_c"] = Ylp;
        ref[:branch][l]["tap"] = 1;
        ref[:branch][l]["shift"] = 0;

        ref[:arcs][l] = (l, i, j);
       
    end
	
    ref[:bus] = Dict{Int64, Any}();
	for i = 1:B # bound on voltage angle rad
        ref[:bus][i] = Dict{String, Any}();
        ref[:bus][i]["angmax"] = 3;
        ref[:bus][i]["angmin"] = -3;
        ref[:bus][i]["angtart"] = 0;
        ref[:bus][i]["vstart"]  = 1;

		ref[:bus][i]["gs"] = 0;
        ref[:bus][i]["bs"] = 0;

        ref[:bus][i]["vmin"] = 0.9;
        ref[:bus][i]["vmax"] = 1.1;

    end

    ref[:ref_buses] = Dict{Int64, Any}();
    ref[:ref_buses][1] = ref[:bus][1];
    
    linePreB = Vector{Int64}(undef, B); 
    linePreLine = Vector{Int64}(undef, L);
    fromb =  Vector{Int64}(undef, L);
    to_b  =  Vector{Int64}(undef, L);

    linePreLine[1] = 0;
    linePreB[1] = 0;
    for (k, (l, i, j)) in enumerate(ref[:arcs])
        fromb[l] = i;
        to = j;
        to_b[l] = j;
        linePreB[to] = i;
    end
    for l=2:L
        linePreLine[l] = linePreB[fromb[l]];
    end

    
    ref, Nc, Ng, Np = Agent_from_random(ref, nAgent, B, linePreB);
   

    for key in keys(ref)
        ref[key]= sort(ref[key], by =  x->x[1])
    end

  


    arcdict = Dict(a => k for (k, a) in enumerate(ref[:arcs])) # (id_branche, bus1, bus2) => idCores et (id_branche, bus2, bus1) => idCores
    busdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:bus]))
    gendict = Dict(k => i for (i, (k, v)) in enumerate(ref[:gen]))
    branchdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:branch]))

  

	return (
        bus = [
            (
                i = busdict[k],
                gs = v["gs"],
                bs = v["bs"]
                #vmax   = v["vmax"],
                #vmin   = v["vmin"] 
                )
             for (k, v) in ref[:bus]
        ],
        arc = [
        (
            i = k,
            rate_a = ref[:branch][l]["rate_a"],
            bus_f = busdict[i],
            bus_t = busdict[j],
        ) for (k, (l, i, j)) in enumerate(ref[:arcs])
        ],
        agent = [
            (
                i = gendict[k], 
                oldi = v["oldi"],
                cost1 = v["cost"][1],
                cost2 = v["cost"][2],
                cost1q = v["costg"][1],
                cost2q = v["costg"][2],
                lb = v["lb"],
                ub = v["ub"],
                lbq = v["lbq"],
                ubq = v["ubq"],
                Mn = v["nVoisin"],
                gamma = v["gamma"],
                bus = busdict[v["gen_bus"]], # busdict[v["gen_bus"]]
                line = v["line_pre"]>0 ? branchdict[v["line_pre"]] : 0  #branchdict[v["line_pre"]]
            ) for (k, v) in ref[:gen]
        ],
        branch = [
            begin
                f_idx = arcdict[i, branch["f_bus"], branch["t_bus"]] # identification f_idx, t_idx sont les indices dans arcdict selon ordre des bus
                g = branch["br_g"]
                b = branch["br_b"]
                tr = branch["tap"] .* cos.(branch["shift"])
                ti = branch["tap"] .* sin.(branch["shift"])
                ttm = tr^2 + ti^2
                g_fr = 0 # je ne sais pas ce que c'est mais cela vaut 0 ? Peut être transfo
                b_fr = 0
                g_to = 0
                b_to = 0
                c1 = (-g * tr - b * ti) / ttm
                c2 = (-b * tr + g * ti) / ttm
                c3 = (-g * tr + b * ti) / ttm
                c4 = (-b * tr - g * ti) / ttm
                c5 = (g + g_fr) / ttm
                c6 = (b + b_fr) / ttm
                c7 = (g + g_to)
                c8 = (b + b_to)
                (
                    i = branchdict[i], #branchdict
                    f_idx = f_idx, 
                    bus_f = busdict[branch["f_bus"]], # busdict
                    bus_t = busdict[branch["t_bus"]], # busdict
                    x = branch["br_x"],
                    r = branch["br_r"],
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    c4 = c4,
                    c5 = c5,
                    c6 = c6,
                    c7 = c7,
                    c8 = c8,
                    rate_a_sq = branch["rate_a"]^2,
                    linePre = linePreLine[i] > 0 ? branchdict[linePreLine[i]] : 0 #branchdict[linePreLine[i]]
                )
            end for (i, branch) in ref[:branch]
        ],
        ref_buses = [busdict[i] for (i, k) in ref[:ref_buses]],
        vmax = [v["vmax"] for (k, v) in ref[:bus]],
        vmin = [v["vmin"] for (k, v) in ref[:bus]],
        pmax = [v["pmax"] for (k, v) in ref[:gen]],
        pmin = [v["pmin"] for (k, v) in ref[:gen]],
        qmax = [v["qmax"] for (k, v) in ref[:gen]],
        qmin = [v["qmin"] for (k, v) in ref[:gen]],
        pstart = [v["pstart"] for (k, v) in ref[:gen]],
        qstart = [v["qstart"] for (k, v) in ref[:gen]],
        rate_a = [ref[:branch][l]["rate_a"] for (k, (l, i, j)) in enumerate(ref[:arcs])],
        angmax = [b["angmax"] for (i, b) in ref[:bus]],
        angmin = [b["angmin"] for (i, b) in ref[:bus]],
        Nc = Nc, 
        Ng = Ng,
        Np = Np,
        ConsoAll = MatConso,
        MatAgent = Mat,
        step = 1,
    )
end

function Agent_from_random(ref, N, B, linePreB)
    _Sbase = 1;
    Pconso = 0.5;
    dP = 0.1; #10%
    dPconso = 0.1;
    dQconso = 0.005;
    Pprod = 2;
    dPprod = 0.1;
    bProd = 1;
    dbProd = 0.1;
    propCons = 0.50;
    propGenNfle = 0.25;

    gam = 2 + 2 * (rand() - 0.5);
    
    nCons = floor(Int64, propCons * (N-1)) + 1;
    nGenNFle = floor(Int64, propGenNfle * (N-1));
    nGenSup = N - nCons - nGenNFle;
    nGen = N - nCons;


    nPro = 0;
    
    a = Vector{Float64}(undef, 2*N);
    b = Vector{Float64}(undef, 2*N);
    
    pstart = Vector{Float64}(undef, 2*N);
    pmin = Vector{Float64}(undef, 2*N);
    pmax = Vector{Float64}(undef, 2*N);
    Lb = Vector{Float64}(undef, 2*N);
    Ub = Vector{Float64}(undef, 2*N);
    nVoisin = Vector{Float64}(undef, 2*N);
    gamma = Vector{Float64}(undef, N);
    
    lineAgent = Vector{Int64}(undef, N);
    oldIndice = Vector{Int64}(undef, N);
    bus = Vector{Int64}(undef, N);

     # loss agent
    oldIndice[1] = 1;
    bus[1] = 1;
    lineAgent[1] = linePreB[bus[1]];

	Lb[1] =  -10000; # pour ne pas avoir besoin de le modifier
    Ub[1] =  0;
	Ub[N + 1] =   10000; # idem
    Lb[N + 1] =  -10000;

    a[1] = 0;
    b[1] = 0;
    a[N + 1] = 0;
    b[N + 1] = 0;

    gamma[1] = 0;

    pstart[1] = 0;
    pmax[1] = 0;
    pmin[1] = 0;
    pstart[N + 1] = 0;
    pmax[N + 1] = 0;
    pmin[N + 1] = 0;
	nVoisin[1] = nGen + nPro;
	nVoisin[N + 1] = N - 1;
   

    for id=(2:N)
        if (id <= nCons) # consumer
            P0 = Pconso + dPconso * 2 * (rand() - 0.5);
            Q0 = dQconso * 2 * (rand() - 0.5);
            pLim1 = -(1 + dP) * P0 / _Sbase;
            pLim2 = -(1 - dP) * P0 / _Sbase;
            qLim1 = Q0 / _Sbase * (1 + dP - 2 * dP * (Q0 > 0));
            qLim2 = Q0 / _Sbase * (1 - dP + 2 * dP * (Q0 > 0));
            cost1 = 1 * (_Sbase * _Sbase);
            cost2 = P0 * cost1 /_Sbase;
            costQ1 = 0.1 * (_Sbase * _Sbase);
            costQ2 = -Q0 * costQ1 / _Sbase;
            Mn = nGen + nPro;
            pstart[id] = -P0 / _Sbase;
            pstart[id + N] = Q0;
            
            Ub[id] = 0;
            Lb[id] = pLim1;
            gamma[id] = - gam;

        elseif (id <= nCons + nGenNFle) 
            P0 = Pconso + dPconso * 2 * (rand() - 0.5);
            Q0 = dQconso * 2 * (rand() - 0.5);
            pLim1 = (1 - dP) * P0 / _Sbase;
            pLim2 = (1 + dP) * P0 / _Sbase;
            qLim1 = Q0 / _Sbase * (1 + dP - 2 * dP * (Q0 > 0));
            qLim2 = Q0 / _Sbase * (1 - dP + 2 * dP * (Q0 > 0));
            cost1 = 0.1 * (_Sbase * _Sbase);
            cost2 = - P0 * cost1 / _Sbase;
            costQ1 = 0.1 * (_Sbase * _Sbase);
            costQ2 = -Q0 * costQ1 / _Sbase;
            Mn = nCons + nPro;
            pstart[id] = P0 / _Sbase;
            pstart[id + N] = Q0;
            Ub[id] = pLim2;
            Lb[id] = 0;
            gamma[id] = gam;

        else  # generator
            P0 = Pconso + dPconso * 2 * (rand() - 0.5);
            P02 = Pprod + dPprod  * 2 * (rand() - 0.5);

            if (P0 < P02) 
                P0 = P02;
            end
            Q0 = 0;
            pLim1 = 0;
            pLim2 = P0 / _Sbase;

            qLim1 = - 2*dQconso;
            qLim2 =   2*dQconso;
            cost1 = 0.1 * (_Sbase * _Sbase);
            cost2 = bProd * _Sbase + dbProd * 2 * (rand() - 0.5) * _Sbase;
            costQ1 = 0.1 * (_Sbase * _Sbase);
            costQ2 = -Q0 * costQ1 / _Sbase;
            Mn = nCons + nPro;
            Ub[id] = pLim2;
            Lb[id] = 0;
            pstart[id] = P0 / _Sbase;
            pstart[id + N] = Q0;
            gamma[id] = gam;
        end
        bus[id] = floor(Int64, rand() * B) + 1;
        lineAgent[id] = linePreB[bus[id]];
        oldIndice[id] = id;
        nVoisin[ id ] = Mn;
	
		pmin[ id ] = pLim1;
		pmax[ id ] = pLim2;

		a[ id ] = cost1;
		b[ id ] = cost2;

		a[ id + N ] = costQ1;
		b[ id + N ] = costQ2;
		
		pmin[ id + N ] = qLim1;
		pmax[ id + N ] = qLim2;
		
		nVoisin[ id + N ] = N - 1;
		Ub[ id + N ] = qLim2 * (qLim2 > 0);
		Lb[ id + N ] = qLim1 * (qLim1 < 0);
    end
    ref[:gen] = Dict{Int64, Any}();
    for n in (1:N)
        ref[:gen][n] = Dict{String, Any}();
        ref[:gen][n]["cost"]     = [a[n] b[n] 0];
        ref[:gen][n]["gen_bus"]  = bus[n];
        ref[:gen][n]["line_pre"] = lineAgent[n]; 
        ref[:gen][n]["pmin"]     = pmin[n];
        ref[:gen][n]["pmax"]     = pmax[n];
        ref[:gen][n]["pstart"]   = pstart[n];
        ref[:gen][n]["costg"]    = [a[n + N] b[n + N] 0];
        ref[:gen][n]["qmin"]     = pmin[n + N];
        ref[:gen][n]["qmax"]     = pmax[n + N];
        ref[:gen][n]["qstart"]   = pstart[n + N];
        ref[:gen][n]["oldi"]     = oldIndice[n];

        ref[:gen][n]["nVoisin"]  = nVoisin[n];
        ref[:gen][n]["lb"]       = Lb[n];
        ref[:gen][n]["ub"]       = Ub[n];
        ref[:gen][n]["lbq"]      = Lb[n + N];
        ref[:gen][n]["ubq"]      = Ub[n + N];
        ref[:gen][n]["gamma"]    = gamma[n];
    end

    return ref, nCons, nGen, nPro

end


function updateData(data, step)  
    data.step = step
    for a in data.agent 
        i = a.i
        oldi = a.oldi
		factor = MatAgent[oldi - 1][2];
		PF =  MatAgent[oldi - 1][3]
        Pobj = ConsoAll[oldi - 1][step];

        S0 = factor * Pobj;
        P0 = - S0 * PF
       
       	a.cost2 = - P0;

		pLim1 = (1 + dP) * P0;
		pLim2 = (1 - dP) * P0;

        a.lb = pLim1;

        data.pstart[i] = P0;
        data.pmax[i] = pLim2;
        data.pmin[i] = pLim1;


        Q0 = S0 * sqrt(1 - PF * PF);
        signe = 2 * (oldi % 5) - 1; # un sur 5 ?
		Q0 = signe * Q0; 
        
        data.qstart[i] = Q0;

        qLim1 = Q0 * (1 - dP + 2 * dP * (Q0 < 0));
        qLim2 = Q0 * (1 + dP - 2 * dP * (Q0 < 0));

		a.costq2 = - 0.1 * Q0;
		
        a.lbq = qLim1;
		
		data.pmin[ i ] = qLim1;
		data.pmax[ i ] = qLim2;
		
		a.ubq = qLim2 * (qLim2 > 0);
		a.lbq = qLim1 * (qLim1 < 0);

	end

    return data
end

convert_data(data::N, backend) where {names,N<:NamedTuple{names}} =
    NamedTuple{names}(ExaModels.convert_array(d, backend) for d in data)

parse_ac_power_data(filename, backend) =
    convert_data(parse_ac_power_data(filename), backend)


filename = "ExaModelsPower.jl-main/data/case10babis.m" #case10babis.m testFeeder

#println(data)
function testMatpower(filename)
    PowerModels.silence()
    #data2 = PowerModels.parse_file(filename)
    #data = parse_ac_power_data(filename, nothing);
    #data = Data_from_AC("ACGrid/", "TestFeeder")
    data = Data_from_random(20, 20, 20, 20);

    Mustsolve = true;
    affichage = false;
    
    #println(keys(data));
   
    B = length(data.bus);
    L = length(data.branch);
    N = length(data.agent);
    println("N = $(N)");

    indiceLine = Vector{Int64}(undef, 0); # sauf les lignes sans antécedant
    FirstLine  = Vector{Int64}(undef, 0);
    
    indexLine = 1;
    for lin in data.branch
        #println(lin);
        if lin.linePre == 0
            push!(FirstLine, indexLine);
        else
            push!(indiceLine, indexLine);
        end
        indexLine+= 1;
        #println("-------")
    end

    indiceAgent   = Vector{Int64}(undef, 0); # sauf agent sur premier bus
    agentFirstBus = Vector{Int64}(undef, 0);
    indiceCompare = Vector{Int64}(undef, N); # passer de mon code a ce code H[i] = k

    indexAgent = 1;
    for agent in data.agent
        #println(agent);
        if agent.line > 0
            push!(indiceAgent, indexAgent);
        else
            push!(agentFirstBus, indexAgent);
        end
        indiceCompare[agent.i] = agent.oldi;
        indexAgent += 1;
    end
    
    vmin = Vector{Float64}(undef, B);
    vmax = Vector{Float64}(undef, B);

    for bus in data.bus
        #println(bus)
        vmin[bus.i] = data.vmin[bus.i] * data.vmin[bus.i];
        vmax[bus.i] = data.vmax[bus.i] * data.vmax[bus.i];
    end

    #println(indiceLine);
    #println(FirstLine);

    #println(indiceAgent)
    #println(data.agent[agentFirstBus])

    #println(data.pmin);
    #println(data.pmax);
    println(indiceCompare)
 
    ## DEBUT
    w = ExaCore()

    ## VARIABLE
    pg = variable(w, N; lvar=data.pmin, uvar = data.pmax, start = data.pstart); # lvar = pmin[1:N], uvar = pmax[1:N], start = pstart[1:N]); # 
    qg = variable(w, N; lvar=data.qmin, uvar = data.qmax, start = data.qstart); #lvar = pmin[N+1:end], uvar = pmax[N+1:end], start = pstart[N+1:end]);
    p  = variable(w, L); #; lvar = -data.rate_a, uvar = data.rate_a
    q  = variable(w, L); #; lvar = -data.rate_a, uvar = data.rate_a
    v  = variable(w, B; start = fill!(similar(data.bus, Float64), 1.0) , lvar = vmin, uvar = vmax)
    l  = variable(w, L; lvar = 0);
    
    ## OBJECTIF
    objective(w, (0.5*g.cost1) * pg[g.i]^2 + g.cost2 * pg[g.i] + (0.5 * g.cost1q) * qg[g.i]^2 + g.cost2q * qg[g.i] for g in data.agent);
    
    ## contraintes
    constraint(w, l[lin.i] * v[lin.bus_t] - (p[lin.i] * p[lin.i] + q[lin.i] * q[lin.i]) for lin in data.branch; lcon = 0, ucon = Inf)
    constraint(w, v[lin.bus_f] - v[lin.bus_t] + 2 * (lin.r *p[lin.i] + lin.x * q[lin.i]) - (lin.r*lin.r + lin.x*lin.x)*l[lin.i] for lin in data.branch);
    
    # constraint on all line
    c1 = constraint(w, -p[lin.i] for lin in data.branch);
    c3 = constraint(w, -q[lin.i] for lin in data.branch);

    # add flow of the child (except for line from the ref bus)
    constraint!(w, c1, lin.linePre => p[lin.i] - lin.r * l[lin.i] for lin in data.branch[indiceLine]);
    constraint!(w, c3, lin.linePre => q[lin.i] - lin.x * l[lin.i] for lin in data.branch[indiceLine]);

    # add agent power on the bus except the agent on first bus...
    constraint!(w, c1, ag.line => pg[ag.i] for ag in data.agent[indiceAgent])
    constraint!(w, c3, ag.line => qg[ag.i] for ag in data.agent[indiceAgent])

    #constraint on the first bus
    c2 = constraint(w, 1); 
    c4 = constraint(w, 1); 

    constraint!(w, c2, 1 => p[lin.i] - lin.r * l[lin.i] for lin in data.branch[FirstLine]);
    constraint!(w, c4, 1 => q[lin.i] - lin.x * l[lin.i] for lin in data.branch[FirstLine]);

    constraint!(w, c2, 1 => pg[ag.i] for ag in data.agent[agentFirstBus]);
    constraint!(w, c4, 1 => qg[ag.i] for ag in data.agent[agentFirstBus]);


    #==#
    ## RESULTATS
    if(Mustsolve)
        m = ExaModel(w)
        result = ipopt(m)

        # 5 4 6 7 2 10 9 8 3 1
        SolutionP11 = [0 -0.178851 -0.0924302 -0.172834 -0.153096 -0.153712 -0.0703927 -0.107012 -0.0895029 -0.155118 1.24157 0];
        SolutionQ11 = [0 -0.0476596 -0.0354596 -0.0458067 -0.184968 -0.0608347 -0.01155 -0.0063 -0.0135995 -0.0205341 0.416734 0.101781];

        SolutionP50 = [0 -0.000208 -0.003232 -0.0063537 -0.0024 -0.00224 -0.0134927 -0.013641 -0.00064 -0.00064 -0.00406252 -0.00548409 (
        -0.00540993) -8e-05 -0.0108358 -0.000424 -0.00225744 -0.00112 -0.00112 -0.00208 -0.00208 -0.00112 -0.00156 -0.00048 -0.00208 (
        -0.00208) -0.00192 -0.00192 -9.6e-05 -0.00048 -0.00314176 -0.00314924 -0.00673133 -0.0373626 -0.0373917 -0.00324 -0.000288 (
        -0.000344) -0.002112 -0.00192 -0.0087873 -0.123145 -0.00256 -0.0213971 -0.00472 -0.00144 -0.00144 -0.00224 -0.00224 0.381397 0];
        
        SolutionQ50 = [0 -0.000209 -0.00298305 -0.00535306 -0.00209554 -0.00187425 -0.0103475 -0.0103768 -0.000575033 -0.0005775 -0.00307496 (
        -0.0035635) -0.00353744 -6.3e-05 -0.00813927 -0.0003675 -0.00203093 -0.0010393 -0.00103787 -0.001767 -0.001767 -0.000957453 (
        -0.00137856) -0.000385211 -0.001767 -0.001767 -0.0016224 -0.001615 -0.000103562 -0.0004515 -0.00265234 -0.00263493 -0.00552149 (
        -0.0273718) -0.0273836 -0.00278783 -0.0002565 -0.0003325 -0.001805 -0.001634 -0.00704829 -0.0886158 -0.002185 -0.0159936 -0.00399456 (
        -0.00125424) -0.00125675 -0.00198275 -0.00198114 0.269557 0.0103033]

        SolutionP60 = [0 -0.0612944 -0.0297706 -0.0287561 -0.0487483 -0.028224 -0.028224 -0.042336 -0.133635 -0.0512898 -0.0515523 -0.0309709 (
        -0.0310553) -0.0310925 -0.0516153 -0.0293314 -0.0286379 -0.0492699 -0.0492005 -0.0284963 -0.028521 -0.0112 -0.0286596 -0.0492496 (
        -0.0492225) -0.0491902 -0.0285476 -0.0285399 -0.028539 -0.0286324 -0.0286365 -0.0286364 -0.0112 -0.0296962 -0.049418 -0.028672 (
        -0.0493918) -0.0493868 -0.0112 -0.0486353 -0.0482984 -0.0479474 -0.0479052 -0.0112 -0.0475018 -0.0471892 -0.028224 -0.0473702 (
        -0.0471736) -0.028224 -0.0471021 -0.0112 -0.0488216 -0.028224 -0.0485766 -0.0484914 -0.028224 -0.0112 -0.028224 2.4751 0];
        
        SolutionQ60 = [0 -0.059988 -0.0341932 -0.0341933 -0.0550303 -0.0341932 -0.0341932 -0.0377924 -0.119976 -0.0547304 -0.0546042 (
        -0.0341932) -0.0341932 -0.0341932 -0.0545414 -0.0341932 -0.0341932 -0.0548729 -0.0548634 -0.0341932 -0.0341932 -0.0135688 (
        -0.0341939) -0.0548217 -0.0547996 -0.0548942 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0135688 (
        -0.0353132) -0.0554443 -0.0342743 -0.0554207 -0.0554011 -0.0135688 -0.0555326 -0.0558424 -0.0561716 -0.0561681 -0.0136579 (
        -0.0569065) -0.0574853 -0.0364788 -0.0571589 -0.0574538 -0.0363115 -0.0576154 -0.0140572 -0.0550796 -0.0360171 -0.0549047 (
        -0.0548569) -0.0341932 -0.0135688 -0.0341932 2.5207 0.143395];


        SolutionP56 = [0 -0.03078 -0.03762 -0.04617 -0.04104 -0.029925 -0.04446 -0.047025 -0.036765 -0.050445 -0.048735 -0.041895 -0.07866 -0.04788 -0.041895 -0.045315 -0.074385 -0.04104 -0.041895 -0.045315 -0.043605 -0.038475 -0.04617 -0.047025 -0.04617 -0.04275 -0.02907 -0.043605 -0.02394 -0.04617 -0.045315 -0.045315 -0.04788 -0.048735 -0.038475 -0.036765 -0.047025 -0.038475 -0.04446 -0.04104 -0.036765 -0.04446 -0.04617 -0.04275 -0.024795 -0.045315 -0.05301 -0.038475 -0.03933 -0.04617 -0.04446 -0.045315 -0.04446 -0.040185 -0.04104 -0.047025 2.39143];
        SolutionQ56 = [0 0.0123651 0.0431979 0.0853917 0.105345 -0.0120216 0.0178607 0.0531953 0.0686578 0.128271 -0.0195781 0.0168303 0.0872372 0.0884508 0.107451 -0.0182042 0.0289232 0.0467987 0.0777981 0.115537 -0.0175172 0.0154564 0.052288 0.0869319 0.11735 -0.0171737 0.0116781 0.0495438 0.0455985 0.11735 -0.0182042 0.0182042 0.054106 0.0899634 0.0989976 -0.0147694 0.0188911 0.0440996 0.0823245 0.105345 -0.0147694 0.0178607 0.052288 0.079289 0.0649743 -0.0182042 0.0209879 0.0440996 0.0732857 0.11735 -0.0178607 0.0182042 0.0504611 0.0748093 0.105345 -0.0188911 -2.79631];

                #println(result);
        println("Status: $(result.status)")
        println("Number of iterations: $(result.iter)")

       
        solP = solution(result, pg)
        
        println("result                        : ", solP);
        
        if N==11 || N == 12
            println("vrai reponse (cas10ba) : ", SolutionP11[indiceCompare]);
        elseif N==56 || N == 57
            println("vrai reponse (cas TestFeeder) : ", SolutionP56[indiceCompare]);
        elseif  N==50 || N == 51
            println("vrai reponse (cas69) : ", SolutionP50[indiceCompare])
        elseif  N==60 || N == 61
            println("vrai reponse (case85) : ", SolutionP60[indiceCompare])
        end

        solp = solution(result, p)
        #println(solp)
        solQ = solution(result, qg)
        #println(solQ)
        solq = solution(result, q)
        #println(solq)
        soll = solution(result, l)
        #println(soll)
        solv = solution(result, v)
        
        if affichage
            for i=(1:B)
                print(sqrt(solv[i]), " ")
            end
            println(" ")
        end

        sum = 0;
        for g in data.agent
            sum += (0.5 * g.cost1) * solP[g.i]^2 + g.cost2 * solP[g.i] + (0.5 * g.cost1q) * solQ[g.i]^2 + g.cost2q * solQ[g.i] 
        end

        sum2 = 0;

        if N == 11 || N == 12
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP11[g.oldi]^2 + g.cost2 * SolutionP11[g.oldi] + (0.5 * g.cost1q) * SolutionQ11[g.oldi]^2 + g.cost2q * SolutionQ11[g.oldi] 
            end
        elseif N == 56 || N == 57
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP56[g.oldi]^2 + g.cost2 * SolutionP56[g.oldi] + (0.5 * g.cost1q) * SolutionQ56[g.oldi]^2 + g.cost2q * SolutionQ56[g.oldi] 
            end
        elseif N == 50 || N == 51
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP50[g.oldi]^2 + g.cost2 * SolutionP50[g.oldi] + (0.5 * g.cost1q) * SolutionQ50[g.oldi]^2 + g.cost2q * SolutionQ50[g.oldi] 
            end
        elseif N == 60 || N == 61
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP60[g.oldi]^2 + g.cost2 * SolutionP60[g.oldi] + (0.5 * g.cost1q) * SolutionQ60[g.oldi]^2 + g.cost2q * SolutionQ60[g.oldi] 
            end
        end
        println("Cost function : $(sum)")
        println("Cost function : $(sum2)")

        

        const1 = zeros(L);
        #c1 = constraint(w, -p[lin.i] for lin in data.branch); 
        # constraint!(w, c1, lin.linePre => p[lin.i] - lin.r * l[lin.i] for lin in data.branch[indiceLine]);
        # constraint!(w, c1, ag.line => pg[ag.i] for ag in data.agent[indiceAgent])

        index = 1;
        for lin in data.branch
            const1[index] += -solp[lin.i];
            index += 1;
        end
        for lin in data.branch[indiceLine]
            index = lin.linePre;
            const1[index] += (solp[lin.i] - lin.r * soll[lin.i]);
        end
        for ag in data.agent[indiceAgent]
            index = ag.line
            const1[index] += solP[ag.i];
        end

        if affichage
            println("contraite : ")
            for lin in data.branch
                print(lin.i, " ", soll[lin.i] * solv[lin.bus_t] - (solp[lin.i] * solp[lin.i] + solq[lin.i] * solq[lin.i]) )
                print(" ", const1[lin.i]);
                println(" ", solv[lin.bus_f] - solv[lin.bus_t] + 2 * (lin.r *solp[lin.i] + lin.x * solq[lin.i]) - (lin.r*lin.r + lin.x*lin.x)*soll[lin.i])
            end
        end
    end
end


function getIndice(i, j, Nc, Ng, Np)
    OffsetG = Nc * (Ng + Np);
    OffsetP = Nc * (Ng + Np) + Ng * (Nc + Np);
    Mc = Ng + Np;
    Mg = Nc + Np;
    Mp = Nc + Ng;


    if i <= Nc # i est un conso
        if j <= Ng # j est un gen
            return OffsetG + Mg * (j - 1) + i, j + Nc
        else
            return OffsetP + Mp * (j - Ng - 1) + i, j + Nc
        end
    elseif i <= Nc + Ng # i est un gen
        if j <= Nc # j est un conso
            return Mc * (j - 1) + (i - Nc), j
        else
            return OffsetP + Mp * (j - Nc - 1) + (i - Nc), j + Ng
        end
    else # i est un prod
        if j <= Nc # j est un conso
            return Mc * (j - 1) + (i -(Nc + Ng)), j
        else
            return OffsetG + Mp * (j - Nc - 1) + (i - (Nc + Ng)), j
        end
    end

end


function testMarket(filename)
    PowerModels.silence()
    #data2 = PowerModels.parse_file(filename)
    #data = parse_ac_power_data(filename, nothing);
    data = Data_from_AC("ACGrid/", "case85")
    #data = Data_from_random(200, 2, 2, 2);

    Mustsolve = true;
    
    N = length(data.agent);
    Nc = data.Nc;
    Ng = data.Ng;
    Np = data.Np;

    M = Nc * (Ng + Np) + Ng * (Nc + Np) + Np * (Nc + Ng);
    Mq = (N - 1) * N
    
    println("N : $(N) Nc : $(Nc), Ng : $(Ng), Np : $(Np), M : $(M), Mq : $(Mq)")
    for a in data.agent
        #print(a.Mn, " ")
    end
    println(" ")

    Ubp = Vector{Float64}(undef, M);
    Lbp = Vector{Float64}(undef, M);
    Ubq = Vector{Float64}(undef, Mq);
    Lbq = Vector{Float64}(undef, Mq);
    tradesp = Vector{@NamedTuple{i::Int64, n::Int64, v::Int64, m::Int64, g::Float64}}(undef, M);
    tradesq = Vector{@NamedTuple{i::Int64, n::Int64, v::Int64, m::Int64}}(undef, Mq);

    indicep = 1;
    indiceq = 1;
    indiceCompare = Vector{Int64}(undef, N); # passer de mon code a ce code H[i] = k

    for a in data.agent
        n = a.i;
        indiceCompare[a.i] = a.oldi;
        for j=(1:a.Mn)
            Ubp[indicep] =  a.ub
            Lbp[indicep] =  a.lb
            v, m = getIndice(n,j, Nc, Ng, Np);
            tradesp[indicep] = (i=indicep, n = n, v = v, m = m, g = a.gamma);
            indicep += 1;
        end
        for j=(1:N-1)
            Ubq[indiceq] =  a.ubq
            Lbq[indiceq] =  a.lbq
            m = 1;
            v = 1;
            if (n <= j)
                m = j + 1;
                v = (m - 1) * (N - 1) + n;
            else
                m = j;
                v = (m - 1) * (N - 1) + (n - 1);
            end
            tradesq[indiceq] = (i=indiceq, n = n, v = v, m = m);
            indiceq += 1;
        end
    end

    #for n=(1:N)
    #   println("n : $(n) -> pmin : $(data.pmin[n]), pobj : $(data.pstart[n]), pmax : $(data.pmax[n]), qmin : $(data.qmin[n]), qobj : $(data.qstart[n]), qmax : $(data.qmax[n])")
    #end


    data.pmin[1] = -1e-12; 
    data.qmin[1] = -1e-12;
    data.qmax[1] =  1e-12;


    #println(sort(tradesp, by=(x ->x.n)));  

    #w = ExaCore(Float64, backend = CUDABackend())
    w = ExaCore(Float64, backend = nothing)
   
    Tp = variable(w, M;  lvar=Lbp, uvar = Ubp); # M
    Tq = variable(w, Mq; lvar=Lbq, uvar = Ubq); # Mq
    pg = variable(w, N; lvar=data.pmin, uvar = data.pmax, start = data.pstart); # N 
    qg = variable(w, N; lvar=data.qmin, uvar = data.qmax, start = data.qstart); # N 

    objective(w, (0.5*g.cost1) * pg[g.i]^2 + g.cost2 * pg[g.i] + (0.5 * g.cost1q) * qg[g.i]^2 + g.cost2q * qg[g.i] for g in data.agent);
    objective(w, t.g * Tp[t.i] for t in tradesp);

    #p = sum(T)
    c1 = constraint(w, -pg[ag.i] for ag in data.agent); # N
    constraint!(w, c1, t.n => Tp[t.i] for t in tradesp);

    c2 = constraint(w, -qg[ag.i] for ag in data.agent); # N
    constraint!(w, c2, t.n => Tq[t.i] for t in tradesq);


    # tij = -tij
    constraint(w, Tp[t.i] + Tp[t.v] for t in tradesp); # 40 ?
    constraint(w, Tq[t.i] + Tq[t.v] for t in tradesq); # 132 ?

    if (Mustsolve)
        m = ExaModel(w)
        result = madnlp(m, tol = 1e-4); #ipopt(m) madnlp

        # 5 4 6 7 2 10 9 8 3 1
        SolutionP11 = [0 -0.184 -0.098 -0.179 -0.1598 -0.161 -0.078 -0.115 -0.098 -0.164 1.2368 0];
        SolutionQ11 = [0 -0.0483 -0.0357 -0.04683 -0.1932 -0.063 -0.01155 -0.0063 -0.01365 -0.021 0.39767 0.0451542];

        SolutionP50 = [0 -0.000311999 -0.00411889 -0.00757786 -0.00308694 -0.00287958 -0.0145634 -0.0145634 -0.000805472 -0.000805472 -0.00462926 -0.00607959 -0.00607959 -0.00012 -0.0114904 -0.000461994 -0.00287958 -0.00149912 -0.00149912 -0.00269867 -0.00269867 -0.00149912 -0.00204442 -0.000533764 -0.00269867 -0.00269867 -0.0024912 -0.0024912 -0.000144 -0.000533764 -0.00399849 -0.00399849 -0.00797228 -0.0385337 -0.0385337 -0.00412892 -0.000431998 -0.000515997 -0.00271752 -0.0024912 -0.0100905 -0.124469 -0.00328858 -0.0227645 -0.00597969 -0.00189253 -0.00189253 -0.00287958 -0.00287958 0.380215 0];
         
        
        SolutionQ50 = [0 -0.000231 -0.00315 -0.00567 -0.00231 -0.001995 -0.0108378 -0.0108378 -0.0005775 -0.0005775 -0.00315 -0.003675 -0.003675 -6.3e-05 -0.008505 -0.0003675 -0.0021 -0.00105 -0.00105 -0.001953 -0.001953 -0.00105 -0.00147 -0.00042 -0.001953 -0.001953 -0.001785 -0.001785 -0.000105 -0.0004515 -0.0027615 -0.0027615 -0.005922 -0.0280524 -0.0280524 -0.0029715 -0.0002835 -0.0003675 -0.001995 -0.001806 -0.00756 -0.0883956 -0.002415 -0.0167882 -0.00441 -0.001365 -0.001365 -0.0021 -0.0021 0.268828 0.00977632];

        SolutionP60 = [0 -0.0612944 -0.0297706 -0.0287561 -0.0487483 -0.028224 -0.028224 -0.042336 -0.133635 -0.0512898 -0.0515523 -0.0309709 (
        -0.0310553) -0.0310925 -0.0516153 -0.0293314 -0.0286379 -0.0492699 -0.0492005 -0.0284963 -0.028521 -0.0112 -0.0286596 -0.0492496 (
        -0.0492225) -0.0491902 -0.0285476 -0.0285399 -0.028539 -0.0286324 -0.0286365 -0.0286364 -0.0112 -0.0296962 -0.049418 -0.028672 (
        -0.0493918) -0.0493868 -0.0112 -0.0486353 -0.0482984 -0.0479474 -0.0479052 -0.0112 -0.0475018 -0.0471892 -0.028224 -0.0473702 (
        -0.0471736) -0.028224 -0.0471021 -0.0112 -0.0488216 -0.028224 -0.0485766 -0.0484914 -0.028224 -0.0112 -0.028224 2.4751 0];
        
        SolutionQ60 = [0 -0.059988 -0.0341932 -0.0341933 -0.0550303 -0.0341932 -0.0341932 -0.0377924 -0.119976 -0.0547304 -0.0546042 (
        -0.0341932) -0.0341932 -0.0341932 -0.0545414 -0.0341932 -0.0341932 -0.0548729 -0.0548634 -0.0341932 -0.0341932 -0.0135688 (
        -0.0341939) -0.0548217 -0.0547996 -0.0548942 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0135688 (
        -0.0353132) -0.0554443 -0.0342743 -0.0554207 -0.0554011 -0.0135688 -0.0555326 -0.0558424 -0.0561716 -0.0561681 -0.0136579 (
        -0.0569065) -0.0574853 -0.0364788 -0.0571589 -0.0574538 -0.0363115 -0.0576154 -0.0140572 -0.0550796 -0.0360171 -0.0549047 (
        -0.0548569) -0.0341932 -0.0135688 -0.0341932 2.5207 0.143395];


        SolutionP56 = [0 -0.03078 -0.03762 -0.04617 -0.04104 -0.029925 -0.04446 -0.047025 -0.036765 -0.050445 -0.048735 -0.041895 -0.07866 -0.04788 -0.041895 -0.045315 -0.074385 -0.04104 -0.041895 -0.045315 -0.043605 -0.038475 -0.04617 -0.047025 -0.04617 -0.04275 -0.02907 -0.043605 -0.02394 -0.04617 -0.045315 -0.045315 -0.04788 -0.048735 -0.038475 -0.036765 -0.047025 -0.038475 -0.04446 -0.04104 -0.036765 -0.04446 -0.04617 -0.04275 -0.024795 -0.045315 -0.05301 -0.038475 -0.03933 -0.04617 -0.04446 -0.045315 -0.04446 -0.040185 -0.04104 -0.047025 2.39143];
        
        SolutionQ56 = [0 0.0123651 0.0431979 0.0853917 0.105345 -0.0120216 0.0178607 0.0531953 0.0686578 0.128271 -0.0195781 0.0168303 0.0872372 0.0884508 0.107451 -0.0182042 0.0289232 0.0467987 0.0777981 0.115537 -0.0175172 0.0154564 0.052288 0.0869319 0.11735 -0.0171737 0.0116781 0.0495438 0.0455985 0.11735 -0.0182042 0.0182042 0.054106 0.0899634 0.0989976 -0.0147694 0.0188911 0.0440996 0.0823245 0.105345 -0.0147694 0.0178607 0.052288 0.079289 0.0649743 -0.0182042 0.0209879 0.0440996 0.0732857 0.11735 -0.0178607 0.0182042 0.0504611 0.0748093 0.105345 -0.0188911 -2.79631];

        #println(result);
        println("Status: $(result.status)")
        println("Number of iterations: $(result.iter)")

       
        solPGPU = solution(result, pg)
        solP = Vector{Float64}(undef, N);
        copyto!(solP, solPGPU);
        println("result                        : ", solP);
        
        if N==11 || N == 12
            println("vrai reponse (cas10ba) : ", SolutionP11[indiceCompare]);
        elseif N==56 || N == 57
            println("vrai reponse (cas TestFeeder) : ", SolutionP56[indiceCompare]);
        elseif  N==50 || N == 51
            println("vrai reponse (cas69) : ", SolutionP50[indiceCompare])
        elseif  N==60 || N == 61
            println("vrai reponse (case85) : ", SolutionP60[indiceCompare])
        end

        solQGPU = solution(result, qg)
        solQ = Vector{Float64}(undef, N);
        copyto!(solQ, solQGPU);
        #println(solQ)
        

        sum = 0;
        for g in data.agent
            sum += (0.5 * g.cost1) * solP[g.i]^2 + g.cost2 * solP[g.i] + (0.5 * g.cost1q) * solQ[g.i]^2 + g.cost2q * solQ[g.i] 
        end

        sum2 = 0;
        if N == 11 || N == 12
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP11[g.oldi]^2 + g.cost2 * SolutionP11[g.oldi] + (0.5 * g.cost1q) * SolutionQ11[g.oldi]^2 + g.cost2q * SolutionQ11[g.oldi] 
            end
        elseif N == 56 || N == 57
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP56[g.oldi]^2 + g.cost2 * SolutionP56[g.oldi] + (0.5 * g.cost1q) * SolutionQ56[g.oldi]^2 + g.cost2q * SolutionQ56[g.oldi] 
            end
        elseif N == 50 || N == 51
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP50[g.oldi]^2 + g.cost2 * SolutionP50[g.oldi] + (0.5 * g.cost1q) * SolutionQ50[g.oldi]^2 + g.cost2q * SolutionQ50[g.oldi] 
            end
        elseif N == 60 || N == 61
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP60[g.oldi]^2 + g.cost2 * SolutionP60[g.oldi] + (0.5 * g.cost1q) * SolutionQ60[g.oldi]^2 + g.cost2q * SolutionQ60[g.oldi] 
            end
        end
        println("Cost function (julia): $(sum)")
        println("Cost function (C)    : $(sum2)")
    end

end

function testEndoMarketDirect(name, sizeAgent)
    PowerModels.silence()
    
    data = Data_from_AC("ACGrid/", "case85", sizeAgent)
    #data = Data_from_random(100, 30, 20, 20);

    Mustsolve = true;
    affichage = false;
    backend = CUDABackend() #CUDABackend();
    
    #println(keys(data));
   
    B = length(data.bus);
    L = length(data.branch);
    N = length(data.agent);

    Nc = data.Nc;
    Ng = data.Ng;
    Np = data.Np;

    M = Nc * (Ng + Np) + Ng * (Nc + Np) + Np * (Nc + Ng);
    Mq = (N - 1) * N

    Ubp = Vector{Float64}(undef, M);
    Lbp = Vector{Float64}(undef, M);
    Ubq = Vector{Float64}(undef, Mq);
    Lbq = Vector{Float64}(undef, Mq);
    tradesp = Vector{@NamedTuple{i::Int64, n::Int64, v::Int64, m::Int64, g::Float64}}(undef, M);
    tradesq = Vector{@NamedTuple{i::Int64, n::Int64, v::Int64, m::Int64}}(undef, Mq);

    indicep = 1;
    indiceq = 1;
    indiceCompare = Vector{Int64}(undef, N); # passer de mon code a ce code H[i] = k

    for a in data.agent
        n = a.i;
        indiceCompare[a.i] = a.oldi;
        for j=(1:a.Mn)
            Ubp[indicep] =  a.ub
            Lbp[indicep] =  a.lb
            v, m = getIndice(n,j, Nc, Ng, Np);
            tradesp[indicep] = (i=indicep, n = n, v = v, m = m, g = a.gamma);
            indicep += 1;
        end
        for j=(1:N-1)
            Ubq[indiceq] =  a.ubq
            Lbq[indiceq] =  a.lbq
            m = 1;
            v = 1;
            if (n <= j)
                m = j + 1;
                v = (m - 1) * (N - 1) + n;
            else
                m = j;
                v = (m - 1) * (N - 1) + (n - 1);
            end
            
            tradesq[indiceq] = (i=indiceq, n = n, v = v, m = m);
            indiceq += 1;
        end
    end

    data.pmin[1] = -1000; # pas de limite pour l'agent des pertes
    data.qmin[1] = -1000; # pas de limite pour l'agent des pertes
    data.qmax[1] = 1000; # pas de limite pour l'agent des pertes
    
    indiceLine = Vector{Int64}(undef, 0); # sauf les lignes sans antécedant
    FirstLine  = Vector{Int64}(undef, 0);
    
    indexLine = 1;
    for lin in data.branch
        #println(lin);
        if lin.linePre == 0
            push!(FirstLine, indexLine);
        else
            push!(indiceLine, indexLine);
        end
        indexLine+= 1;
        #println("-------")
    end

    indiceAgent   = Vector{Int64}(undef, 0); # sauf agent sur premier bus
    agentFirstBus = Vector{Int64}(undef, 0);
    indiceCompare = Vector{Int64}(undef, N); # passer de mon code a ce code H[i] = k

    indexAgent = 1;
    for agent in data.agent
        #println(agent);
        if agent.line > 0
            push!(indiceAgent, indexAgent);
        elseif (agent.oldi != 1)
            push!(agentFirstBus, indexAgent);
        end
        indiceCompare[agent.i] = agent.oldi;
        indexAgent += 1;
    end
    
    vmin = Vector{Float64}(undef, B);
    vmax = Vector{Float64}(undef, B);

    for bus in data.bus
        #println(bus)
        vmin[bus.i] = data.vmin[bus.i] * data.vmin[bus.i];
        vmax[bus.i] = data.vmax[bus.i] * data.vmax[bus.i];
    end

    ## DEBUT
    w = ExaCore(Float64, backend = backend)

    ## VARIABLE
    Tp = variable(w, M;  lvar=Lbp, uvar = Ubp); # M
    Tq = variable(w, Mq; lvar=Lbq, uvar = Ubq); # Mq
    pg = variable(w, N; lvar=data.pmin, uvar = data.pmax, start = data.pstart); # lvar = pmin[1:N], uvar = pmax[1:N], start = pstart[1:N]); # 
    qg = variable(w, N; lvar=data.qmin, uvar = data.qmax, start = data.qstart); #lvar = pmin[N+1:end], uvar = pmax[N+1:end], start = pstart[N+1:end]);
    p  = variable(w, L); #; lvar = -data.rate_a, uvar = data.rate_a
    q  = variable(w, L); #; lvar = -data.rate_a, uvar = data.rate_a
    v  = variable(w, B; start = fill!(similar(data.bus, Float64), 1.0) , lvar = vmin, uvar = vmax)
    l  = variable(w, L; lvar = 0);
    
    ## OBJECTIF
    objective(w, (0.5*g.cost1) * pg[g.i]^2 + g.cost2 * pg[g.i] + (0.5 * g.cost1q) * qg[g.i]^2 + g.cost2q * qg[g.i] for g in data.agent);
    objective(w, t.g * Tp[t.i] for t in tradesp);

    ## contraintes
    constraint(w, l[lin.i] * v[lin.bus_t] - (p[lin.i] * p[lin.i] + q[lin.i] * q[lin.i]) for lin in data.branch; lcon = 0, ucon = Inf)
    constraint(w, v[lin.bus_f] - v[lin.bus_t] + 2 * (lin.r *p[lin.i] + lin.x * q[lin.i]) - (lin.r*lin.r + lin.x*lin.x)*l[lin.i] for lin in data.branch);
    
    # constraint on all line
    c1 = constraint(w, -p[lin.i] for lin in data.branch);
    c3 = constraint(w, -q[lin.i] for lin in data.branch);

    # add flow of the child (except for line from the ref bus)
    constraint!(w, c1, lin.linePre => p[lin.i] - lin.r * l[lin.i] for lin in data.branch[indiceLine]);
    constraint!(w, c3, lin.linePre => q[lin.i] - lin.x * l[lin.i] for lin in data.branch[indiceLine]);

    # add agent power on the bus except the agent on first bus...
    constraint!(w, c1, ag.line => pg[ag.i] for ag in data.agent[indiceAgent])
    constraint!(w, c3, ag.line => qg[ag.i] for ag in data.agent[indiceAgent])

    #constraint on the first bus
    c2 = constraint(w, 1); 
    c4 = constraint(w, 1); 

    constraint!(w, c2, 1 => p[lin.i] - lin.r * l[lin.i] for lin in data.branch[FirstLine]);
    constraint!(w, c4, 1 => q[lin.i] - lin.x * l[lin.i] for lin in data.branch[FirstLine]);

    constraint!(w, c2, 1 => pg[ag.i] for ag in data.agent[agentFirstBus]);
    constraint!(w, c4, 1 => qg[ag.i] for ag in data.agent[agentFirstBus]);


     #p = sum(T)
    c5 = constraint(w, -pg[ag.i] for ag in data.agent); # N
    constraint!(w, c5, t.n => Tp[t.i] for t in tradesp);

    c6 = constraint(w, -qg[ag.i] for ag in data.agent); # N
    constraint!(w, c6, t.n => Tq[t.i] for t in tradesq);


    # tij = -tij
    constraint(w, Tp[t.i] + Tp[t.v] for t in tradesp); # 40 ?
    constraint(w, Tq[t.i] + Tq[t.v] for t in tradesq); # 132 ?

    #==#
    ## RESULTATS
    if(Mustsolve)
        m = ExaModel(w)
        t = @elapsed(result = madnlp(m, tol = 1e-4)) #ipopt(m) madnlp

        # 5 4 6 7 2 10 9 8 3 1
        SolutionP11 = [ -0.0698454 -0.177316 -0.0913216 -0.172273 -0.15302 -0.154027 -0.0709483 -0.107816 -0.0905968 -0.156435 1.2435 0];
        SolutionQ11 = [ -0.0933341 -0.0471442 -0.0351415 -0.0457168 -0.18509 -0.0610423 -0.01155 -0.0063 -0.01365 -0.0209497 0.417456 0.102503];
        
        SolutionP50 = [0 -0.000208 -0.003232 -0.0063537 -0.0024 -0.00224 -0.0134927 -0.013641 -0.00064 -0.00064 -0.00406252 -0.00548409 (
        -0.00540993) -8e-05 -0.0108358 -0.000424 -0.00225744 -0.00112 -0.00112 -0.00208 -0.00208 -0.00112 -0.00156 -0.00048 -0.00208 (
        -0.00208) -0.00192 -0.00192 -9.6e-05 -0.00048 -0.00314176 -0.00314924 -0.00673133 -0.0373626 -0.0373917 -0.00324 -0.000288 (
        -0.000344) -0.002112 -0.00192 -0.0087873 -0.123145 -0.00256 -0.0213971 -0.00472 -0.00144 -0.00144 -0.00224 -0.00224 0.381397 0];
        
        SolutionQ50 = [0 -0.000209 -0.00298305 -0.00535306 -0.00209554 -0.00187425 -0.0103475 -0.0103768 -0.000575033 -0.0005775 -0.00307496 (
        -0.0035635) -0.00353744 -6.3e-05 -0.00813927 -0.0003675 -0.00203093 -0.0010393 -0.00103787 -0.001767 -0.001767 -0.000957453 (
        -0.00137856) -0.000385211 -0.001767 -0.001767 -0.0016224 -0.001615 -0.000103562 -0.0004515 -0.00265234 -0.00263493 -0.00552149 (
        -0.0273718) -0.0273836 -0.00278783 -0.0002565 -0.0003325 -0.001805 -0.001634 -0.00704829 -0.0886158 -0.002185 -0.0159936 -0.00399456 (
        -0.00125424) -0.00125675 -0.00198275 -0.00198114 0.269557 0.0103033]

        SolutionP60 = [0 -0.0612944 -0.0297706 -0.0287561 -0.0487483 -0.028224 -0.028224 -0.042336 -0.133635 -0.0512898 -0.0515523 -0.0309709 (
        -0.0310553) -0.0310925 -0.0516153 -0.0293314 -0.0286379 -0.0492699 -0.0492005 -0.0284963 -0.028521 -0.0112 -0.0286596 -0.0492496 (
        -0.0492225) -0.0491902 -0.0285476 -0.0285399 -0.028539 -0.0286324 -0.0286365 -0.0286364 -0.0112 -0.0296962 -0.049418 -0.028672 (
        -0.0493918) -0.0493868 -0.0112 -0.0486353 -0.0482984 -0.0479474 -0.0479052 -0.0112 -0.0475018 -0.0471892 -0.028224 -0.0473702 (
        -0.0471736) -0.028224 -0.0471021 -0.0112 -0.0488216 -0.028224 -0.0485766 -0.0484914 -0.028224 -0.0112 -0.028224 2.4751 0];
        
        SolutionQ60 = [0 -0.059988 -0.0341932 -0.0341933 -0.0550303 -0.0341932 -0.0341932 -0.0377924 -0.119976 -0.0547304 -0.0546042 (
        -0.0341932) -0.0341932 -0.0341932 -0.0545414 -0.0341932 -0.0341932 -0.0548729 -0.0548634 -0.0341932 -0.0341932 -0.0135688 (
        -0.0341939) -0.0548217 -0.0547996 -0.0548942 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0341932 -0.0135688 (
        -0.0353132) -0.0554443 -0.0342743 -0.0554207 -0.0554011 -0.0135688 -0.0555326 -0.0558424 -0.0561716 -0.0561681 -0.0136579 (
        -0.0569065) -0.0574853 -0.0364788 -0.0571589 -0.0574538 -0.0363115 -0.0576154 -0.0140572 -0.0550796 -0.0360171 -0.0549047 (
        -0.0548569) -0.0341932 -0.0135688 -0.0341932 2.5207 0.143395];


        SolutionP56 = [-0.00839374 -0.03078 -0.03762 -0.04617 -0.04104 -0.029925 -0.04446 -0.047025 -0.036765 -0.050445 -0.048735 -0.041895 -0.07866 -0.04788 -0.041895 -0.045315 -0.074385 -0.04104 -0.041895 -0.045315 -0.043605 -0.038475 -0.04617 -0.047025 -0.04617 -0.04275 -0.02907 -0.043605 -0.02394 -0.04617 -0.045315 -0.045315 -0.04788 -0.048735 -0.038475 -0.036765 -0.047025 -0.038475 -0.04446 -0.04104 -0.036765 -0.04446 -0.04617 -0.04275 -0.024795 -0.045315 -0.05301 -0.038475 -0.03933 -0.04617 -0.04446 -0.045315 -0.04446 -0.040185 -0.04104 -0.047025 2.40196];
        SolutionQ56 = [ -0.00148291 0.0101169 0.0370953 0.0780664 0.0983952 -0.0120216 0.0146133 0.0463691 0.0604204 0.12224 -0.0195781 0.0137702 0.0794016 0.0800209 0.0998888 -0.0182042 0.0244492 0.0404676 0.0689796 0.108332 -0.0175172 0.0126461 0.045526 0.0783333 0.110739 -0.0171737 0.00955485 0.0429968 0.0393435 0.109903 -0.0182042 0.0148943 0.0472122 0.0808348 0.0902575 -0.0147694 0.0154564 0.0379384 0.0735685 0.0970819 -0.0147694 0.0146133 0.045526 0.0702562 0.0570481 -0.0182042 0.0174235 0.0379384 0.0646357 0.110396 -0.0178607 0.0148943 0.0438399 0.0660409 0.0970573 -0.0188911 -2.33867];
        
                #println(result);
        println("Status: $(result.status), $(termination_code(result.status))")
        println("Number of iterations: $(result.iter)")
        println("iterations : $(result.counters.k)")
        println("total time : $(result.counters.total_time) and init time $(result.counters.init_time)")
        println("time overall : $(t)")
       
        solPGPU = solution(result, pg)
        solP = Vector{Float64}(undef, N);
        copyto!(solP, solPGPU)
        
        #=println("result                        : ", solP);
        
        if N==11 || N == 12
            println("vrai reponse (cas10ba) : ", SolutionP11[indiceCompare]);
        elseif N==56 || N == 57
            println("vrai reponse (cas TestFeeder) : ", SolutionP56[indiceCompare]);
        elseif  N==50 || N == 51
            println("vrai reponse (cas69) : ", SolutionP50[indiceCompare])
        elseif  N==60 || N == 61
            println("vrai reponse (case85) : ", SolutionP60[indiceCompare])
        end=#

        solp = solution(result, p)

        #println(solp)
        solQGPU = solution(result, qg)
        solQ = Vector{Float64}(undef, N);
        copyto!(solQ, solQGPU)
        #println(solQ)
        solq = solution(result, q)
        #println(solq)
        soll = solution(result, l)
        #println(soll)
        solv = solution(result, v)
        
        if affichage
            for i=(1:B)
                print(sqrt(solv[i]), " ")
            end
            println(" ")
        end

        sum = 0;
        for g in data.agent
            sum += (0.5 * g.cost1) * solP[g.i]^2 + g.cost2 * solP[g.i] + (0.5 * g.cost1q) * solQ[g.i]^2 + g.cost2q * solQ[g.i] 
        end

        sum2 = 0;

        if N == 11 || N == 12
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP11[g.oldi]^2 + g.cost2 * SolutionP11[g.oldi] + (0.5 * g.cost1q) * SolutionQ11[g.oldi]^2 + g.cost2q * SolutionQ11[g.oldi] 
            end
        elseif N == 56 || N == 57
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP56[g.oldi]^2 + g.cost2 * SolutionP56[g.oldi] + (0.5 * g.cost1q) * SolutionQ56[g.oldi]^2 + g.cost2q * SolutionQ56[g.oldi] 
            end
        elseif N == 50 || N == 51
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP50[g.oldi]^2 + g.cost2 * SolutionP50[g.oldi] + (0.5 * g.cost1q) * SolutionQ50[g.oldi]^2 + g.cost2q * SolutionQ50[g.oldi] 
            end
        elseif N == 60 || N == 61
            for g in data.agent
                sum2 += (0.5 * g.cost1) * SolutionP60[g.oldi]^2 + g.cost2 * SolutionP60[g.oldi] + (0.5 * g.cost1q) * SolutionQ60[g.oldi]^2 + g.cost2q * SolutionQ60[g.oldi] 
            end
        end
        #println("Cost function : $(sum)")
        #println("Cost function : $(sum2)")

        return termination_code(result.status), (result.iter) , (result.counters.total_time)

        #open("result.csv", "a") do io
        #    write(io, "$(termination_code(result.status)); $(result.iter); $(result.counters.total_time); $(result.counters.init_time) \n")
        #end

    end
end

function testMultipleEndo(name, K)

    nNAgent = 7;
    nAgent = [10, 20, 50, 70, 100, 120, 150]
    M = K * 3;
    Mat = [[0.0 for i =(1:nNAgent)] for j=(1:M)]

    for i=(6:7)
        sizeAgent = nAgent[i];
        for j=(1:K)
            c, k, t = testEndoMarketDirect(name, sizeAgent);

            Mat[j][i] = c;
            Mat[j + K][i] = k;
            Mat[j + 2*K][i] = t;
        end
    end

    save_file("EndoGPU2.csv", Mat)

end


#testMatpower(filename)
#testMarket(filename)
#testEndoMarketDirect("case85", 20);
testMultipleEndo("case85", 50);


#=
# bus : ["zone", "bus_i", "bus_type", "vmax", "source_id", "area", "vmin", "index", "va", "vm", "base_kv"]
# source_type : matpower
# name : case10babis
# dcline : vide 
# gen : ["ncost", "qc1max", "pg", "model", "shutdown", "startup", "qc2max", "ramp_agc", "qg", "gen_bus",
# "pmax", "ramp_10", "vg", "mbase", "source_id", "pc2", "index", "cost", "qmax", "gen_status", "qmin", "qc1min",
# "qc2min", "pc1", "ramp_q", "ramp_30", "pmin", "apf"]  --> [gen_bus, pmax, index, qmax, qmin, ncost, cost[...] ]
# branch : ["br_r", "shift", "br_x", "g_to", "g_fr", "source_id", "b_fr", "f_bus", "br_status", "t_bus", "b_to", "index",
# "angmin", "angmax", "transformer", "tap"] -- > ["index", "f_bus", "t_bus", "br_r", "br_x", "angmin", "angmax"]
# storage : 
# switch : 
# baseMVA : 10 cas 10ba
# per_unit : false or true
# shunt  :
# load  : ["source_id", "load_bus", "status", "qd", "pd", "index"]
=#


# P et Q
function cas2bus()
    c = ExaCore()
    pmin = [0.0, -30.0, 0.0, 0.0, -1.1, -2.0]
    pmax = [0.0, 0.0, 60.0, 0.0, -0.9, 2.0]
    a = [0.0, 1.0, 1.0, 0.1, 0.1, 0.1]
    b = [0.0, 8.0, 4.0, 0.0, -0.1, 0]

    vmin = [1, 0.81]
    vmax = [1, 1.21]

    println(pmin)
    println(pmax)
    println(typeof(pmin))

    PQ = variable(c, 6, lvar = pmin, uvar = pmax); #

    #objective(c, 100 * (x[i-1]^2 - x[i])^2 + (x[i-1] - 1)^2 for i = 2:N)
    for i = 1:6
        objective(c, (0.5 * a[i] * PQ[i]^2 + b[i] * PQ[i]) )
    end

    pq = variable(c, 2);

    # grid
    x =  0.01;
    r =  0.005;

    l = variable(c, 1);
    v = variable(c, 2, lvar=vmin, uvar = vmax);

    # i = 0 
    constraint(c, (pq[1]- r*l[1]) + PQ[1] + PQ[2] - 0);
    constraint(c, (pq[2]- x*l[1]) + PQ[4] + PQ[5] - 0)

    # i = 1
    constraint(c, PQ[3] - pq[1]);
    constraint(c, PQ[6] - pq[2]);


    constraint(c, v[1] - v[2] + 2 * (r *pq[1] + x * pq[2]) - (r*r + x*x)*l[1]);
    constraint(c, l[1] * v[1] - (pq[1] * pq[1] + pq[2] * pq[2]));

    m = ExaModel(c)
    result = ipopt(m)

    println("Status: $(result.status)")
    println("Number of iterations: $(result.iter)")


    solPQ = solution(result, PQ)
    println(solPQ)
    solpq = solution(result, pq)
    println(solpq)
    soll = solution(result, l)
    println(soll)
    solv = solution(result, v)
    println(solv)

    println("contraite : " , soll[1] * solv[1] - (solpq[1] * solpq[1] + solpq[2] * solpq[2]))
end




