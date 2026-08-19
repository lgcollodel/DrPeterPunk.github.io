include("get_metric.jl")
include("get_oel.jl")
include("radiant_energy.jl")



# Insert Solution, Mass and Spin parameter to get radiant energy flux for sol

function fstr(isol,rh,mass,nt)
    xp,gd,sqd = get_metric(isol,rh);
    xi,o,en,an,sq = get_oel(xp,gd,rh,sqd,nt);
    is,ie = bound(o,xi);
    re = radiant_energy(is,ie,rh,xi,sq,en,an,o);
    fs = re*mass^2;
    rs  = r.(xi[is:ie],rh)/mass;
    return rs,fs
end


# Calculate Max of Peak for each sol

function max_peak(isol,rh,mass)
    if rh != 0
        xp,gd,sqd = get_metric(isol,rh);
    else
        xp,gd,sqd = get_metric_bs(isol);
    end
    xi,o,en,an,sq = get_oel(xp,gd,rh,sqd,1001);
    is,ie = bound(o,xi);
    re = radiant_energy(is,ie,rh,xi,sq,en,an,o);
    peak, ipeak = findmax(re)
    peak *= mass^2
    rpeak = r(xi[ipeak+is-1],rh)/mass
    return rpeak, peak
end


regol = readdlm("reg_01-10.dat");
io = open("romajqetpr_01-10.dat","a");
rₕ = 0.11; rh = "0.11";
for line in 1:length(regol[:,1])
    nrₕ   = regol[line,1];
    if nrₕ!=rₕ
        rₕ = nrₕ;
        length(string(rₕ)) == 4 ? rh = string(rₕ) : rh = string(rₕ,"0");
        Lnewr = true
    else
        Lnewr = false
    end
    ωₛ   = regol[line,2]; lom = length(string(ωₛ));
    lom <= 8 ? om = string(ωₛ,"0"^(8-lom)) : om = string(ωₛ)[1:8];
    mass = regol[line,3];
    angm = regol[line,4];
    spin = regol[line,5];
    q    = regol[line,6];
    en   = regol[line,7];
    T    = regol[line,8];
    isol = glob("rh=$rh/*$om*")[1]
    println(isol)
    rp, mp = max_peak(isol, rₕ, mass)
    println(io,rₕ,"\t",ωₛ,"\t",mass,"\t",angm,"\t",spin,"\t",q,"\t",en,"\t",T,"\t",rp,"\t",mp)
end
close(io)


###################################################################################

regol = readdlm("isco_bs.dat");
io = open("romajqepr_bs.dat","a");
for line in 1:length(regol[:,1])
    h0   = regol[line,1]
    hs   = string(h0)
    ωₛ   = regol[line,2]; lom = length(string(ωₛ));
    lom <= 8 ? om = string(ωₛ,"0"^(8-lom)) : om = string(ωₛ)[1:8];
    mass = regol[line,3];
    angm = regol[line,4];
    spin = angm/mass^2;
    q    = 1;
    en   = regol[line,5];
    files = glob("Solitons/h=*$hs*")
    for isol in files
        println(isol, "\t", spin)
        rp, mp = max_peak(isol, 0.0, mass)
        println(io,h0,"\t",ωₛ,"\t",mass,"\t",angm,"\t",spin,"\t",q,"\t",en,"\t",rp,"\t",mp)
    end
end
close(io)
