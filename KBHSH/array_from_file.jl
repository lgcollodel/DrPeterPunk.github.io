using DelimitedFiles, Dierckx, Glob;
cd("/media/Storage/Stationary_Solutions/SKBH/");

files = glob("n2gc*");
#for ifile in files[2:10]
for ifile in files[11:11]
#    frh_read = readlines(ifile);
    frh_read = readdlm(ifile, Float64);
    for line=1:size(frh_read)[1]
        rₕ   = frh_read[line,1]; rh = string(rₕ,"0");
        ωₛ   = frh_read[line,3]; om = string(ωₛ,"00000")[1:8];
        mass = frh_read[line,4];
        angm = frh_read[line,5];
        surg = frh_read[line,8];
        sol_file = glob("rh=$rh/*$om*")
        for isol in sol_file
            println(rh,"\t",om,"\t",isol)

#isol = "rh=0.08/rc=0.080_om=0.718540695300";
#isol = "rh=0.01/rb=0.01_om=0.993931185502218_h0=0.0000015608700590"
#for file in readdir("test/", join=true)
#    f_read = open(readdlm, file)


f_read =  open(readdlm,isol);
#rₕ  = 0.02;
#ωₛ  = 0.838022;
xx = f_read[:,1];
nt = length(xx);
for i=1:nt;
    if xx[i]==1.0
#        println(i)
        global nx = i
        break
    end
end
ny = Int(nt/nx);
ne = (ny-1)*nx+1;
x  = xx[ne:nt]
y  = f_read[ne:nt,2];
f0 = f_read[ne:nt,3];
f1 = f_read[ne:nt,4];
f2 = f_read[ne:nt,5];
o  = f_read[ne:nt,6];

function r(x)
    sqrt((x - 1.0)^2*rₕ^2 + x^2)/(1.0-x)
end
#r   = zeros(nx);
gtt = zeros(nx);
for i=1:nx
#      r[i] = (sqrt((x[i] - 1.0)^2*rₕ^2 + x[i]^2))/(1.0-x[i]);
    rc = r(xx[i])
    gtt[i] = -(1.0-rₕ/rc)*exp(2*f0[i])+exp(2*f2[i])*o[i]^2;
end

static_ring = false
local_minimum = false
itp = Spline1D(x, gtt);
nᵢ=1000; xᵢ=zeros(nᵢ); gttᵢ=zeros(nᵢ);
xᵢ[1]=0.0; xᵢ[2] = 1.0/(nᵢ-1); gttᵢ[1] = gtt[1]; gttᵢ[2] = itp(xᵢ[2]);
for i=3:nᵢ
    xᵢ[i]   = (i-1)/(nᵢ-1);
    gttᵢ[i] = itp(xᵢ[i]);
#    local_minimum || gttᵢ[i] < 0.0 && gttᵢ[i] > gttᵢ[i-1] && global local_minimum = true;
#    local_minimum && gttᵢ[i] < 0.0 && gttᵢ[i] < gttᵢ[i-1] && (global iₛ=i-1; static_ring=true; break)
    gttᵢ[i] < 0.0 && gttᵢ[i] < gttᵢ[i-1] && gttᵢ[i-1] > gttᵢ[i-2] && (global iₛ=i; static_ring=true; break)
#    println(gttᵢ[i])
end

if static_ring
    io = open("static_ring", "a");
    println(io,rₕ,"\t",ωₛ,"\t",mass,"\t",angm,"\t",surg,"\t",xᵢ[iₛ],"\t",r(xᵢ[iₛ]));
    close(io);
    println("Static Ring Found\n")
    static_ring = false;
else
    println("Static Ring Not Found\n")
end

end
end
end


#=
# @printf("%f \t %e \t %f \t %f \t %f \t %e \t %e \t %e \n",rh,n2gc00rh010[1,1], n2gc00rh010[1,2], n2gc00rh010[1,3], n2gc00rh010[1,4],n2gc00rh010[1,5],n2gc00rh010[1,6],sgc)
# @printf(io,"%f \t %e \t %f \t %f \t %f \t %e \t %e \t %e \n",n2gc[i,1], n2gc[i,2], n2gc[i,3], n2gc[i,4],n2gc[i,5],n2gc[i,6],n2gc[i,7],n2gc[i,8])

for i=1:184
       sgc=sg(rh,n2gc00rh010[i,5],n2gc00rh010[1,6])
       @printf("%f \t %e \t %f \t %f \t %f \t %e \t %e \t %e \n",rh,n2gc00rh010[i,1], n2gc00rh010[i,2], n2gc00rh010[i,3], n2gc00rh010[i,4],n2gc00rh010[i,5],n2gc00rh010[i,6],sgc)
       end

rh=0.20
n2gc=readdlm("n2gc00rh200",Float64)
io = open("n2gc00rh200_new", "a")
sn = size(n2gc)[1]
       for i in 1:sn
            sgc=sg(rh,n2gc[i,5],n2gc[1,6])
            @printf(io,"%f \t %e \t %f \t %f \t %f \t %e \t %e \t %e \n",rh,n2gc[i,1], n2gc[i,2], n2gc[i,3], n2gc[i,4],n2gc[i,5],n2gc[i,6],sgc)
       end
close(io)

# Another way of writing it
for i=1:nᵢ
    xᵢ[i]   = (i-1)/(nᵢ-1);
    gttᵢ[i] = itp(xᵢ[i]);
    if gttᵢ[i] < 0.0 && gttᵢ[i] < gttᵢ[i-1]
        global iₛ=i;
        println("Oops",i);
        break
    end
end
=#


#=
using Plots;
plot(xᵢ,gttᵢ)
scatter!(x,gtt)
=#
