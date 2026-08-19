include("sDifferenceFormula.jl")
include("functions.jl")
include("sBS_funtions.jl");

# Declare Variables for the Run

Nd=1;  Nk=2;
Px=4;
M=[2,2];
Aleft=2 ; Aright=2;
Nlr=[Aleft,Aright];
# Get Initial guess
using DelimitedFiles, LinearAlgebra
#f0="h=0.102262840_ws=0.284000000000000_sigma=0.100";
f0="solitons_sol/h0=0.10122337782145846_om=0.1833_sigma=0.100";
#f0="phi=1.000000000000000_ws=0.812329116238928";
s0=readdlm(f0,Float64);
x=s0[:,1]; u0=s0[:,2:3];#u0=s0[:,4:5]
#p = [1.0,1.0,0.284,0.1]
p = [1.0,1.0,0.18329,0.1]
#p = [1.0,1.0,0.812329116238928,0.0]
potential=pot_solitonic;
#potential=pot_free;

#MAIN FUNCTION STARTS HERE
us = main(Nd,Px,x,u0,p,Nlr)


for i=1:1000
    i%100 == 0 ? println("Run n.",i,"\n ωₛ=",p[3]) : nothing
    us = main(Nd,Px,x,u0,p,Nlr)
    if us == 1
        break
    else
        u0 = us
        h0=u0[1,2]
        om=p[3]
        p[3]=p[3]-5.e-6
        i%100 == 0 ? println("h0/σ=",h0) : nothing

        if i%100 == 0
            dform=DifferenceFormula(x,Px);
            ux=zeros(dform.Nx); uxx=zeros(dform.Nx);
            dform(ux,uxx,u0[:,1])
            M=ux[end]/2
            outf=string("solitons_sol/h0=",h0,"_om=",om,"_sigma=0.100")
            io=open(outf,"a")
            for i=1:length(x)
                println(io,x[i],"\t",us[i,1],"\t",us[i,2])
            end
            close(io)
            println("M=",M)
        end
    end
end
