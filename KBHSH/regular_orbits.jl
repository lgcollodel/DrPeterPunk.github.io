using Plots, Glob, DelimitedFiles, Dierckx;
files = glob("s*0");
#ii = open("degenerate_prograde_orbits.dat", "r")
io = open("regular_orbits_10-30.dat", "a");
#de = readdlm(ii, Float64);
## frh_read=de
rₕ=0.01; rh="0.01";
## #=
for ifile in files
#for ifile in files[2:31]
   frh_read = readdlm(ifile, Float64);
   for line=1:size(frh_read)[1]
## =#
##   for line=1:length(de[:,1])
        nrₕ   = frh_read[line,1];
        if nrₕ!=rₕ
            rₕ = nrₕ;
            length(string(rₕ)) == 4 ? rh = string(rₕ) : rh = string(rₕ,"0");
            Lnewr = true
        else
            Lnewr = false
        end
        ωₛ   = frh_read[line,3]; lom = length(string(ωₛ));
        lom <= 8 ? om = string(ωₛ,"0"^(8-lom)) : om = string(ωₛ)[1:8];
        mass = frh_read[line,4];
        angm = frh_read[line,5];
        surg = frh_read[line,8];
#        Ldeg=false
        #=
        for lined=1:length(de[:,1])
            [rₕ,ωₛ]==de[lined,1:2] ? (Ldeg=true;break) : nothing
        end
        if !Ldeg
        =#
            sol_file = glob("rh=$rh/*$om*")
            for isol in sol_file
                println(rh,"\t",om,"\t",mass,"\t",isol)
                x, θ, mv, jv, mh, jh = charge_densities(isol,rₕ,ωₛ)
                q = jv/(jv+jh); j=angm/mass^2
                data = [rₕ,ωₛ,mass,angm,j,q,8*surg]
                Lnewr ? println(io,"\n") : nothing
                writedlm(io,[data])
            end
        #end
    end
end
#close(ii)
close(io)



# Add energy at ISCO
ii1 = open("regular_orbits_10-30.dat","r")
ii2 = open("isco_flip_10-30.dat","r")
io  = open("romajqe_regular_flip_10-30.dat","a")
read1=readdlm(ii1,Float64);
read2=readdlm(ii2,Float64);
rₕ=0.01; rh="0.01";
for line=1:size(read1)[1]
    nrₕ   = read1[line,1]
    if nrₕ!=rₕ
        rₕ = nrₕ;
        Lnewr = true
    else
        Lnewr = false
    end
    ωₛ=read1[line,2]
    mass=read1[line,3]
    angm=read1[line,4]
    j=read1[line,5]
    q=read1[line,6]
    T=read1[line,7]
    Lirreg=false
    for line2=1:length(read2[:,1])
        if [rₕ,ωₛ]==read2[line2,1:2]
            en = read2[line2,8]
            Lnewr ? println(io,"\n") : nothing
            println(rₕ,"\t",j,"\t",q,"\t",en)
            println(io,rₕ,"\t",ωₛ,"\t",mass,"\t",angm,"\t",j,"\t",q,"\t",en,"\t",T)
        end
    end
end
close(io)
close(ii1)
close(ii2)
