using DelimitedFiles, Dierckx, Glob, Plots;

cd("/media/Warehouse/Work/SKBH/");
#files = glob("n2gc*");
files = open("s*0.dat", "a");

#io = open("isco.dat", "a");
#for ifile in files[22:30]
io = open("isco_flip_10-30.dat", "a")
rₕ=0.01; rh="0.01";
for ifile in files
# ifile = files[1]
    frh_read = readdlm(ifile, Float64);
#   frh_read = readdlm("om980000.temp",Float64);
#    frh_read = readdlm("isco_bh_flip.dat",Float64)
    for line=1:size(frh_read)[1]
#    for line=20:25
        nrₕ   = frh_read[line,1];
        if nrₕ!=rₕ
            rₕ = nrₕ;
            length(string(rₕ)) == 4 ? rh = string(rₕ) : rh = string(rₕ,"0");
            Lnewr = true
        else
            Lnewr = false
        end
        ωₛ   = frh_read[line,3]; om = string(ωₛ,"00000")[1:8];
#        ωₛ=0.980; om=string(ωₛ);
        mass = frh_read[line,4];
        angm = frh_read[line,5];
        surg = frh_read[line,8];
        sol_file = glob("rh=$rh/*$om*")
#        sol_file = glob("om=0.980/*$rh*")

        for isol in sol_file
            println(rh,"\t",om,"\t",mass,"\t",isol)
#=
isol="rh=0.01/rh=0.01_om=0.892500303203238_h0=0.1059955127926511"
rₕ=0.01
mass=0.617369
angm=0.369041
ωₛ=0.892500303203238
#isol="rh=0.01/rh=0.01_om=0.892765222811783_h0=0.1070357476209636"
=#

f_read = readdlm(isol,Float64);
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
ω  = -f_read[ne:nt,6];

xp=x[1:end-1];
gd   = zeros(nx-1,3);
function r(x)
    sqrt((x - 1.0)^2*rₕ^2 + x^2)/(1.0-x)
end
for i=1:nx-1
#      r[i] = (sqrt((x[i] - 1.0)^2*rₕ^2 + x[i]^2))/(1.0-x[i]);
    rc = r(x[i]);
    e2f2 = exp(2*f2[i]);
    gd[i,1] = -(1.0-rₕ/rc)*exp(2*f0[i])+e2f2*ω[i]^2;
    gd[i,2] = -e2f2*ω[i]*rc;
    gd[i,3] =  e2f2*rc^2;
end

sp1 = Spline1D(x[1:end-1],gd[:,1]); sp2 = Spline1D(x[1:end-1],gd[:,2]); sp3 = Spline1D(x[1:end-1],gd[:,3]);

#=
# test omega < 0
for i=2:nx-1
    x0 = x[i]
    dx_dr   = sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*(x0 - 1.)^2/x0;
    gtt_r   = derivative(sp1,x0)*dx_dr;
    gtp_r   = derivative(sp2,x0)*dx_dr;
    gpp_r   = derivative(sp3,x0)*dx_dr;
    Δ       = gtp_r^2-gtt_r*gpp_r;
    if Δ > 0
        Ωp = (-gtp_r+sqrt(Δ))/gpp_r;
        Ωm = (-gtp_r-sqrt(Δ))/gpp_r;
        Ωp <= 0.0 ? println(i,"\t",Ωp,"\t",Ωm) : nothing
    end
end
=#


nᵢ = 10000;
xᵢ   = zeros(nᵢ);
en   = zeros(nᵢ,2);
an   = zeros(nᵢ,2);
#er   = zeros(nᵢ,2);
#ar   = zeros(nᵢ,2);
Ω    = zeros(nᵢ,2);
v_rr = zeros(nᵢ,2);
global Lvp  = false
global Lvm  = false
global Lflip = false
global ip   = nᵢ-1
global im   = nᵢ-1


etol = 1.e-5
for i=2:nᵢ-1
#for i=nᵢ-10:-1:100
    xᵢ[i]   = (i-1)/(nᵢ-1);
    x0 = xᵢ[i]
    dx_dr   = sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*(x0 - 1.)^2/x0;
    dx_dr2  = (((2*x0 + 1.)*(x0 - 1.)^2*rₕ^2 + 2*x0^3)*(x0 - 1.))/(sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*x0^2);
    gtt     = evaluate(sp1,x0);
    gtt_r   = derivative(sp1,x0)*dx_dr;
    gtt_rr  = derivative(sp1,x0,nu=2)*dx_dr^2+gtt_r*dx_dr2;
    gtp     = evaluate(sp2,x0);
    gtp_r   = derivative(sp2,x0)*dx_dr;
    gtp_rr  = derivative(sp2,x0,nu=2)*dx_dr^2+gtp_r*dx_dr2;
    gpp     = evaluate(sp3,x0);
    gpp_r   = derivative(sp3,x0)*dx_dr;
    gpp_rr  = derivative(sp3,x0,nu=2)*dx_dr^2+gpp_r*dx_dr2;
    Δ       = gtp_r^2-gtt_r*gpp_r;
    if Δ>=0
        Ω[i,1] = (-gtp_r+sqrt(Δ))/gpp_r;
        Ω[i,2] = (-gtp_r-sqrt(Δ))/gpp_r;
        Γ1      = -gtt-2*Ω[i,1]*gtp-Ω[i,1]^2*gpp;
        Γ2      = -gtt-2*Ω[i,2]*gtp-Ω[i,2]^2*gpp;
        if Γ1>0
            en[i,1] = -(gtt+gtp*Ω[i,1])/sqrt(Γ1);
            an[i,1] = (gtp+gpp*Ω[i,1])/sqrt(Γ1);
#=
            er[i,1] = (2*(gtp*gtt_r - gtp_r*gtt)*gtp^2*gtp_r - gpp_r^2*gtt^3 - (3*gtp*gtt_r - 4*
            gtp_r*gtt)*gpp_r*gtp*gtt + 2*sqrt( - gpp_r*gtt_r + gtp_r^2)*(gpp*gtt - gtp^2)*
            (gtp*gtt_r - gtp_r*gtt) - (gtp^2*gtt_r^2 - 2*gtp*gtp_r*gtt*gtt_r + 2*gtp_r^2*
            gtt^2 - gpp_r*gtt^2*gtt_r)*gpp)/(((gpp*gtt_r - 2*gpp_r*gtt)*gtt_r - 4*(gtp*
            gtt_r - gtp_r*gtt)*gtp_r)*gpp + (4*(gtp*gtt_r - gtp_r*gtt)*gtp + gpp_r*gtt^2)*gpp_r)

            ar[i,1] = ((2*(gtp*gtt_r - gtp_r*gtt)*gtp^2*gtp_r - gpp_r^2*gtt^3 - (3*gtp*gtt_r - 4*
            gtp_r*gtt)*gpp_r*gtp*gtt + 2*sqrt( - gpp_r*gtt_r + gtp_r^2)*(gpp*gtt - gtp^2)*
            (gtp*gtt_r - gtp_r*gtt) - (gtp^2*gtt_r^2 - 2*gtp*gtp_r*gtt*gtt_r + 2*gtp_r^2*
            gtt^2 - gpp_r*gtt^2*gtt_r)*gpp)*((gpp_r*gtt - gtp*gtp_r)*gtp + (gtp*gtt_r -
            gtp_r*gtt)*gpp + sqrt( - gpp_r*gtt_r + gtp_r^2)*(gpp*gtt - gtp^2))^2)/((((gpp
            *gtt_r - 2*gpp_r*gtt)*gtt_r - 4*(gtp*gtt_r - gtp_r*gtt)*gtp_r)*gpp + (4*(gtp*
            gtt_r - gtp_r*gtt)*gtp + gpp_r*gtt^2)*gpp_r)*((gtp*gtt_r - 2*gtp_r*gtt)*gtp + gpp_r*gtt^2)^2)
=#
            v_rr[i,1] = en[i,1]^2*gpp_rr+2*en[i,1]*an[i,1]*gtp_rr+an[i,1]^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
#            !Lvp && abs(v_rr[i,1]) <= etol && en[i,1] > 0.0 && (global Lvp=true; global ip=i);
            !Lvp && v_rr[i,1] <= 0.0 && en[i,1] > 0.0 && (global Lvp=true; global ip=i);

        else
            en[i,1]=NaN; an[i,1]=NaN; v_rr[i,1]=NaN;
        end
        if Γ2>0
            en[i,2] = -(gtt+gtp*Ω[i,2])/sqrt(Γ2);
            an[i,2] = (gtp+gpp*Ω[i,2])/sqrt(Γ2);
#=
            er[i,2] = (2*(gtp*gtt_r - gtp_r*gtt)*gtp^2*gtp_r - gpp_r^2*gtt^3 - (3*gtp*gtt_r - 4*
            gtp_r*gtt)*gpp_r*gtp*gtt - 2*sqrt( - gpp_r*gtt_r + gtp_r^2)*(gpp*gtt - gtp^2)*
            (gtp*gtt_r - gtp_r*gtt) - (gtp^2*gtt_r^2 - 2*gtp*gtp_r*gtt*gtt_r + 2*gtp_r^2*
            gtt^2 - gpp_r*gtt^2*gtt_r)*gpp)/(((gpp*gtt_r - 2*gpp_r*gtt)*gtt_r - 4*(gtp*
            gtt_r - gtp_r*gtt)*gtp_r)*gpp + (4*(gtp*gtt_r - gtp_r*gtt)*gtp + gpp_r*gtt^2)*gpp_r)


            ar[i,2] =  ((2*(gtp*gtt_r - gtp_r*gtt)*gtp^2*gtp_r - gpp_r^2*gtt^3 - (3*gtp*gtt_r - 4*gtp_r*gtt)*gpp_r*gtp*gtt -
            2*sqrt( - gpp_r*gtt_r + gtp_r^2)*(gpp*gtt - gtp^2)*(gtp*gtt_r - gtp_r*gtt) - (gtp^2*gtt_r^2 -
            2*gtp*gtp_r*gtt*gtt_r + 2*gtp_r^2*gtt^2 - gpp_r*gtt^2*gtt_r)*gpp)*((gpp_r*gtt - gtp*gtp_r)*gtp +
            (gtp*gtt_r - gtp_r*gtt)*gpp - sqrt( - gpp_r*gtt_r + gtp_r^2)*(gpp*gtt - gtp^2))^2)
            /((((gpp*gtt_r - 2*gpp_r*gtt)*gtt_r - 4*(gtp*gtt_r - gtp_r*gtt)*gtp_r)*gpp + (4*(gtp*gtt_r - gtp_r*gtt)*gtp +
             gpp_r*gtt^2)*gpp_r)*((gtp*gtt_r - 2*gtp_r*gtt)*gtp + gpp_r*gtt^2)^2)
=#
            v_rr[i,2] = en[i,2]^2*gpp_rr+2*en[i,2]*an[i,2]*gtp_rr+an[i,2]^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
#            !Lvm && abs(v_rr[i,2]) <= etol && en[i,2] > 0.0 &&  (global Lvm=true; global im=i);
            !Lvm && v_rr[i,2] <= 0.0 && en[i,2] > 0.0 &&  (global Lvm=true; global im=i);

            (Ω[i,2] > 0.0 && im > ip) && global Lvm=false;
            (Ω[i,2] > 0.0 && im < ip) && (global Ω[i,1] = Ω[i,2]; global ip = im; global Lflip = true; global Lvp=true; global Lvm=false);
        else
            en[i,2]=NaN; an[i,2]=NaN; v_rr[i,2]=NaN;
        end
    else
        Ω[i,1]=NaN; Ω[i,2]=NaN; en[i,1]=NaN; en[i,2]=NaN; an[i,1]=NaN; an[i,2]=NaN; v_rr[i,1]=NaN; v_rr[i,2]=NaN;
    end

    (Lvm && Lvp) ? (println(i),break) : nothing
#    println(gtp,"\t",gtp_r,"\t",gtp_rr)
#    println(x0,"\t",gpp,"\t",gpp_r,"\t",gpp_rr)
#    println(x0,"\t\t",gtt,"\t\t",gtp,"\t\t",gpp,"\t\t",gtp_r^2,"\t\t",gtt_r*gpp_r,"\t\t",gtp_r^2-gtt_r*gpp_r)
#    println(x0,"\t\t",en[i,1],"\t\t",er[i,1],"\t\t",en[i,2],"\t\t",er[i,2])
#=
    v1   = (en[i,1]^2*gpp+2*en[i,1]*an[i,1]*gtp+an[i,1]^2*gtt)/(gtp^2-gtt*gpp)
    v2   = (en[i,2]^2*gpp+2*en[i,2]*an[i,2]*gtp+an[i,2]^2*gtt)/(gtp^2-gtt*gpp)
    v1_r = (en[i,1]^2*gpp_r+2*en[i,1]*an[i,1]*gtp_r+an[i,1]^2*gtt_r)/(gtp^2-gtt*gpp)-(en[i,1]^2*gpp+2*en[i,1]*an[i,1]*gtp+an[i,1]^2*gtt)*(2*gtp*gtp_r-gtt_r*gpp-gtt*gpp_r)/(gtp^2-gtt*gpp)^2
    v2_r = (en[i,2]^2*gpp_r+2*en[i,2]*an[i,2]*gtp_r+an[i,2]^2*gtt_r)/(gtp^2-gtt*gpp)-(en[i,2]^2*gpp+2*en[i,2]*an[i,2]*gtp+an[i,2]^2*gtt)*(2*gtp*gtp_r-gtt_r*gpp-gtt*gpp_r)/(gtp^2-gtt*gpp)^2
    println(x0,"\t\t",v1,"\t\t",v1_r)
=#
end
# println(io,rₕ,"\t",ωₛ,"\t",mass,"\t",angm,"\t",surg,"\t",xᵢ[ip],"\t",r(xᵢ[ip]),"\t",v_rr[ip,1],"\t",en[ip,1],"\t",an[ip,1],"\t",xᵢ[im],"\t",r(xᵢ[im]),"\t",v_rr[im,1],"\t",en[im,1],"\t",an[im,1])
#r_data= [rₕ,ωₛ,mass,angm,surg,xᵢ[ip],r(xᵢ[ip]),v_rr[ip,1],Ω[ip,1],en[ip,1],an[ip,1],xᵢ[im],r(xᵢ[im]),v_rr[im,2],Ω[im,2],en[im,2],an[im,2]]
Lflip ? (iflip=1; v_rr[ip,1]=v_rr[ip,2]; Ω[ip,1]=Ω[ip,2];en[ip,1]=en[ip,2];an[ip,1]=an[ip,2]) : iflip=0
r_data= [rₕ,ωₛ,mass,angm,xᵢ[ip],v_rr[ip,1],Ω[ip,1],en[ip,1],an[ip,1],xᵢ[im],v_rr[im,2],Ω[im,2],en[im,2],an[im,2],iflip]
Lnewr ? println(io,"\n") : nothing
writedlm(io,[r_data]);
println(ip,"\t",im)
end
end
end
close(io);
