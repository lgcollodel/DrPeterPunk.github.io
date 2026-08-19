using Plots, DelimitedFiles, Dierckx;

function pot(gm,E,L)
    return -1.0 .+(E^2*gm[:,3]+2*E*L*gm[:,2]+L^2*gm[:,1])./(gm[:,2].^2-gm[:,1].*gm[:,3])
end
function potL(gm,L,ϵ)
    Δ = gm[:,2].^2-gm[:,1].*gm[:,3];
    Δ[1] = 0.0
    return (-L*gm[:,2]+ϵ*sqrt.(Δ.*(gm[:,3].+L^2)))./gm[:,3]
end
function get_metric(isol)
    f_read = readdlm(isol,Float64);
    xx = f_read[:,1];
    nt = length(xx);
    for i=1:nt;
        if xx[i]==1.0
            global nx = i
            break
        end
    end
    ny = Int(nt/nx);
    ne = (ny-1)*nx+1;
    x  = xx[ne:nt]
    xp=x[1:end-1];
    y  = f_read[ne:nt,2];
    f0 = f_read[ne:nt,3];
    f1 = f_read[ne:nt,4];
    f2 = f_read[ne:nt,5];
    ω  = -f_read[ne:nt,6];

    gd   = zeros(nx-1,3);
    function r(x)
        sqrt((x - 1.0)^2*rₕ^2 + x^2)/(1.0-x)
    end
    for i=1:nx-1
        rc = r(x[i]);
        e2f2 = exp(2*f2[i]);
        gd[i,1] = -(1.0-rₕ/rc)*exp(2*f0[i])+e2f2*ω[i]^2;
        gd[i,2] = -e2f2*ω[i]*rc;
        gd[i,3] =  e2f2*rc^2;
    end
    return xp, gd
end
function get_isco_data(gm,n)
    sp1 = Spline1D(xp,gm[:,1]); sp2 = Spline1D(xp,gm[:,2]); sp3 = Spline1D(xp,gm[:,3]);
    #n=10000
    Ω=zeros(n,2)
    en   = zeros(n,2);
    an   = zeros(n,2);
    v_rr = zeros(n,2);
    xi=zeros(n);
    global Lvp  = false
    global Lvm  = false
    global Lflip = false
    global ip   = n-1
    global im   = n-1

    for i=2:n
        xi[i]=(i-1)/(n+10)
        x0 = xi[i]
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
                v_rr[i,1] = en[i,1]^2*gpp_rr+2*en[i,1]*an[i,1]*gtp_rr+an[i,1]^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
                !Lvp && v_rr[i,1] <= 0.0 && en[i,1] > 0.0 && (global Lvp=true; global ip=i);
            else
                en[i,1]=NaN; an[i,1]=NaN; v_rr[i,1]=NaN;
            end
            if Γ2>0
                en[i,2] = -(gtt+gtp*Ω[i,2])/sqrt(Γ2);
                an[i,2] = (gtp+gpp*Ω[i,2])/sqrt(Γ2);
                v_rr[i,2] = en[i,2]^2*gpp_rr+2*en[i,2]*an[i,2]*gtp_rr+an[i,2]^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
                !Lvm && v_rr[i,2] <= 0.0 && en[i,2] > 0.0 &&  (global Lvm=true; global im=i; println(i,"\t","Condition 1"));
                (Lvm && Ω[i,2] > 0.0 && im > ip) && (global Lvm=false; println(i,"\t","Condition 2"));
                (Lvm && Ω[i,2] > 0.0 && im < ip) && (global Ω[i,1] = Ω[i,2]; global ip = im; global Lflip = true; global Lvp=true; global Lvm=false; println(i,"\t","Condition 3"));
            else
                en[i,2]=NaN; an[i,2]=NaN; v_rr[i,2]=NaN;
            end
        else
            Ω[i,1]=NaN; Ω[i,2]=NaN; en[i,1]=NaN; en[i,2]=NaN; an[i,1]=NaN; an[i,2]=NaN; v_rr[i,1]=NaN; v_rr[i,2]=NaN;
        end

    #(Lvm && Lvp) ? (println(i),break) : nothing

    end
    isco_data = [xi[ip], v_rr[ip,1], Ω[ip,1], en[ip,1], an[ip,1],
                 xi[im], v_rr[im,2], Ω[im,2], en[im,2], an[im,2]]
    return  isco_data, v_rr[:,2], en[:,2], Ω[:,2]
end
# CHECK INTERPOLATION
isol1="rh=0.06/rh=0.060_om=0.692210003400";
en1=0.2735658108742825;an1=-0.2852723891434697;x1=0.19641964196419642;
isol2="rh=0.06/rh=0.060_om=0.691070002900";
en2=0.17955533015522307;an2=-2.6496788261580173;x2=0.42554255425542553;
rₕ=0.06








V=pot(gd,0.3113842777611763,1.878457201815943)


plot(x[40:70],V[40:70])

sp1 = Spline1D(x[1:end-1],gd[:,1]); sp2 = Spline1D(x[1:end-1],gd[:,2]); sp3 = Spline1D(x[1:end-1],gd[:,3]);
gtt_r=derivative(sp1,x[1:end-1]);
gtp_r=derivative(sp2,x[1:end-1]);
gpp_r=derivative(sp3,x[1:end-1]);
Δr = gtp_r.^2-gtt_r.*gpp_r;
Δr[1]=0.0
omp=(-gtp_r[2:end-10]+sqrt.(Δr[2:end-10]))./gpp_r[2:end-10]
omm=(-gtp_r[2:end-10]-sqrt.(Δr[2:end-10]))./gpp_r[2:end-10]

Δp  = -gd[2:end-10,1]-2*gd[2:end-10,2].*omp-gd[2:end-10,3].*omp.^2
Δm  = -gd[2:end-10,1]-2*gd[2:end-10,2].*omm-gd[2:end-10,3].*omm.^2

plot(xp[2:end-10],omp)
plot!(xp[2:end-10], Δp)


enp=-(gd[:,1]+gd[:,2].*omp)./sqrt.(-gd[:,1]-2*gd[:,2].*omp-gd[:,3].*omp.^2)
enm=-(gd[:,1]+gd[:,2].*omm)./sqrt(-gd[:,1]-2*gd[:,2].*omm-gd[:,3].*omm.^2)
anp=(gd[:,2]+gd[:,3].*omp)./sqrt(-gd[:,1]-2*gd[:,2].*omp-gd[:,3].*omp.^2)
anm=(gd[:,2]+gd[:,3].*omm)./sqrt(-gd[:,1]-2*gd[:,2].*omm-gd[:,3].*omm.^2)


nl = 101
vl = zeros(length(gd[:,1]),nl,2);
for i=1:nl
       li = (i-51)*10/(nl-1)
       global vl[:,i,1]=potL(gd,li,1)
       global vl[:,i,2]=potL(gd,li,-1)
end
vl_rr=zeros(length(gd[:,1]),nl);
for i=1:nl
       spv = Spline1D(xp,vl[:,i,1])
       global vl_rr[:,i]=derivative(spv,xp,nu=2)
end
plot(xp,vl[:,:,1], label=false)
plot!(xp,vl[:,:,2])


vp = -1.0 .+(en[:,1].^2 .*gd[:,3]+2*en[:,1].*an[:,1].*gd[:,2]+an[:,1].^2 .*gd[:,1])./Δ
vm = -1.0 .+(en[:,2].^2 .*gd[:,3]+2*en[:,2].*an[:,2].*gd[:,2]+an[:,2].^2 .*gd[:,1])./Δ

function ve(gm,E,L)
    Δ = gm[:,2].^2-gm[:,1].*gm[:,3];
    return -1.0 .+(E^2 *gm[:,3]+2*E*L*gm[:,2]+L^2 .*gm[:,1])./Δ
end
