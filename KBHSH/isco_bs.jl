using Plots, Glob, DelimitedFiles, Dierckx;

include("get_metric.jl")



function get_isco_data(xp,gm,n)
    sp1 = Spline1D(xp,gm[:,1]); sp2 = Spline1D(xp,gm[:,2]); sp3 = Spline1D(xp,gm[:,3]);
    ϵ=0.0
    Ω=zeros(n,2)
    en   = zeros(n,2);
    an   = zeros(n,2);
    v_rr = zeros(n,2);
    xi=zeros(n);
    global Lvp_first = false
    global Lvm_first = false
    global Ldiscont  = false
    global L2stable  = false
    global Ldomega   = false
    global Lfirst    = false
    global ip = n
    global im = n
    global Lflip = false
    for i=2:n
        global Lvp  = false
        global Lvm  = false
        xi[i]=(i-1)/(n+10)
        x0 = xi[i]
        x0 > 0.9 ? break : nothing
        dx_dr   = (x0 - 1.0)^2;
        dx_dr2  = 2*(x0 - 1.0);
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
            global Lfirst = true
            Ω[i,1] = (-gtp_r+sqrt(Δ))/gpp_r;
            Ω[i,2] = (-gtp_r-sqrt(Δ))/gpp_r;
            Γ1      = -gtt-2*Ω[i,1]*gtp-Ω[i,1]^2*gpp;
            Γ2      = -gtt-2*Ω[i,2]*gtp-Ω[i,2]^2*gpp;
            if Γ1>0
                en[i,1] = -(gtt+gtp*Ω[i,1])/sqrt(Γ1);
                an[i,1] = (gtp+gpp*Ω[i,1])/sqrt(Γ1);
                v_rr[i,1] = en[i,1]^2*gpp_rr+2*en[i,1]*an[i,1]*gtp_rr+an[i,1]^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
                if en[i,1] >= 0.0
                    if v_rr[i,1] < ϵ
                        !Lvp_first ? (global Lvp_first = true; global ip = i) : nothing
                        global Lvp = true;
                    else
                        Lvp_first ? (global Ldiscont = true ) : nothing
                    end
                elseif Lvp_first
                    global Ldiscont=true
                    #println("energy super luminal")
                end
            else
                en[i,1]=NaN; an[i,1]=NaN; v_rr[i,1]=NaN;
                Lvp_first ? (global Ldiscont=true) : nothing
            end
            if Γ2>0
                en[i,2] = -(gtt+gtp*Ω[i,2])/sqrt(Γ2);
                an[i,2] = (gtp+gpp*Ω[i,2])/sqrt(Γ2);
                v_rr[i,2] = en[i,2]^2*gpp_rr+2*en[i,2]*an[i,2]*gtp_rr+an[i,2]^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
                if (en[i,2] >= 0.0 && v_rr[i,2] <= ϵ && Ω[i,2] >= 0 )
                    if !Lvp_first
                        global Lvp_first = true
                        global ip = i
                        global Lflip = true
                    elseif Lvp
                        global L2stable = true
                    else
                        global Ldiscont = true
                    end
                elseif (en[i,2] >= 0.0 && v_rr[i,2] <= ϵ && Ω[i,2] <= 0 && !Lvm_first)
                    global Lvm_first = true
                    global im = i
                end

            else
                en[i,2]=NaN; an[i,2]=NaN; v_rr[i,2]=NaN;
            end
        else
            Ω[i,1]=NaN; Ω[i,2]=NaN; en[i,1]=NaN; en[i,2]=NaN; an[i,1]=NaN; an[i,2]=NaN; v_rr[i,1]=NaN; v_rr[i,2]=NaN;
            Lvp_first ? (global Ldiscont = true) : nothing
            Lfirst ? (global Ldomega = true; println("Discontinuity Found")) : nothing
        end
    end
    Lflip ? en[ip,1] = en[ip,2] : nothing
    return  Ldiscont, L2stable, Ldomega,Lflip, en[ip,1], en[im,2]
end




cd("/media/Warehouse/Work/SKBH/");
file="n2gc00rh000"
f_read = readdlm(file, Float64);
io = open("isco_bs.dat","a")
for line=2:size(f_read)[1]
    h0   = f_read[line,1];    sh0 = string(h0)
    ωₛ   = f_read[line,2];    lom = length(string(ωₛ));
    lom <= 8 ? om = string(ωₛ,"0"^(8-lom)) : om = string(ωₛ)[1:8];
    mass = f_read[line,3]
    angm = f_read[line,4]
    sol_file = glob("Solitons/*$h0*")
    length(sol_file) > 1 ? println("MORE THAN ONE FILE","\t",sh0,"\t",sol_file) :  println(sh0,"\t",sol_file[1])
    xp,gd,sqdetg=get_metric(sol_file[1])
    Ldiscont, L2stable, Ldomega, Lflip, ep, em = get_isco_data(xp,gd,1000)
    data = [h0, ωₛ, mass, angm, ep, em, Ldiscont, L2stable, Ldomega, Lflip]
    writedlm(io,[data])
end
close(io)
