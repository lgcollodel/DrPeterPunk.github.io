using DelimitedFiles, Dierckx, Glob, Plots;

cd("/media/Warehouse/Work/SKBH/");
io = open("srisco.dat", "a");
fsr_read = readdlm("static_ring",Float64);
for line=1:size(fsr_read)[1]
   rₕ   = fsr_read[line,1];
   rₕ==0.1 ? rh = string(rₕ,"0") : rh = string(rₕ);
   ωₛ   = fsr_read[line,2]; om = string(ωₛ,"00000")[1:8];
   mass = fsr_read[line,3];
   angm = fsr_read[line,4];
   surg = fsr_read[line,5];
   xsr  = fsr_read[line,6];
   rsr  = fsr_read[line,7];
   sol_file = glob("rh=$rh/*$om*")
   for isol in sol_file
       println(rh,"\t",om,"\t",isol)
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
       global ip   = nᵢ-1
       global im   = nᵢ-1



       for i=2:nᵢ-1
       #for i=nᵢ-10:-1:100
           xᵢ[i]   = (i-1)/(nᵢ-1);
           x0 = xᵢ[i]
           dx_dr   = sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*(x0 - 1.)^2/x0;
           dx_dr2  = (((2*x0 + 1.)*(x0 - 1.)^2*rₕ^2 + 2*x0^3)*(x0 - 1.))/(sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*x0^2);
           gtt     = evaluate(sp1,xᵢ[i]);
           gtt_r   = derivative(sp1,xᵢ[i])*dx_dr;
           gtt_rr  = derivative(sp1,xᵢ[i],nu=2)*dx_dr^2+gtt_r*dx_dr2;
           gtp     = evaluate(sp2,xᵢ[i]);
           gtp_r   = derivative(sp2,xᵢ[i])*dx_dr;
           gtp_rr  = derivative(sp2,xᵢ[i],nu=2)*dx_dr^2+gtp_r*dx_dr2;
           gpp     = evaluate(sp3,xᵢ[i]);
           gpp_r   = derivative(sp3,xᵢ[i])*dx_dr;
           gpp_rr  = derivative(sp3,xᵢ[i],nu=2)*dx_dr^2+gpp_r*dx_dr2;
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
           println(x0,"\t\t",v1_r,"\t\t",v2_r)
       =#
       end
#       r_data= [rₕ,ωₛ,mass,angm,surg,xsr,rsr,xᵢ[ip],r(xᵢ[ip]),v_rr[ip,1],Ω[ip,1],en[ip,1],an[ip,1],xᵢ[im],r(xᵢ[im]),v_rr[im,2],Ω[im,2],en[im,2],an[im,2]]
       r_data= [rₕ,ωₛ,mass,angm,xsr,xᵢ[ip],Ω[ip,1],xᵢ[im],Ω[im,2]]
       writedlm(io,[r_data]);
       println(ip,"\t",im)
       end
       end
       close(io);
