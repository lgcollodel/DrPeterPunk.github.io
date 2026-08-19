files = glob("00*_reg")
files2 = glob("n2gc00rh0*0")
io = open("reg_01-10.dat","a")
i=1
for file in files
       i+=1
       frh_read=readdlm(file,Float64)
       for line in 1:length(frh_read[:,1])
          nrₕ   = frh_read[line,1];
          rₕ = nrₕ;
          length(string(rₕ)) == 4 ? rh = string(rₕ) : rh = string(rₕ,"0");
          ωₛ=frh_read[line,2]
          mass = frh_read[line,3]
          angm = frh_read[line,4]
          j = frh_read[line,5]
          q = frh_read[line,6]
          en = frh_read[line,7]
          f2 = readdlm(files2[i],Float64)
          for line2 in 1:length(f2[:,1])
             o2 = f2[line2,3]
             m2 = f2[line2,4]
             T  = 8*f2[line2,8]
             if (o2==ωₛ && m2==mass)
                println(io,rₕ,"\t",ωₛ,"\t",mass,"\t",angm,"\t",j,"\t",q,"\t",en,"\t",T)
                println(rₕ,"\t",ωₛ,"\t",mass)
                break
             end
          end
       end
end
close(io)
