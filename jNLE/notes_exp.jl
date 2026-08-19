#=
NOTE:
IN JULIA, ON CANNOT CHANGE ROWS OR COLUMNS OF A MULTIDIMENSIONAL ARRAY IN PLACE BY SENDING JUST
THIS ROW OR COLUMN (LIKE A[:,1]). EITHER ONE SENDS THE WHOLE ARRAY, OR A SUBARRAY WITH THE MACRO @VIEW.
=#

function (derivatives!::DifferenceFormula)(vx::Array{Float64,1},vxx::Array{Float64,1},u)
    for i0=1:derivatives!.Nx
    # Compute ux(x0), uxx(x0)
        u0x=0.0;u0xx=0;
        for j=1:derivatives!.Ns[i0]
            for i=1:j
                u0x  += u[derivatives!.s[i0][i]]*derivatives!.w[i0][i,j]*derivatives!.bx[i0][j]
                u0xx  += u[derivatives!.s[i0][i]]*derivatives!.w[i0][i,j]*derivatives!.bxx[i0][j]
            end
        end
        vx[i0] = u0x
        vxx[i0] = u0xx
        println(vx[i0])
    end
end
function (derivatives!::DifferenceFormula)(vx,vxx,u)
    for i0=1:derivatives!.Nx
    # Compute ux(x0), uxx(x0)
        u0x=0.0;u0xx=0;
        for j=1:derivatives!.Ns[i0]
            for i=1:j
                u0x  += u[derivatives!.s[i0][i]]*derivatives!.w[i0][i,j]*derivatives!.bx[i0][j]
                u0xx  += u[derivatives!.s[i0][i]]*derivatives!.w[i0][i,j]*derivatives!.bxx[i0][j]
            end
        end
        vx[i0] = u0x
        vxx[i0] = u0xx
        println(vx[i0])
    end
end
vux=zeros(dform.Nx,Nk);vuxx=zeros(dform.Nx,Nk);
# Doesn't Work
for k=1:Nk
    dform(vux[:,k],vuxx[:,k],u0[:,k])
end
# With Subarray
for k=1:Nk
    sax=@view vu.ux[:,k];saxx=@view vu.uxx[:,k]
    dform(sax,saxx,u0[:,k])
end




function changea!(vx::Array{Float64,1},a::Float64)
    for i=1:length(vx)
        vx[i]=(i-1)*a
    end
    return vx
end

function changea!(vx::Array{Float64,2},a::Float64)
    for j=1:length(vx[1,:])
        for i=1:length(vx[:,1])
            vx[i,j]=(i-1)*a
        end
    end
end


mutable struct CA
    a :: Float64
end
function changef!(c::Float64)
    c=3.0
end
