struct DifferenceFormula <: Function
    q   :: Int64
    Nx  :: Int64
    Ns  :: Array{Int64,1}
    w   :: Array{Array{Float64,2},1}
    bx  :: Array{Array{Float64,1},1}
    bxx :: Array{Array{Float64,1},1}
    s   :: Array{Array{Int64,1},1}

    function DifferenceFormula(x,q)
        Nx = length(x)
        Ns  = zeros(Int64,Nx)
        ss   = Array{Int64}[]
        sw   = Array{Float64,2}[]
        sbx  = Array{Float64}[]
        sbxx = Array{Float64}[]

        for i0=1:Nx
            # Initialize local arrays
            w         = zeros(q+2,q+2)
            bx        = zeros(q+2)
            bxx       = zeros(q+2)
            bx[1] = 0.0; bx[2] = 1.0;
            bxx[1] = 0.0; bxx[2] = 0.0; bxx[3] = 2.0;

            x0=x[i0]
            # Fill In the Stencil Vector
            si0=[i0]
            # Fully Centered
            if i0 > q/2 && (Nx-i0) >= q/2
                for i=1:Int(q/2)
                    push!(si0,i0+i)
                    push!(si0,i0-i)
                end
            # Skewed to the right
            elseif i0 <= q/2
                for i=1:i0-1
                    #push!(s,i0+i)
                    push!(si0,i0-i)
                end
                for i=1:Int(q+2-i0)
                    push!(si0,i0+i)
                end
            # Skewed to the left
            elseif (Nx-i0) < q/2
                for i=1:(Nx-i0)
                    push!(si0,i0+i)
                end
                for i=1:Int(q+1-Nx+i0)
                    push!(si0,i0-i)
                end
            end

            Ns[i0]=length(si0)


        # Derivative of Basis Function
            for j=3:Ns[i0]
                bx[j]=0.0
                for k=1:j-1
                    fact=1.0
                    for i=1:j-1
                        if i != k
                            fact *= (x0-x[si0[i]])
                        end
                    end
                    bx[j]+=fact
                end
            end

        # Second Derivative Basis Function
            for j=4:Ns[i0]
                bxx[j]=0.0
                for l=1:j-1
                    for k=1:j-1
                        if k != l
                            fact=1.0
                            for i=1:j-1
                                if i != k && i != l
                                    fact *= (x0-x[si0[i]])
                                end
                            end
                            bxx[j]+=fact
                        end
                    end
                end
            end

        # Divided Difference
            for i=1:Ns[i0]
                for j=1:Ns[i0]
                    w[i,j]=1.0
                    for p=1:j
                        if p != i
                            w[i,j] = w[i,j]/(x[si0[i]]-x[si0[p]])
                        end
                    end
                end
            end
            push!(sbx,bx)
            push!(ss,si0)
            push!(sbxx,bxx)
            push!(sw,w)
        end
        new(q,Nx,Ns,sw,sbx,sbxx,ss)
    end
end



# For a regular grid, the stencil and basis functions only need to be
# computed once in each direction

function (derivatives!::DifferenceFormula)(Mx,Mxx)
    for i0=1:derivatives!.Nx
        for j=1:derivatives!.Ns[i0]
            for i=1:j
                j0=derivatives!.s[i0][i]
                Mx[i0,j0]+=derivatives!.w[i0][i,j]*derivatives!.bx[i0][j]
                Mxx[i0,j0]+=derivatives!.w[i0][i,j]*derivatives!.bxx[i0][j]
            end
        end
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
    end
end


################################################################################
################################################################################
# 1 D
# Structure for the parameters, size and pre-allocated arrays
abstract type Parameter <: Integer end
struct ArraySize <: Parameter
    Nd :: Int64
    Nk :: Int64
    Nx :: Array{Int64,1}
end
mutable struct SolutionVector
    u0 :: Array{Float64,2}
    u1 :: Array{Float64,2}
end
mutable struct DerivArrays1D <: Function
    ux  :: Array{Float64,2}
    uxx :: Array{Float64,2}
    function DerivArrays1D(params::ArraySize)
        Nk=params.Nk
        Nx=params.Nx[1]
        self=new()
        self.ux  = Array{Float64,2}(undef,(Nx,Nk))
        self.uxx = Array{Float64,2}(undef,(Nx,Nk))
        return self
    end
end
mutable struct DerivMatrices1D <: Function
    dx  :: Array{Float64,2}
    dxx :: Array{Float64,2}
    function DerivMatrices1D(params::ArraySize)
        Nx=params.Nx[1]
        self=new()
        self.dx  = zeros(Nx,Nx)
        self.dxx = zeros(Nx,Nx)
        return self
    end
end
mutable struct EqArray1D <: Function
    Pu  :: Array{Float64,1}
    Res :: Array{Float64,1}
    function EqArray1D(params::ArraySize)
        Nk = params.Nk[1]
        Nx = params.Nx[1]
        self=new()
        self.Pu  = Array{Float64,1}(undef,Nk)
        self.Res = Array{Float64,1}(undef,Nk*Nx)
        return self
    end
end
mutable struct JacMatrices1D <: Function
    Ju   :: Array{Float64,3}
    Jux  :: Array{Float64,3}
    Juxx :: Array{Float64,3}
    Jac  :: Array{Float64,2}
    function JacMatrices1D(params::ArraySize)
        Nx = params.Nx[1]
        Nk = params.Nk[1]
        self=new()
        self.Ju   = zeros(Nx,Nk,Nk)
        self.Jux  = zeros(Nx,Nk,Nk)
        self.Juxx = zeros(Nx,Nk,Nk)
        self.Jac  = zeros(Nx*Nk,Nx*Nk)
        return self
    end
end

################################################################################
################################################################################
# 2 D
# Structure for the parameters, size and pre-allocated arrays
mutable struct DerivArrays2D <: Function
    ux  :: Array{Float64,3}
    uxx :: Array{Float64,3}
    function DerivArrays2D(params::ArraySize)
        Nk=params.Nk
        Nx=params.Nx
        self=new()
        self.ux=Array{Float64,3}(undef,(Nx[1],Nx[2],Nk))
        self.uxx=Array{Float64,3}(undef,(Nx[1],Nx[2],Nk))
        return self
    end
end
mutable struct DerivArrays3D <: Function
    ux  :: Array{Float64,4}
    uxx :: Array{Float64,4}
    function DerivArrays3D(params::ArraySize)
        Nk=params.Nk
        Nx=params.Nx
        self=new()
        self.ux=Array{Float64,3}(undef,(Nx[1],Nx[2],Nx[3],Nk))
        self.uxx=Array{Float64,3}(undef,(Nx[1],Nx[2],Nx[3],Nk))
        return self
    end
end


######
function newvaluate(x,u,q,x0)
    lx  = length(x)
    dx0 = x[end]-x[1]
    i0  = lx
    for i=1:lx
        dx = abs(x[i]-x0)
        dx < dx0 ? (dx0=dx; i0=i) : nothing
        dx > dx0 ? break : nothing
    end

    s=[i0]
    if i0 > q/2 && (lx-i0) > q/2
        for i=1:Int(q/2)
            push!(s,i0+i)
            push!(s,i0-i)
        end
    elseif i0 < q/2 + 1
        for i=1:i0-1
            push!(s,i0+i)
            push!(s,i0-i)
        end
        for i=i0:Int(q/2+2)
            push!(s,i0+i)
        end
    elseif (lx-i0) < q/2 + 1
        for i=1:(lx-i0)
            push!(s,i0+i)
            push!(s,i0-i)
        end
        for i=(lx-i0+1):Int(q/2+1)
            push!(s,i0-i)
        end
    end

    ls = length(s)
    b    = zeros(ls)
    b[1] = 1.0
    for j=2:ls
        b[j] = b[j-1]*(x0 - x[s[j-1]])
    end

    w = zeros(ls,ls)
    for i=1:ls
        for j=1:ls
            w[i,j]=1.0
            for p=1:j
                if p != i
                    w[i,j] = w[i,j]/(x[s[i]]-x[s[p]])
                end
            end
        end
    end


    u0=0.0
    for j=1:ls
        for i=1:j
            u0   += u[s[i]]*w[i,j]*b[j]
        end
    end

    return u0
end
