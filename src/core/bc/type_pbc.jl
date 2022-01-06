abstract type BC end

struct NoPBC <:BC end

struct PBC <:BC
    dim::Int64
    Lx::Float64
    Ly::Float64
    Lz::Float64
    xhi::Float64
    xlo::Float64
    yhi::Float64
    ylo::Float64
    zhi::Float64
    zlo::Float64
    function PBC(   ;xhi=0.0,
                    xlo=0.0,
                    yhi=0.0,
                    ylo=0.0,
                    zhi=0.0,
                    zlo=0.0)
        Lx=xhi-xlo
        Ly=yhi-ylo
        Lz=zhi-zlo
        dim=3
        new(dim::Int64,
            Lx::Float64,
            Ly::Float64,
            Lz::Float64,
            xhi::Float64,
            xlo::Float64,
            yhi::Float64,
            ylo::Float64,
            zhi::Float64,
            zlo::Float64

        )
    end
end
