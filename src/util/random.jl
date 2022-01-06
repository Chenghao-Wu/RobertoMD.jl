function GaussianRandom(sys::System)
    return rand(sys.normal)
end

function RandomNumber()
    return rand()*2-1.0
end