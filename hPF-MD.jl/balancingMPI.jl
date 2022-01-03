
function SerializeArray(array::Array{Float64,2})
    out=zeros(size(array)[1]*3)
    for i in 0:size(array)[1]-1
        out[i*3+1]=array[i+1,1]
        out[i*3+2]=array[i+1,2]
        out[i*3+3]=array[i+1,3]
    end
    return out
end

function SerializeArray(array::Array{Int,2})
    out=zeros(Int64,size(array)[1]*3)
    for i in 0:size(array)[1]-1
        out[i*3+1]=array[i+1,1]
        out[i*3+2]=array[i+1,2]
        out[i*3+3]=array[i+1,3]
    end
    return out
end

function DeSerializeArray(array::Nothing)
end

function DeSerializeArray(array::Array{Float64,1})
    out=zeros(div(size(array)[1],3),3)
    for i in 0:div(size(array)[1],3)-1
        out[i+1,1]=array[i*3+1]
        out[i+1,2]=array[i*3+2]
        out[i+1,3]=array[i*3+3]
    end
    return out
end

function DeSerializeArray(array::Array{Int64,1})
    out=zeros(Int64,div(size(array)[1],3),3)
    for i in 0:div(size(array)[1],3)-1
        out[i+1,1]=array[i*3+1]
        out[i+1,2]=array[i*3+2]
        out[i+1,3]=array[i*3+3]
    end
    return out
end

function DeSerializeBondArray(array::Array{Int64,1})
    out=zeros(Int64,div(size(array)[1],3),3)
    for i in 0:div(size(array)[1],3)-1
        out[i+1,1]=array[i*3+1]
        out[i+1,2]=array[i*3+2]
        out[i+1,3]=array[i*3+3]
    end
    return out
end

function first_fit(items::Array{Int64,1},bincount::Int64)
    sort_idex = sortperm(items,rev=true)
    items = items[sort_idex] # New - improves first fit.
    bins     = [Int64[] for c =1:bincount]
    index_bins     = [Int64[] for c =1:bincount]
    binsizes = [0 for c =1:bincount]
    for (index,item) in enumerate(items)
        minbinindex = findall(x->x==minimum(binsizes),binsizes)[1]
        push!(bins[minbinindex],item)
        push!(index_bins[minbinindex],sort_idex[index])
        binsizes[minbinindex] += item
    end
    average = sum(binsizes) / bincount
    maxdeviation = maximum([abs(average - bs) for bs in binsizes])
    return bins, index_bins,binsizes, average, maxdeviation
end

function greedy_swap(columns, colsize, average, margin)
    # See if you can do a swap to smooth the heights
    a=0
end

# Target: distribute atoms on tasks as even as possible
struct BalancingMPI
    comms::Array{Array{Int64,1},1}
    mole_index_comms::Array{Array{Int64,1},1}
    commsizes::Array{Int64,1}
    average::Float64
    maxdeviation::Float64
    bond_commsizes::Array{Int64,1}
    function BalancingMPI(config::Configuration,comm_size::Int64)
        comms, mole_index_comms,commsizes, average, maxdeviation=first_fit(config.moles_length,comm_size)
        bond_commsizes=zeros(Int64,comm_size)
        if config.moles_length!=zeros(Int64,0)
            @info "CPU Loading Information: CPU_index Num_Atoms"
            for (index,tasksize_i) in enumerate(commsizes)
                @info "CPU Loading Information: $(index) $(tasksize_i)"
            end
            @info "CPU Loading Average: $(average) Atoms with Deviation: $(maxdeviation)"
            for (index,moles_comm) in enumerate(mole_index_comms)
                n_bonds=0
                for mole_i in moles_comm
                    n_bonds+=size(config.moles[mole_i].bonds)[1]
                end
                bond_commsizes[index]=n_bonds
            end
        end

        new(
            comms::Array{Array{Int64,1},1},
            mole_index_comms::Array{Array{Int64,1},1},
            commsizes::Array{Int64,1},
            average::Float64,
            maxdeviation::Float64,
            bond_commsizes
        )
    end
end