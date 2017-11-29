using JuMP, Gurobi, Cbc

function degerry(
    votes,
    contiguity_matrix, 
    number_districts,
    common_size_threshold = 0.2; 
    solver = "Cbc"
    )
    
    _V = votes
    _C = contiguity_matrix

    blocks = size(_V,1)
    districts = number_districts
    total_vote = _V * ones(2,1)

    # Do some checks
    if any(_C != _C')
        throw(ArgumentError("Contiguity matrix is not valid. It must be symmetric."))
    end
    
    if !(solver in ["Cbc","Gurobi"])
        throw(ArgumentError(string(solver, " is not a valid solver choice. Must be either Cbc or Gurobi")))
    end
    
    if solver == "Cbc"
        m = Model(solver = CbcSolver())
    else
        m = Model(solver = GurobiSolver(Presolve=0))
    end
    
    ## Variables

    @variable(m, 0 <= D[i=1:blocks,j=1:districts] <= 1 , Bin)
    
    ## Constraints  

    # each block can be in only one district
    @constraint(m, D * ones(districts,1) .== 1)  
    
    # Each district must have at least one block
    # @constraint(m, (D' * V) * [1;1] .>= 1)

    # These constraints set wasted_u to the number of wasted votes for the losing party
    @variable(m, 0 <= w[i=1:districts] <= 1, Bin)
    @variable(m, wasted_u[i=1:districts, j=1:2])
    M = blocks * sum(total_vote) 
    @constraint(m, wasted_u .>= 0)
    @constraint(m, wasted_u[:,1] .>= (D' * _V)[:,1] - M * w)
    @constraint(m, wasted_u[:,2] .>= (D' * _V)[:,2] - M * (1-w))

    # These constraints set wasted_o to the number of wasted votes for the winning party
    @variable(m, wasted_o[i=1:districts, j=1:2])
    @variable(m, votes_to_win[i=1:districts])
    @constraint(m, votes_to_win .== (D' * _V) * [1;1] / 2)
    @constraint(m, wasted_o .>= 0)
    @constraint(m, wasted_o .>= (D' * _V) - votes_to_win * [1 1])

    # These constraints calculate the efficiency gap
    @variable(m, eff_gap)
    @variable(m, abs_eff_gap)
    @constraint(m, eff_gap .== ones(1,districts) * (wasted_u + wasted_o) * [1;-1])
    @constraint(m, abs_eff_gap >= eff_gap)
    @constraint(m, abs_eff_gap >= - eff_gap)

    # These constraints enforce roughly equal sizes. 
    @variable(m, common_size) # this approach is too slow
    fixed_common_size = sum(_V) / districts
    # @constraint(m, (D' * _V) * [1;1] .>= fixed_common_size * (1-common_size_threshold)) # we don't really need this
    @constraint(m, (D' * _V) * [1;1] .<= fixed_common_size * (1+common_size_threshold))

    # These constraints enforce contiguity, but we need to allow districts with only one block
    @variable(m, 0 <= multi_block_districts[i=1:districts] <= 1, Bin)
    @constraint(m, multi_block_districts .>= 0 )
    @constraint(m, M * multi_block_districts' .>= ones(1,blocks) * D - ones(1,districts) )
    # @constraint(m, _C * D .>= 2 * D) this by itself is not enough
    @constraint(m, _C * D .>= 2 * D - M * (1-repmat(multi_block_districts',blocks)))
    

    ## Objective

    @objective(m, Min, abs_eff_gap  + sum(wasted_u) + sum(wasted_o) + sum(multi_block_districts) ) 
    
    @time begin
        status = solve(m)
    end
    
    res = Dict([("Model",m),
        ("Solve Status", status), 
        ("Efficiency Gap", getvalue(abs_eff_gap) ),
        ("Wasted Over Votes", getvalue(wasted_o)),
        ("Wasted Under Votes", getvalue(wasted_u)),
        ("Total Wasted Votes [D R]", ones(1,districts) * ( getvalue(wasted_u) + getvalue(wasted_o))),
        ("Votes By District", getvalue(D)' * _V), 
        ("Common Size", getvalue(common_size)), 
        ("Fixed Common Size", fixed_common_size), 
        ("District Assignments", getvalue(D)), 
        ("Total Vote Share", sum(getvalue(D)' * _V,1) ), 
        ("Total Seat Share", sum( getvalue(D)' * _V .>= repmat(maximum((getvalue(D)' * _V),2),1,2), 1)  ), 
        ("Number of blocks in each district", sum(getvalue(D),1) ), 
        ("Total votes in each district", getvalue(D)' * _V * [1;1] )
    ])
    
    return res
end

# Load data
using CSV
WI_votes = CSV.read("/home/tpatricksullivan/gerry/data/Gerrymander County_election_data.csv")
WI_contiguity = CSV.read("/home/tpatricksullivan/gerry/data/Gerrymander County_contiguity.csv", rows = 73)

WI_V = convert(Array, WI_votes[:,3:4])
WI_C = convert(Array, WI_contiguity[:,2:73])
