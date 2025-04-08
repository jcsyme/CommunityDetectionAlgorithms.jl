"""
Implement the Leiden Algorithm (Traag et al. 2019). See pseudocode in 
    supplementary material for souce of algorithm. The function names used 
    herein correspond directly with those denoted in the original paper.
"""
#module Leidens

using Graphs

#export aggregate_graph
#export leiden
#export mode_nodes_fast
#export merge_nodes_subset
#export refine_partition
#export singleton_partition



"""
Aggregate nodes into a community graph

"""
function aggregate_graph()
end



"""
Implement the Constant Potts Hamiltonian from 2011 Traag et al. (https://arxiv.org/pdf/1104.3083)

# Constructs


##  Function Arguments

- `vec_membership`: ranged vector with denoting community elements. Unique elements must be 1:max()


##  Keyword Arguments


"""
function hamiltonian_constant_potts(
    vec_membership::Vector{Int64},
    gamma::Float64,
    graph::AbstractSimpleWeightedGraph;
)
    
    # verify the membership vector is appropriate
    n = nv(graph)
    (length(vec_membership) != n) && error("Invalid length of membership vector $(length(membership)): must have $(n) vertices.")
    
    # init
    vec_comm_size = zeros(Int64, maximum(vec_membership))
    hamiltonian = 0.0
    
    # iterate over the edges of the graph
    for (i, e) in enumerate(edges(graph))
        
        from, to, w = src(e), dst(e), e.weight     
        
        (vec_membership[from] != vec_membership[to]) && continue
        
        # increase the size of the element
        vec_comm_size[vec_membership[from]] += 1
        hamiltonian += w
    end
    
    # drop the community totals 
    hamiltonian -= gamma*sum(vec_comm_size .^ 2)
    
    return hamiltonian
end



"""
Exceute the leiden algorithm


"""
function leiden()
end
     


"""
Move nodes into candidate communinities.

"""
function mode_nodes_fast(
    graph::AbstractGraph,
    partition::AbstractPartition;
    queue::Union{RandomQueue{Int64}, Nothing} = nothing,
)
    # 
    vec = collect(1:20);
    Random.randcycle!(vec);

end
   


"""
Move nodes into candidate communinities.

"""
function merge_nodes_subset(
    graph::AbstractGraph,
    partition::Vector{Int64},
    subset::Vector,
    gamma::Float64,
)

end



"""
DESCRIPTION
"""
function refine_partition(
    graph::AbstractGraph,
    partition::Vector,
)
    partition_refined = singleton_partition(graph)

    for (comm, ind) in enumerate(partition)
        partition_refined = merge_nodes_subset(
            graph,
            partition_refined,
            comm,
        )
    end

    return partition_refined
end 


"""
Generate an initial state, or singleton partition, for the algoritm.
"""
function singleton_partition(
    graph::AbstractGraph,
)
    out = collect(vertices(graph))
    out = [[x] for x in out]

    return out
end 

#

# end the module
end