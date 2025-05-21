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
    nothing
end



"""Calculate the change in the hamiltonian
"""
function delta_hamilton(
    cluster_edge_weights::Float64,
    vertex_weight::Float64,
    cluster_weight::Float6s4,
    gamma::Float64,
)
    out = cluster_edge_weights - vertex_weight*cluster_weight*gamma

    return out
end



"""
Implement the Constant Potts Hamiltonian from 2011 Traag et al. 
    (https://arxiv.org/pdf/1104.3083). The Hamiltonian is set as 

    \$\\mathcal{H}(P) = - \\sum_c e_c - \\gamma n_c2\$


# Constructs


##  Function Arguments

- `vec_membership`: ranged vector with denoting community elements. Unique 
    elements must be 1:max()


##  Keyword Arguments


"""
@inline function hamiltonian_constant_potts_unsafe!(
    vec_membership::Vector{Int64},
    vec_comm_size::Vector{Float64},
    vec_degrees_inv::Vector{Float64},
    gamma::Float64,
    graph::AbstractSimpleWeightedGraph;
)
    
    # verify the membership vector is appropriate
    n = nv(graph)
    m = maximum(vec_membership)

    # init--note; important to standardize comms with size 1 since the first addition in the loop should count for two
    fill!(vec_comm_size, 0)
    #vec_comm_size[1:m] .= 1

    hamiltonian = 0.0
    
    # iterate over the edges of the graph
    @inbounds for (i, e) in enumerate(edges(graph))
        
        from, to, w = src(e), dst(e), e.weight     
        
        # increase the size of the element
        vec_comm_size[vec_membership[from]] += vec_degrees_inv[from]
        vec_comm_size[vec_membership[to]] += vec_degrees_inv[to]

        (vec_membership[from] != vec_membership[to]) && continue
        
        hamiltonian += w
    end

    #vec_comm_size[1:m] = max.(vec_comm_size[1:m], 1.0)
    
    # drop the community totals 
    # Use Equation 6 in Traag 2011
    hamiltonian -= gamma*sum(vec_comm_size[1:m] .^ 2)
    # hamiltonian -= gamma*sum(binomial.(vec_comm_size[1:m], 2))
    
    return hamiltonian
end



function hamiltonian_constant_potts(
    vec_membership::Vector{Int64},
    gamma::Float64,
    graph::AbstractSimpleWeightedGraph;
)
    
    # verify the membership vector is appropriate
    n = nv(graph)
    (length(vec_membership) != n) && error("Invalid length of membership vector $(length(membership)): must have $(n) vertices.")
    
    # run unsafe herei

    return hamiltonian
end



"""
Exceute the leiden algorithm


"""
function leiden()
end



"""Implement the FastMoveNodes part of the Leiden Algorithm.

Returns a bool--change--describing whether or not the graph changed.


# Constructs

```
leiden_fast_move_nodes!(
    partition::GraphPartition1, 
    graph::SimpleWeightedGraph, # graph
    incidence::Vector,
    gamma::Float64,
    change::Bool,
)
```


##  Function Arguments

- `partition`: A GraphPartition storing membership, element size, and element
    weight
- `graph`: Graph object storing vertices, edges, and weights
- `incidence`: Graph incidence list
- `gamma`: resolution parameter
- `change`: bool that is passed 


##  Keyword Arguments


"""
function leiden_fast_move_nodes!(
    partition::GraphPartition1, 
    graph::SimpleWeightedGraph, # graph
    incidence::Vector,
    gamma::Float64,
    change::Bool,
)

    # 
    n_v, n_e = partition.graph_shape
    
    
    ##  INITIALIZE VECTORS
    
    vec_edge_weights_per_cluster = zeros(Float64, n_v, );
    vec_weights_clusters = zeros(Int64, n_v, );
    vec_stack_index = zeros(Int64, n_v, )
    
    vec_cluster_neighbors = zeros(Int64, n_v, )
    vec_cluster_neighbor_added = BitVector(ones(n_v))
    
    
    # vector indicating whether or not a vertex is stable
    vec_vertices_stable = zeros(Int64, size(partition.membership));
    queue_unstable = SimpleQueue{Int64}(
        partition.graph_shape[1]
    )
    fill_and_randomize_queue!(queue_unstable)
    
    # build a stack of empty partitions
    stack_empty_clusters = FixedStack1{Int64}(
        nv(graph); 
        index = vec_stack_index,
    )
    
    for (i, s) in enumerate(partition.element_size)
        (s > 0) && continue
        stack_push!(stack_empty_clusters, i)
    end


    ##  ITERATE OVER QUEUE

    while !isempty(queue_unstable)
        
        v = heap_pop!(queue_unstable)
        best_cluster = current_cluster = partition.membership[v]
    
        # remove the node from the current cluster and add to queue if ends up empty
        partition.element_weight[current_cluster] -= partition.vertex_weight[v]
        partition.element_size[current_cluster] -= 1
        (partition.element_size[current_cluster] == 0) && stack_push!(stack_empty_clusters, current_cluster)
    
        # attempt a "neighboring" cluster
        c = stack_top(stack_empty_clusters)
        vec_cluster_neighbors[1] = c
        vec_cluster_neighbor_added[c] = 1
        n_neighbor_clusters = 1
    
        # get edge weight
        adj_vertices = incidence[v]
        degree = length(adj_vertices)
        
        for (i, u) in enumerate(adj_vertices)
            weight = graph.weights[v, u]
    
            # get adjacent community, verify that it's been added
            c_adj = partition.membership[u]
            
            if !Bool(vec_cluster_neighbor_added[c_adj])
                n_neighbor_clusters += 1
                vec_cluster_neighbor_added[c_adj] = 1
                vec_cluster_neighbors[n_neighbor_clusters] = c_adj
            end
    
            vec_edge_weights_per_cluster[c_adj] += weight        
        end
    
    
        ##  NEXT, GET THE CLUSTER WITH BEST IMPROVEMENT
    
        max_diff = delta_hamilton(
            vec_edge_weights_per_cluster[current_cluster],
            partition.vertex_weight[v],
            partition.element_weight[current_cluster],
            gamma,
        )
        
        # iterate over neighboring clusters
        for i in 1:n_neighbor_clusters
            c_adj = vec_cluster_neighbors[i]
    
            # check the difference if swapping to c_adj
            diff = delta_hamilton(
                vec_edge_weights_per_cluster[c_adj],
                partition.vertex_weight[v],
                partition.element_weight[c_adj],
                gamma,
            )
            
            # skip if no improvement
            (diff <= max_diff) && continue
    
            best_cluster = c_adj
            max_diff = diff
            vec_edge_weights_per_cluster[c_adj] = 0
            vec_cluster_neighbor_added[c_adj] = 0
        end
    
    
        ##  UPDATE 
    
        # update cluster weights
        partition.element_weight[best_cluster] += partition.vertex_weight[v]
        partition.element_size[best_cluster] += 1
        vec_vertices_stable[v] = 1
    
        # check if best cluster is on the stack
        (stack_top(stack_empty_clusters) == best_cluster) && stack_pop!(stack_empty_clusters)
    
        # check membership and add any stable neighbors NOT in the new cluster to the queue
        if best_cluster != current_cluster
            change = true
            
            partition.membership[v] = best_cluster
    
            for u in adj_vertices
                (!Bool(vec_vertices_stable[u]) | (partition.membership[u] == best_cluster)) && continue
    
                heap_push!(queue_unstable, u)
                vec_vertices_stable[u] = 0
            end
        end
    end
    
    # clean up the indexing
    reindex_membership!(partition)

    return change
end
     


"""
Move nodes into candidate communinities.


- `graph`: Graphs.graph to operate on
- `arr_vertices`: array storing vertex membership (vertex ind is in column) 
    in a partition (row)
- `weights_nodes`: weights to place on vertices
- `gamma`: resolution parameter
- `partition`: input partition storing membership
- `changed`: whether or not move nodes has changed community alignment
"""
function mode_nodes_fast(
    graph::GGraph,
    arr_vertices::Matrix,
    weights_vertices::Vector{Int64},
    gamma::Float64,
    partition::GraphPartition1,
    changed::Bool;
    #partition::AbstractPartition;
    # queue::Union{RandomQueue{Int64}, Nothing} = nothing,
)

    # 
    vec = collect(1:20);
    Random.randcycle!(vec);


end
   


"""
Move nodes into candidate communinities.

"""
function merge_nodes_subset(
    graph::GGraph,
    partition::Vector{Int64},
    subset::Vector,
    gamma::Float64,
)

end



"""
DESCRIPTION
"""
function refine_partition(
    graph::GGraph,
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
    graph::GGraph,
)
    out = collect(vertices(graph))
    out = [[x] for x in out]

    return out
end 

#

# end the module
#end