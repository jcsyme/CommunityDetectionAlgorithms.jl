"""Partition a graph; stores information about the partition, including

    * membership vector (membership)

"""




##  KEY CLASSES

"""Partition that can be updated.


# Constructs


##  Initialization Arguments 

- `graph`: graph used to initialize the partition


##  Optional Initialization Arguments

- `element_size`: vector storing the sizes of each partition element
- `element_weight`: vector storing the aggregate weight of each partition 
    element
- `max_size`: optional specification of maximum graph size. If nothing on 
    initialization, sets to graph size.  
- `membership`: vector storing each graph vertex's partition (v_i = P_i). 
    Requires that indices are 1:n
- `membership_aux`: vector used for auxiliary membership operations
- `partition_aux`: vector used for auxiliary partition operations
- `vertex_weight`: vector storing vertex optional weights. If not specified, 
    all weights default to one

# Properties

- `element_size`::Vector{Int64}
    Vector (1:max_vertex_size) storing the size of each partition element u in
    position u for the first `partition_size` elements 
- `element_weight`::Vector{Float64}
    Vector (1:max_vertex_size) storing the size of each partition element u in
    position u for the first `partition_size` elements 
- `graph_shape`::Vector{Int64}
    Vector storing the number of vertices (n_v) and number of edges (n_e) in 
    ordered vector
- `max_vertex_size`::Int64
    Maximum size of vectors stored in partition. Only set if size could exceed 
    that of initial graph; otherwise, defaults to number of vertices n in
    initialization graph
- `membership`::Vector{Int64}
    Vector (1:max_vertex_size) storing the partition id of each vertex i in 
    position i for the first `graph_shape[1]` elements 
- `membership_aux`::Vector{Int64}
    Vector (1:max_vertex_size) used for auxiliary membership operations when 
    needed
- `partition_aux`::Vector{Int64} 
    Vector (1:max_vertex_size) used for auxiliary operations when needed
- `partition_size`::Int64 
    Number of elements in the partition
"""
struct GraphPartition1
    element_size::Vector{Int64}        # size of each partition; max size is the same size as membership
    element_weight::Vector{Float64}    # aggregate weight of the partitions (if vertex weights are none, then equal to element_size)
    graph_shape::Vector{Int64}         # number of vertices and edges associated with the partition
    max_vertex_size::Int64             # 
    membership::Vector{Int64}          # 
    membership_aux::Vector{Int64}      # used in reindex and other operations 
    partition_aux::Vector{Int64}       # int auxiliary
    partition_aux_fl::Vector{Float64}  # float auxiliary
    partition_size::Vector{Int64}      # length 1 vector storing mutable partition length
    vertex_weight::Vector{Float64}     # vector of vertex weights to store
    
    
    function GraphPartition1(
        graph::Union{AbstractGraph, SimpleWeightedGraph};
        element_size::Union{Vector{Int64}, Nothing} = nothing,
        element_weight::Union{Vector{Float64}, Nothing} = nothing,
        max_vertex_size::Union{Int64, Nothing} = nothing,
        membership::Union{Vector{Int64}, Nothing} = nothing,
        membership_aux::Union{Vector{Int64}, Nothing} = nothing,
        partition_aux::Union{Vector{Int64}, Nothing} = nothing,
        partition_aux_fl::Union{Vector{Int64}, Nothing} = nothing,
        vertex_weight::Union{Vector{Float64}, Nothing} = nothing,
    )
        
        # initialize some chacteristics from graph
        n_v, n_e = nv(graph, ), ne(graph, )
        graph_shape = [n_v, n_e]
        max_vertex_size = (
            isa(max_vertex_size, Int64) 
            ? max(n_v, max_vertex_size)
            : n_v
        )
        
        
        ##  INITIALIZE MEMBERSHIP AND ELEMENT SIZES

        # initialize the membership vector
        if isa(membership, Nothing)
            membership = zeros(Int64, max_vertex_size)
            membership[1:n_v] = collect(1:n_v)
        else
            (max_vertex_size < length(membership)) && error("Length membership cannot be less than $(max_vertex_size)")
            resize!(membership, max_vertex_size)
        end
        
        # initialize size of each element in the partition 
        element_size = initialize_vector(
            element_size,
            max_vertex_size,
            Int64,
        )

        # initialize weight of each element in the partition
        element_weight = initialize_vector(
            element_weight,
            max_vertex_size,
            Float64,
        )

        # initialize size of partition by element (auxiliary)
        partition_aux = initialize_vector(
            partition_aux,
            max_vertex_size,
            Int64,
        )

        # initialize size of partition by element (auxiliary)
        partition_aux_fl = initialize_vector(
            partition_aux_fl,
            max_vertex_size,
            Float64,
        )

        # initialize membership aux vector for faster operations
        membership_aux = initialize_vector(
            membership_aux,
            max_vertex_size,
            Int64,
        )
        
        vertex_weight = initialize_vertex_weight(
            vertex_weight,
            max_vertex_size,
            n_v,
        )
    

        # initialize partition size
        partition_size = count_partition_sizes!(
            element_size,
            element_weight,
            vertex_weight,
            membership,
            n_v,
        )

        partition_size = [partition_size]

        
        return new(
            element_size,
            element_weight,
            graph_shape,
            max_vertex_size,
            membership,
            membership_aux,
            partition_aux,
            partition_aux_fl,
            partition_size,
            vertex_weight,
        )
    end

    
end



#########################
#    START FUNCTIONS    #
#########################


"""Build a membership array (long by partition index, wide by vertex) storing 
    vertices associated with partition; i.e., element i, j is the jth vertex 
    (not necessarily ordered) associated with partition i. 

Maximum size is p*(n - p + 1) 
"""
function build_cluster_array(
    partition::GraphPartition1;
)
    m = maximum(partition.element_size)
    p = partition.partition_size
    counts = ones(Int64, p)
    
    # initialize
    arr_out = zeros(
        Int64,
        (p, m)
    )

    
    # fill in vectors with membership
    for i in 1:partition.graph_shape[1]
        # get community
        c = partition.membership[i]

        # index in vector
        ind = counts[c]
        arr_out[c, ind] = i
        counts[c] += 1
    end

    return arr_out
end



"""Build a membership array (long by partition index, wide by vertex) storing 
    vertices associated with partition; i.e., element i, j is the jth vertex 
    (not necessarily ordered) associated with partition i. 
"""
function build_cluster_dict(
    partition::GraphPartition1,
)
    m = maximum(partition.element_size)
    p = partition.partition_size
    counts = ones(Int64, p)
    
    # initialize
    dict_out = Dict(
        x => zeros(Int64, m)
        for x in 1:p
    )

    
    # fill in vectors with membership
    for i in 1:partition.graph_shape[1]
        # get community
        c = partition.membership[i]

        # index in vector
        ind = counts[c]
        dict_out[c][ind] = i
        counts[c] += 1
    end

    return dict_out
end



"""Verify a partition id.
"""
function check_parition_id!(
    c::Int64,
    n::Int64,
)
    # set some errors
    (c < 0) && error("""Invalid partition member $(c) specified; 
    memberships must be non-negative."""
    )

    (c > n) && error("""Invalid partition member $(c) specified; 
        membership specifictations cannot exceed number of vertices.
        """
    )

end



"""Count the size of each partition element individually, then sum to give the
    total number of elements in the partition (returned).
"""
function count_partition_sizes!(
    element_size::Vector{Int64},
    element_weight::Vector{Float64},
    vertex_weights::Vector{Float64},
    membership::Vector{Int64},
    n_v::Int64,
)

    # initialize element_size and get partition size
    fill!(element_size, zero(Int64))
    fill!(element_weight, zero(Float64))

    partition_size = maximum(membership)
    
    # iterate to n_v since membership can be longer than n_v
    for i in 1:n_v
        comm = membership[i]
        element_size[comm] += 1
        element_weight[comm] += vertex_weights[i]
    end

    return partition_size
end



"""Initialize the vertex weight vector based on inputs
"""
function initialize_vertex_weight(
    vertex_weight::Union{Vector{Float64}, Nothing},
    max_size::Int64,
    n_v::Int64;
    fill_value::Float64 = 1.0,
)::Vector{Float64}
    # initialize vertex weights
    vertex_weight_new = zeros(Float64, max_size)

    if isa(vertex_weight, Vector)
        (length(vertex_weight) != n_v) && error("Invalid length of vertex_weight $(length(vertex_weight)): must be length $(n_v)")
        vertex_weight_new[1:n_v] .= vertex_weight
        vertex_weight = vertex_weight_new
    else
        vertex_weight = vertex_weight_new
        vertex_weight[1:n_v] .= fill_value
    end

    return vertex_weight
end



"""Reindex membership to ensure that it is always 1:K for some K ≤ n where n is 
    the length of the membership vector.

NOTE: Based on igraph: community_misc.c (igraph_reindex_membership). O(n) 
    implementation.


# Constructs

```
reindex_membership!(
    partition::GraphPartition1;
)
```


##  Function Arguments

- `partition`: partition to reindex membership in


##  Keyword Arguments

"""
function reindex_membership!(
    partition::GraphPartition1;
)
    # verify length
    n = partition.max_vertex_size
    fill!(partition.membership_aux, 0)
    fill!(partition.partition_aux, 0)
    fill!(partition.partition_aux_fl, 0)

    ##  FIRST, ITERATE OVER MEMBERSHIP NUMBERS AND MOVE TO AUX

    ind_clusters = 1
    
    for (i, c) in enumerate(partition.membership)

        # check id first, then see if it's been defined already
        check_parition_id!(c, n)
        (partition.membership_aux[c] != 0) && continue

        # reassign once
        partition.membership_aux[c] = ind_clusters
        partition.partition_aux[ind_clusters] = partition.element_size[c]
        partition.partition_aux_fl[ind_clusters] = partition.element_weight[c]

        ind_clusters += 1
    end


    ##  NEXT, ASSIGN BACK TO ORIGINAL  

    # update element info
    partition.element_size .= partition.partition_aux
    partition.element_weight .= partition.partition_aux_fl
    fill!(partition.partition_aux, 0)
    fill!(partition.partition_aux_fl, 0)

    for (i, c) in enumerate(partition.membership)
        partition.membership[i] = partition.membership_aux[c]
    end
end



"""Update a partition with new membership vector.

# Constructs

```
update_partition!(
    partition::GraphPartition1,
    membership::Vector{Int64},
)
```


##  Function Arguments

- `partition`: Partition storing membership, element size
- `membership`: new membership vector


##  Keyword Arguments

- `reset_membership_index`: reset the membership index? If False, ensure that 
    the membership vector has the proper ordering
"""
function update_partition!(
    partition::GraphPartition1,
    membership::Vector{Int64};
    reset_membership_index::Bool = true,
)
    # verify length
    n_v = partition.graph_shape[1]
    (length(membership) != n_v) && error("Invalid membership size; should be of length $(partition.n_v)")

    # update membership
    reset_membership_index && reset_index!(membership)
    partition.membership[1:n_v] .= membership

    # update element size in place
    partition_size = count_partition_sizes!(
        partition.element_size,
        partition.element_weights,
        partition.vertex_weights,
        partition.membership,
        n_v,
    )

    partition.partition_size[1] = partition_size

    return true
end