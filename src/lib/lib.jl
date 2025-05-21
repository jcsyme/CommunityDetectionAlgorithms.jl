"""Library functions to support generic actions.
"""



#########################
#    START FUNCTIONS    #
#########################

"""Retrieve an incidence list for a graph. Create a vector that is independent
    of changes to the graph itself.

# Constructs

```
get_incidence_list(graph::GGraph, )
```


"""
function get_incidence_list(
    graph::GGraph,
)
    inc_list = collect(
        neighbors.((graph, ), vertices(graph))
    )

    return inc_list
end




"""Reindex membership to ensure that it is always 1:K for some K ≤ n where n is 
    the length of the membership vector.

NOTE: Based on igraph: community_misc.c (igraph_reindex_membership). O(n) 
    implementation.


# Constructs

```
reindex_membership!(
    membership::Vector{Int64};
    membership_aux::Union{Vector{Int64}, Nothing} = nothing,
)
```


##  Function Arguments

- `membership`: membership vector to reindex 


##  Keyword Arguments

- `membership_aux`: optional auxiliary vector to specify for operations


"""
function reindex_membership!(
    membership::Vector{Int64};
    membership_aux::Union{Vector{Int64}, Nothing} = nothing,
)
    # verify length
    n = length(membership)
    membership_aux = initialize_vector(
        membership_aux,
        n,
        Int64,
    )

    ##  FIRST, ITERATE OVER MEMBERSHIP NUMBERS AND MOVE TO AUX

    ind_clusters = 1
    
    for (i, c) in enumerate(membership)

        # check id first, then see if it's been defined already
        check_parition_id!(c, n)
        (membership_aux[c] != 0) && continue

        membership_aux[c] = ind_clusters
        ind_clusters += 1
    end


    ##  NEXT, ASSIGN BACK TO ORIGINAL

    for (i, c) in enumerate(membership)
        membership[i] = membership_aux[c]
    end
end



