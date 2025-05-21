"""Fixed length stack.
"""
struct FixedStack1{T<:Real}
    empty_value::T
    index::Vector{T}
    max_size::Int64
    size::Vector{Int64}
    
    function FixedStack1{T}(
        max_size::Int64;
        empty_value::T = zero(T),
        index::Union{Vector{T}, Nothing} = nothing,
    ) where {T}

        # check that max size is properly specified
        (max_size <= 0) && error("Stack must be able to hold more than 0 elements.")
        
        # initialize the index--if it already exists, just fills it
        index = initialize_vector(
            index,
            max_size,
            T,
        )

        size = [0]
        
        return new{T}(
            empty_value,
            index,
            max_size,
            size,
        )
    end
end



"""Clear a stack.

# Constructs

```
clear!(stack::FixedStack1{T}, element::T; )
```

"""
function clear!(
    stack::FixedStack1,
)
    fill!(stack.index, stack.empty_value)
    stack.size[1] = 0
end



"""Push to a fixed-length stack. Can push individual elements or a vector. 

# Constructs

```
stack_push!(stack::FixedStack1{T}, element::T; )
```

```
stack_push!(stack::FixedStack1{T}, elements::Vector{T}; )
```
"""
function stack_push!(
    stack::FixedStack1{T},
    element::T,
) where {T}
    # check size
    if (stack.size[1] == stack.max_size)
        @info("Stack is full with $(stack.max_size) elements. Cannot push.")
        return nothing
    end

    ind = stack.size[1] + 1
    stack.index[ind] = element
    stack.size[1] += 1;
end

function stack_push!(
    stack::FixedStack1{T},
    elements::Vector{T},
) where {T}
    # check size
    if (stack.size[1] == stack.max_size)
        @info("Stack is full with $(stack.max_size) elements. Cannot push.")
        return nothing
    end

    ind_0 = stack.size[1] + 1
    ind_1 = min(stack.size[1] + length(elements), stack.max_size)
    n = ind_1 - ind_0 + 1
    
    stack.index[ind_0:ind_1] .= elements[1:n]
    stack.size[1] += n; 
end



"""Retrieve the last element pushed to the stack.

# Constructs

```
stack_pop!(stack::FixedStack1{T}; )
```

"""
function stack_pop!(
    stack::FixedStack1{T},
) where {T}
    # check size
    (stack.size[1] == 0) && (return nothing)

    # get the size
    ind = stack.size[1]
    out = stack.index[ind]
    stack.index[ind] = stack.empty_value
    stack.size[1] -= 1

    return out
    
end



"""Show the top element of the stack.

# Constructs

```
stack_top!(stack::FixedStack1{T}; )
```

"""
function stack_top(
    stack::FixedStack1{T},
) where {T}
    # check size
    (stack.size[1] == 0) && (return nothing)

    # get the size
    value = stack.index[stack.size[1]]

    return value
    
end