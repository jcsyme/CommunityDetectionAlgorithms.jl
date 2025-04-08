# TEMP
"""
Initialize a vector or Matrix. NOTE: if passing a matrix, then the operation
    is unsafe; it cannot be resized.

    This
    
##  Constructs

```
initialize_vector(
    vec::Union{Vector{T}, Nothing},
    size::Int64,
    U::DataType = Float64;
    allow_sz_geq::Bool = true,
    fill_func::Function = zero,
) where {T<:Real}
```

```
initialize_vector(
    vec::Nothing,
    len::Int64,
    U::DataType = Float64;
    allow_sz_geq::Bool = true,
    coerce_vector::Bool = false,
    fill_func::Function = zero,
) where T<:Real
```


##  Function Arguments

- `vec`: Vector object to initialize. If nothing, returns a new vector size
    `len` with element types `U`
- `len`: length of the vector
- `U`: (only required if vec is Nothing) element type for the Vector


##  Keyword Arguments

- `allow_sz_geq`: if True, allows a vector that's passed that's greater than or
    equal than len to be preserved. Most cases allow this to be true
- `coerce_vector`; set to true to force the input is coerced to a vector
- `fill_func`: 
    * `one` to fill with ones
    * `typemax` to set using maximum value for type `T`
    * `zero` to fill with zeros
"""
function initialize_vector(
    vec::Nothing,
    len::Int64,
    U::DataType = Float64;
    allow_sz_geq::Bool = true, # does nothing in this case
    coerce_vector::Bool = false,
    fill_func::Function = zero,
)

    vec = zeros(U, len)
    (fill_func == zero) && (return vec)

    fill!(vec, fill_func(U))
    return vec
end


function initialize_vector(
    vec::Vector{T},
    len::Int64,
    U::DataType = T;
    allow_sz_geq::Bool = true,
    coerce_vector::Bool = false,
    fill_func::Function = zero,
) where T<:Real

    # do we need to resize the vector?
    sz = length(vec)
    resize_q = allow_sz_geq ? (len > sz) : (len != sz)

    resize_q && resize!(vec, len)
    fill!(vec, fill_func(eltype(vec)))

    return vec
end



function initialize_vector(
    vec::Matrix{T},
    len::Int64,
    U::DataType = T;
    allow_sz_geq::Bool = true,
    coerce_vector::Bool = false,
    fill_func::Function = zero,
) where T<:Real

    # convert matrix or subarray to vector?
    (!isa(vec, Vector) & coerce_vector) && (vec = vec[:, 1])
    sz = size(vec)
    

    # throw an error if there's a mismatch
    error_row = allow_sz_geq ? (len > sz[1]) : (len != sz[1])
    error_row && error("Invalid size of vec in initialize_vector(); nrow(vec) == $(sz[1]), not the specified size $(len)")
    
    (sz[2] != 1) && error("Invalid size of vec in initialize_vector(); ncol(vec) != 1. If a matrix is passed, it must have only one column.")
    
    fill!(vec, fill_func(eltype(vec)))

    return vec
end



function initialize_vector(
    vec::SubArray{T},
    len::Int64,
    U::DataType = T;
    allow_sz_geq::Bool = true,
    coerce_vector::Bool = false,
    fill_func::Function = zero,
) where T<:Real

    sz = size(vec)

    # throw an error if there's a mismatch
    error_row = allow_sz_geq ? (len <= sz[1]) : (len != sz[1])
    error_row && error("Invalid size of vec in initialize_vector(); nrow(vec) == $(sz[1]), not the specified size $(sz)")
    
    if length(sz) > 1
        (sz[2] != 1) && error("Invalid size of vec in initialize_vector(); ncol(vec) != 1. If a matrix is passed, it must have only one column.")
    end

    # THIS WILL BREAK, BUT THAT'S OK--NEED TO FIX BEHAVIOR
    fill!(vec, fill_func(eltype(vec)))

    return vec
end
