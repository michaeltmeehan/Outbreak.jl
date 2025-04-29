struct Node
    id::Int
    left::Union{Nothing, Int}
    right::Union{Nothing, Int}
    time::Float64
end


isleaf(node::Node)::Bool = isnothing(node.left) && isnothing(node.right)