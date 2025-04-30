struct Node
    id::Int
    left::Union{Nothing, Int}
    right::Union{Nothing, Int}
    time::Float64
end


isleaf(node::Node)::Bool = isnothing(node.left) && isnothing(node.right)

isbinary(node::Node)::Bool = !isnothing(node.left) && !isnothing(node.right)




Base.show(io::IO, node::Node) = print(io, summary(node))

function Base.summary(node::Node)
    t = round(node.time; sigdigits=3)
    if isleaf(node)
        return "Node $(node.id): (leaf, time=$(t))"
    else
        return "Node $(node.id): (left=$(node.left), right=$(node.right), time=$(t))"
    end
end


issorted(tree::Vector{Node})::Bool = Base.issorted([node.time for node in tree], rev=true)


function get_height(tree::Vector{Node})::Float64
    @assert issorted(tree) "Tree is not sorted in descending order of time."
    return tree[1].time - tree[end].time
end


function Base.show(io::IO, ::MIME"text/plain", tree::Vector{Node})
    n_nodes = length(tree)
    n_leaves = count(isleaf, tree)

    println(io, "Tree with $n_nodes nodes")
    
    if issorted(tree)
        height = round(get_height(tree); sigdigits=3)
        println(io, "  Tree height: $height")
    else
        println(io, "  ⚠️  Unsorted tree (expected descending order of time)")
    end

    print(io, "  Number of leaf nodes: $n_leaves")
end
