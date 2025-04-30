import Base: issorted, show, summary

abstract type AbstractNode end

struct LeafNode <: AbstractNode
    id::Int
    time::Float64
end


struct BinaryNode <: AbstractNode
    id::Int
    left::Int
    right::Int
    time::Float64
end


struct RootNode <: AbstractNode
    id::Int
    left::Union{Nothing, Int}
    right::Int
    time::Float64
end

isleaf(node::AbstractNode)::Bool = node isa LeafNode
isbinary(node::AbstractNode)::Bool = node isa BinaryNode || (node isa RootNode && !isnothing(node.left) && !isnothing(node.right))
isroot(node::AbstractNode)::Bool = node isa RootNode

ischild(node_id::Int, parent::AbstractNode)::Bool = node_id == parent.left || node_id == parent.right


function Base.summary(node::LeafNode)
    t = round(node.time; sigdigits=3)
    return "LeafNode $(node.id): (time=$t)"
end

function Base.summary(node::BinaryNode)
    t = round(node.time; sigdigits=3)
    return "BinaryNode $(node.id): (left=$(node.left), right=$(node.right), time=$t)"
end

function Base.summary(node::RootNode)
    t = round(node.time; sigdigits=3)
    lstr = isnothing(node.left) ? "∅" : string(node.left)
    return "RootNode $(node.id): (left=$lstr, right=$(node.right), time=$t)"
end


Base.show(io::IO, node::AbstractNode) = print(io, summary(node))


issorted(tree::Vector{<:AbstractNode})::Bool = Base.issorted([node.time for node in tree], rev=true)


function get_height(tree::Vector{<:AbstractNode})::Float64
    @assert issorted(tree) "Tree is not sorted in descending order of time."
    return tree[1].time - tree[end].time
end

get_root_id(tree::Vector{<:AbstractNode})::Int = [node.id for node in tree if isroot(node)][1]
get_node_ids(tree::Vector{<:AbstractNode})::Vector{Int} = [node.id for node in tree]
get_node_times(tree::Vector{<:AbstractNode})::Vector{Float64} = [node.time for node in tree]


function Base.show(io::IO, ::MIME"text/plain", tree::Vector{<:AbstractNode})
    n = length(tree)
    n_leaves = count(isleaf, tree)
    n_roots  = count(isroot, tree)
    n_binary = count(isbinary, tree)

    println(io, "Tree with $n nodes")
    println(io, "  Root nodes:   $n_roots")
    println(io, "  Binary nodes: $n_binary")
    println(io, "  Leaf nodes:   $n_leaves")

    n_roots == 1 && print(io, "  Tree height:  $(round(get_height(tree),sigdigits=3))")
end


function Base.in(node_id::Int, tree::Vector{<:AbstractNode})::Bool
    return any(node -> node.id == node_id, tree)
end


function get_branch_length(parent::AbstractNode, child_id::Int, tree::Vector{<:AbstractNode})::Float64
    @assert ischild(child_id, parent) "Node $child_id is not a child of node $(parent.id)."
    @assert child_id in tree "Node $child_id is not in the tree."
    @assert issorted(tree) "Tree is not sorted in descending order of time."
    return tree[child_id].time - parent.time
end
