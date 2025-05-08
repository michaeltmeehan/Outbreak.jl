module BEASTPlugIns

import EzXML: readxml, Document, ElementNode, AttributeNode, link!
import ..SeqSim

function build_elementnode(dict::Dict{T, T}) where T<:AbstractString
    # Create an ElementNode from a dictionary of attributes
    !haskey(dict, "name") && error("The key 'name' is required in the dictionary.")
    element = ElementNode(dict["name"])
    for (key, value) in dict
        link!(element, AttributeNode(key, value))
    end
    return element
end


function to_elementnode(sequence::SeqSim.Sequence)
    dict = Dict("name" => "sequence",
                "spec" => "Sequence",
                "id" => string(sequence.taxon),
                "taxon" => string(sequence.taxon),
                "count" => "4",
                "value" => sequence.value)
    return build_elementnode(dict)
end


function to_elementnode(alignment::SeqSim.Alignment; alignment_id::Union{Nothing, String}=nothing)
    dict = Dict("name" => "alignment",
                "spec" => "Alignment",
                "id" => alignment_id !== nothing ? alignment_id : "alignment")
    alignment_element = build_elementnode(dict)
    for sequence in alignment
        sequence_element = to_elementnode(sequence)
        link!(alignment_element, sequence_element)
    end
    return alignment_element
end


function make_datetrait_element(alignment::SeqSim.Alignment; alignment_id::Union{Nothing, String}=nothing)
    datetrait_element = Dict("name" => "datetrait",
                             "spec" => "base.base.evolution.tree.TraitSet",
                             "id" => "datetrait",
                             "traitname" => "date",
                             "value" => join([string(seq.taxon)*"="*string(seq.time) for seq in alignment], ",")) |> build_elementnode

    taxa_element = Dict("name" => "taxa",
                        "id" => "TaxonSet",
                        "spec" => "TaxonSet") |> build_elementnode

    link!(taxa_element, build_elementnode(Dict("name" => "alignment", "idref" => alignment_id !== nothing ? alignment_id : "alignment")))
    link!(datetrait_element, taxa_element)

    return datetrait_element
end


function make_tree_element(alignment::SeqSim.Alignment; tree_type::T="beast.base.evolution.tree.Tree", alignment_id::Union{Nothing, String}=nothing) where T<:AbstractString
    tree_element = Dict("name" => "stateNode",
                        "spec" => tree_type,
                        "id" => "Tree") |> build_elementnode

    link!(tree_element, make_datetrait_element(alignment, alignment_id=alignment_id))
    return tree_element
end


function make_beast_xml(template::Document,
                        alignment::SeqSim.Alignment;
                        tree_type::T="beast.base.evolution.tree.Tree", 
                        alignment_id::Union{Nothing, String}=nothing) where T<:AbstractString

    # Make a copy of the template document
    # doc = copy(template)

    # Create alignment element and link it to the root
    alignment_element = to_elementnode(alignment, alignment_id=alignment_id)
    link!(template.root, alignment_element)


    # Create the tree element and link it to the root element
    tree_element = make_tree_element(alignment, tree_type=tree_type, alignment_id=alignment_id)
    link!(findfirst("//run/state", template), tree_element)

    return template

end



export build_elementnode, to_elementnode, make_datetrait_element, make_tree_element, make_beast_xml

end