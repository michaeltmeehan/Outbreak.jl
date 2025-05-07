module BEASTPlugIns

import EzXML: readxml, ElementNode, AttributeNode, link!
import ..SeqSim


function to_elementnode(sequence::SeqSim.Sequence)
    # Convert a Sequence object to an ElementNode for XML serialization
    sequence_element = ElementNode("sequence")
    sequence_label = sequence.taxon isa AbstractString ? sequence.taxon : string(sequence.taxon)
    link!(sequence_element, AttributeNode("id", sequence_label))
    link!(sequence_element, AttributeNode("spec", "Sequence"))
    link!(sequence_element, AttributeNode("taxon", sequence_label))
    link!(sequence_element, AttributeNode("count", "4"))
    link!(sequence_element, AttributeNode("value", sequence.value))
    return sequence_element
end


function to_elementnode(alignment::SeqSim.Alignment)
    # Convert an Alignment object to an ElementNode for XML serialization
    alignment_element = ElementNode("alignment")
    link!(alignment_element, AttributeNode("id", "alignment"))
    link!(alignment_element, AttributeNode("spec", "Alignment"))
    for sequence in alignment
        sequence_element = to_elementnode(sequence)
        link!(alignment_element, sequence_element)
    end
    return alignment_element
end


export to_elementnode

end