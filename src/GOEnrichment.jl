module GOEnrichment

using Printf
using CSV
using LightGraphs
#using LightGraphs
### Gene ontology enrichments
export parseontology
include("parse_ontology.jl")
include("annotations.jl")

ontol = parseontology("c:\\home\\resource\\geneontology\\go-basic.obo")
annot = load_annot_long("c:\\home\\resource\\Xt\\9.1\\go_b2g_xb_lf.map", ontol)
gographs = ontologygraphs(ontol)
exannot = extendannotation(annot, gographs)

rexannot = reverseannotation(exannot)
rannot = reverseannotation(annot)
end # module
