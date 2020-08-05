


#### annotation
# P, F, C
# id -> vec of ids

# function load_annot_short()
#
#
#
# end

function load_annot_long(file, ontol ; genecol=:Gene, gocol=:GOID)

    table = CSV.read(file)
    annot = Dict(a => Dict{String, Vector{Int}}() for a in ["P", "F", "C"])
    namedict = Dict(biological_process => "P", molecular_function => "F", cellular_component => "C")
    missingids = 0
    annots = 0
    for (gene, goid) in zip(table[!, genecol], table[!, gocol])
        gid = goid_to_num(goid)

        if !haskey(ontol, gid)
            missingids += 1
            continue
        end
        gid = ontol[gid].id ## remove alt ids
        ns = namedict[ontol[gid].namespace]

        !haskey(annot[ns], gene)  && (annot[ns][gene] = Int[])
        if gid ∉ annot[ns][gene]
            push!(annot[ns][gene], gid)
            annots += 1
        end
    end

    println("Loaded: $(file)\t", annots, " annotations, ", missingids, " not found in ontology.");
    annot

end


reverseannotation(annot) = Dict(o => reverse_ontol_annot(ann) for (o, ann) in annot)

function reverse_ontol_annot(ann)
    r = Dict{Int, Vector{String}}()

    for (gene, ids) in ann
        for id in ids
            r[id] = push!(get(r, id, String[]), gene)
        end
    end
    r
end
