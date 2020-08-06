

function ontologygraphs(ontol)
    Dict("P" => ontologygraph(ontol, biological_process),
         "F" => ontologygraph(ontol, molecular_function),
         "C" => ontologygraph(ontol, cellular_component))
end


function ontologygraph(ontol, ontology=biological_process)
    gids = collect(Iterators.filter(k -> (ontol[k].namespace == ontology) && (ontol[k].id == k), keys(ontol))) ## select from namespace and confirm its not an alt id
    n = length(gids)
    graphindex = Dict(gids[i] => i for i = 1:n)

    gograph = SimpleDiGraph(n)
    for (id, gt) in ontol
        (gt.namespace != ontology) && continue
        (gt.id != id) && continue

        gi = graphindex[id]

        for j in gt.is_a
            add_edge!(gograph, gi, graphindex[j])
        end
        for j in gt.part_of
            add_edge!(gograph, gi, graphindex[j])
        end
        for j in gt.positively_regulates
            add_edge!(gograph, gi, graphindex[j])
        end
        for j in gt.negatively_regulates
            add_edge!(gograph, gi, graphindex[j])
        end
        for j in gt.regulates
            add_edge!(gograph, gi, graphindex[j])
        end
    end
    gographclosed = transitiveclosure(gograph, true)
    (gids=gids, graphindex=graphindex, gograph=gograph, gographclosed=gographclosed)
end


function extendannotation(annot, gographs)

    exannot = Dict{String, Dict{String, Vector{Int}}}()
    for (ont, ann) in annot
        exannot[ont] = Dict{String, Vector{Int}}()
        for (gid, gis) in ann
            try
                exannot[ont][gid] = goancestors(gographs[ont], gis)
            catch
                @show gid
                @show gis
                @show ont

                error("")
            end
        end
    end

    exannot
end

goancestors(gg, gis::AbstractVector)  = mapreduce(gi -> goancestors(gg, gi), vcat, gis) |> unique
goancestors(gg, gi) = gg.gids[neighbors(gg.gographclosed, gg.graphindex[gi])]


function test_ex(gene, exannot, annot)
    a = annot["P"][gene]
    e = exannot["P"][gene]
    @show setdiff(a, e)
    all([x ∈ e for x in a])

end

#test_ex("gene5966|vegt", exannot, annot)
