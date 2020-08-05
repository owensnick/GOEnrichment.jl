

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
    end
    gographclosed = transitiveclosure(gograph)
    (gids=gids, graphindex=graphindex, gograph=gograph, gographclosed=gographclosed)
end
