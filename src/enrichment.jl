
### enrichment tests
# function loadclusters(file="c:\\home\\projects\\pot\\gc\\pot_xen_tpm_gc_corrected_gp_smoothed.tsv")
#     CSV.read(file)
# end
# data = loadclusters()
#
# background = data.Gene
# for_ind = data.ClusterBothBIC .== 6


function enrichtest(genesel, background, exannot, ontol, ontology="P", minannot=2, pthresh=0.1)

    bgc = counter(Int)
    fgc = counter(Int)

    genes_in_bg = 0
    genes_in_fg = 0

    id_gene = Dict{Int, Vector{String}}()

    for (fi, bg) in zip(genesel, background)
        !haskey(exannot[ontology], bg) && continue
        genes_in_bg += 1
        fi && (genes_in_fg += 1)
        ids = exannot[ontology][bg]
        for id in ids
            push!(bgc, id)
            if fi
                push!(fgc, id)
                # !haskey(id_gene, id) && (id_gene[id] = String[])
                # push!(id_gene[id], bg)
            end
        end
    end


    gts = collect(keys(fgc))
    A = getindex.(Ref(fgc), gts)
    B = genes_in_fg .- A
    C = getindex.(Ref(bgc), gts) .- A
    D = genes_in_bg .- genes_in_fg .- C

    filtind = (A .>= minannot) .& .!iszero.(B) .& .!iszero.(C) .& .!iszero.(D)

    FT = FisherExactTest.(A[filtind], B[filtind], C[filtind], D[filtind])
    OR = [ft.ω for ft in FT]
    PFT = pvalue.(FT, tail=:right)

    fgts = gts[filtind]

    df = DataFrame(GT=num_to_goid.(fgts), Term=[ontol[t].name for t in fgts], FG=A[filtind], BG=getindex.(Ref(bgc), fgts), OR=OR, PFT=PFT)
    re = reverseannotation(exannot)
    df[!, :Genes] = [join(reducedgenename.(intersect(re[ontology][gid], background[genesel])), ", ") for gid in fgts]

    sdf = sort(df[df.PFT .< pthresh, :], :PFT)

    sdf
end


function reducedgenename(g)
    fields = split(g, "|")
    if length(fields) == 1
        first(fields)
    else
        last(fields)
    end
end
