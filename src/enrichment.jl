
### enrichment tests
function loadclusters(file="c:\\home\\projects\\pot\\gc\\pot_xen_tpm_gc_corrected_gp_smoothed.tsv")
    CSV.read(file)
end
data = loadclusters()

background = data.Gene
for_ind = data.ClusterBothBIC .== 6


function enrichtest(for_ind, background, exannot, ontol="P")

    bgc = counter(Int)
    fgc = counter(Int)

    genes_in_bg = 0
    genes_in_fg = 0


    for (fi, bg) in zip(for_ind, background)
        !haskey(exannot[ontol], bg) && continue
        genes_in_bg += 1
        fi && (genes_in_fg += 1)
        ids = exannot[ontol][bg]
        for id in ids
            push!(bgc, id)
            fi && push!(fgc, id)
        end

        # (genes_in_fg == 10) && break
    end

    @show genes_in_fg, sum(for_ind)
    @show genes_in_bg, length(background)
    gts = collect(keys(fgc))
    A = getindex.(Ref(fgc), gts)
    B = genes_in_fg .- A
    C = getindex.(Ref(bgc), gts) .- A
    D = genes_in_bg .- genes_in_fg .- C
    DataFrame(GT=gts, A=A, B=B, C=C, D=D, FG=A, BG=getindex.(Ref(bgc), gts), OR=(A./B)./(C./D))
end

Z = enrichtest(for_ind, background, exannot)

Z[Z[:, 1] .== 7498, :]
Z[Z[:, 1] .== 7369, :]
