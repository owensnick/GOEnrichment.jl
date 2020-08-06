ontol = parseontology("c:\\home\\resource\\geneontology\\go-basic.obo")
annot = load_annot_long("c:\\home\\resource\\Xt\\9.1\\go_b2g_lf.map", ontol)
gographs = ontologygraphs(ontol)
exannot = extendannotation(annot, gographs)



### enrichment tests
function loadclusters(file="c:\\home\\projects\\pot\\gc\\pot_xen_tpm_gc_corrected_gp_smoothed.tsv")
    CSV.read(file)
end
data = loadclusters()
#
background = data.Gene
genesel = data.ClusterBothBIC .== 6
enrichtest(data.ClusterBothBIC .== 6, background, exannot, ontol) |> display

rexannot = reverseannotation(exannot)
rannot = reverseannotation(annot)
