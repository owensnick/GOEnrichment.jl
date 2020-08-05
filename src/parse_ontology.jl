

@enum Namespace biological_process=1 cellular_component=2 molecular_function=3 external=4 gene_ontology=5

mutable struct GOTerm
    id::Int
    name::String
    namespace::Namespace
    definition::String

    alt_id::Vector{Int}
    is_a::Vector{Int}

    is_obsolete::Bool
    replaced_by::Int
    consider::Vector{Int}

    part_of::Vector{Int}
    regulates::Vector{Int}
    positively_regulates::Vector{Int}
    negatively_regulates::Vector{Int}


    GOTerm() = new(0, "", gene_ontology, "", Int[], Int[], true, 0, Int[], Int[], Int[], Int[], Int[])
end

function Base.show(io::IO, got::GOTerm)
    if get(io, :compact, false)
        print(io, num_to_goid(got.id))
        print(io, "\t", got.namespace)
        print(io, "\t", got.name)
    else
        println(io, "GO Term")
        println(io, "id                   : ", num_to_goid(got.id))
        println(io, "name                 : ", got.name)
        println(io, "namespace            : ", got.namespace)
        println(io, "definition           : ", got.definition)
        println(io, "alt_id               : ", got.alt_id)
        println(io, "part_of              : ", got.part_of)
        println(io, "is_a                 : ", got.is_a)
        println(io, "is_obsolete          : ", got.is_obsolete)
        println(io, "replaced_by          : ", got.replaced_by)
        println(io, "consider             : ", got.consider)
        println(io, "regulates            : ", got.regulates)
        println(io, "positively_regulates : ", got.positively_regulates)
        println(io, "negatively_regulates : ", got.negatively_regulates)
    end
end

function parse_line(line, go_term, conds)
    for (prefix, parse_fun, field, ignore) in conds
        if startswith(line, prefix)
            if !ignore
                v = parse_fun(line[(length(prefix) + 1):end])
                if typeof(getfield(go_term, field)) <: AbstractVector
                    push!(getfield(go_term, field), v)
                else
                    try
                        setfield!(go_term, field, v)
                    catch e
                        println("Error setting: $field with $v using $parse_fun")
                        display(field)
                        display(v)
                        throw(e)
                    end
                end
            end
            return
        end
    end

    println("field missing: #$line#")
    return
end

@inline validid(go_terms, id) = haskey(go_terms, id) ? go_terms[id].id : -1
@inline isaltid(go_terms, id) = haskey(go_terms, id) && (go_terms[id].id != id)
@inline goid_to_num(id::AbstractString) = parse(Int, id[4:end])::Int
num_to_goid(num::Int) = @sprintf("GO:%07d", num)
function parse_namespace(namespace)
    if     namespace == "biological_process"
            return biological_process
    elseif namespace == "cellular_component"
            return cellular_component
    elseif namespace == "molecular_function"
            return molecular_function
    elseif namespace == "external"
            return external
    elseif namespace == "gene_ontology"
            return gene_ontology
    else
        error("Namespace parse error: $namespace")
    end
end

parse_annot_id(annot_id) = goid_to_num(annot_id[1:10])
parse_is_obsolete(is_obs) = parse(Bool, is_obs)


function get_parse_conditions()

    conds = Vector{Tuple{String, Function, Symbol, Bool}}()

    # cond tuple (prefix, parse_fun, field, ignore)
    push!(conds, ("id: ",          goid_to_num,       :id,          false))
    push!(conds, ("name: ",        identity,          :name,        false))
    push!(conds, ("namespace: ",   parse_namespace,   :namespace,   false))
    push!(conds, ("def: ",         identity,          :definition,  false))
    push!(conds, ("alt_id: ",      goid_to_num,       :alt_id,      false))
    push!(conds, ("is_a: ",        parse_annot_id,    :is_a,        false))
    push!(conds, ("is_obsolete: ", parse_is_obsolete, :is_obsolete, false))
    push!(conds, ("replaced_by: ", goid_to_num,       :replaced_by, false))
    push!(conds, ("consider: ",    goid_to_num,       :consider,    false))


    ### relationshsips
    push!(conds, ("relationship: part_of ",              parse_annot_id, :part_of,              false))
    push!(conds, ("relationship: regulates " ,            parse_annot_id, :regulates,            false))
    push!(conds, ("relationship: positively_regulates ", parse_annot_id, :positively_regulates, false))
    push!(conds, ("relationship: negatively_regulates ", parse_annot_id, :negatively_regulates, false))

    ## ignore fields
    push!(conds, ("synonym: ",  identity, :ignore, true))
    push!(conds, ("xref: ",     identity, :ignore, true))
    push!(conds, ("subset: ",   identity, :ignore, true))
    push!(conds, ("comment: ",  identity, :ignore, true))


end

function parseontology(file)
    ids = Vector{Int}()

    reading_term = false
    reading_typedef = false

    go_terms = Dict{Int, GOTerm}()
    go_term = GOTerm()

    conds = get_parse_conditions()

    for line in eachline(file)
        isempty(line) && continue
        if line == "[Term]"

            reading_term && addterm!(go_terms, go_term)
            reading_term = true
            go_term = GOTerm()
        elseif line == "[Typedef]"
            reading_term && addterm!(go_terms, go_term)
            reading_term = false
        elseif reading_term
            parse_line(line, go_term, conds)
        end
    end
    reading_term && addterm!(go_terms, go_term)
    go_terms
end

function addterm!(go_terms, go_term)
    go_terms[go_term.id] = go_term
    for gt in go_term.alt_id
        go_terms[gt] = go_term
    end
end
