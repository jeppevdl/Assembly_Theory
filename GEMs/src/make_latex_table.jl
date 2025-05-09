using Pkg
if basename(pwd()) != "src"
    cd(@__DIR__)
end
Pkg.activate(".")

using CSV, DataFrames

function make_latex_table(data::String; newline_delim = "\n", newcell_delim = ",", header = true, filter = true)

    if filter # turn double / triple / ... tabs into single tabs
        multi_tab_bs = r"[\t|\n|\s]*\t[\t|\n|\s]*" 
        data = replace(data, multi_tab_bs => "\t", r"\n{2,}" => "\t", "*" => "\\*", "_" => "\\_")
    end

    data_array = split(data, newline_delim) .|> x -> split(x, newcell_delim)

    text = ""

    text *= """
    \\begin{tabular}{||$( ("c "^length(data_array[1]))[1:end-1] )||} 
    \\hline
    """

    if header
        for cell in popfirst!(data_array)
            text *= "$cell & "
        end
        text = text[1:end-2] # remove trailing `& `
        text *= "\\\\\n\\hline\n"
    end
        
    for row in data_array
        for cell in row
            text *= "$cell & "
        end
        text = text[1:end-2] # remove trailing `& `
        text *= "\\\\\n"
    end

    text *= "\\hline\n\\end{tabular}"
    print(text)
end

df = CSV.read("../bin/kegg-small/data/kegg-small.tsv", DataFrame, delim = "\t", header = true)

# If no row names, use just the data
header = join(names(df), ",")

# Create rows (no need for row names)
rows = [join(df[i, :], ",") for i in 1:nrow(df)]

# Combine into string
table_string = join([header; rows], "\n")

println(table_string)

make_latex_table(table_string)