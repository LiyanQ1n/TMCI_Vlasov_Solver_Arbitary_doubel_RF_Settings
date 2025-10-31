push!(LOAD_PATH, "../src/")

using DoubleRFs, Documenter
using Interpolations

makedocs(
    sitename="DoubleRFs Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "主页" => "index.md",
        "示例" => "examples.md",
        "添加函数与文档编译" => "docs.md",
        "函数库" => "api.md",
    ]
)