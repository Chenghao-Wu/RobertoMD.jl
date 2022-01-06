using PkgTemplates

t = Template(
    user="Chenghao Wu",
    dir=".",
    julia=v"1.5.0",
    plugins=[
        License(; name="MIT"),
        Git(),
        GitHubActions(),
        Codecov(),
        Documenter{GitHubActions}(),
    ]
)
t("hPFMD")