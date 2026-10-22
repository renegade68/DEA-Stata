# DEA-Stata

Stata software for the article:

> Lee, C. Forthcoming. *dea: An Update and Extension to Data Envelopment Analysis*.
> Accepted for publication in the [Stata Journal](https://www.stata-journal.com/).

This repository maintains the article's DEA command family and subsequent updates.
The current development version may differ from the software submitted with the
article. The acceptance notice refers to the article; it does not extend the
article's validation claims to later code changes.

Requires Stata 16 or newer. The command index below defines the software scope.
The game-theoretic commands are exploratory diagnostics, as described in the
article, and are not fully specified structural game estimators. Command-specific
assumptions and limitations are documented in the linked help files.

## Start using the commands

Download or clone this repository. Keep **all files in `ado/` together**, including
the two implementation `.do` files. In Stata, add that folder for this session:

```stata
adopath ++ "/your/path/DEA-Stata/ado"
which dea
help dea
```

Replace the example path with the location on your computer. Run `which` to check
which copy Stata uses if you already have another DEA installation. This working
copy uses a local-folder installation; it does not yet provide `net install` metadata.

The Examples section in each help file generates its own teaching data, or uses
Stata's bundled `auto` data. Copy its command lines, omitting leading dots, into
the Do-file Editor and run the complete block from `preserve` through `restore`.
No separate dataset or example script is required. Small simulation samples and
replication counts demonstrate syntax, not an empirical research design.

```stata
preserve
clear
input id x y
1 1 1
2 2 3
3 4 4
4 5 3
end
dea x = y, dmu(id) rts(vrs) orientation(input) generate(eff)
list id eff
restore
```

## Choose a command

Use the help link to inspect syntax and examples on GitHub, or type
`help command_name` in Stata for the formatted version. The
[game-model overview](ado/dea_game.sthlp) compares the game-related commands.

| Family | Command | Purpose | Files |
| --- | --- | --- | --- |
| Core DEA | `dea` | Radial, SBM, economic, directional and super-efficiency models | [code](ado/dea.ado) · [help](ado/dea.sthlp) |
| Core DEA | `dea_congestion` | Input congestion | [code](ado/dea_congestion.ado) · [help](ado/dea_congestion.sthlp) |
| Core DEA | `dea_sbmundesirable` | SBM with undesirable outputs | [code](ado/dea_sbmundesirable.ado) · [help](ado/dea_sbmundesirable.sthlp) |
| Time and productivity | `dea_malmquist` | Malmquist productivity | [code](ado/dea_malmquist.ado) · [help](ado/dea_malmquist.sthlp) |
| Time and productivity | `dea_window` | Window analysis | [code](ado/dea_window.ado) · [help](ado/dea_window.sthlp) |
| Cross-efficiency and networks | `dea_crosseff` | Cross-efficiency | [code](ado/dea_crosseff.ado) · [help](ado/dea_crosseff.sthlp) |
| Cross-efficiency and networks | `dea_network` | Two-stage network DEA | [code](ado/dea_network.ado) · [help](ado/dea_network.sthlp) |
| Inference and uncertainty | `dea_bootstrap` | Bootstrap bias correction | [code](ado/dea_bootstrap.ado) · [help](ado/dea_bootstrap.sthlp) |
| Inference and uncertainty | `dea_sw2007` | Two-stage bootstrap inference | [code](ado/dea_sw2007.ado) · [help](ado/dea_sw2007.sthlp) |
| Inference and uncertainty | `dea_fuzzy` | Fuzzy DEA | [code](ado/dea_fuzzy.ado) · [help](ado/dea_fuzzy.sthlp) |
| Exploratory game diagnostics | `dea_crossgame` | Cross-efficiency strategy comparison | [code](ado/dea_crossgame.ado) · [help](ado/dea_crossgame.sthlp) |
| Exploratory game diagnostics | `dea_nash` | Nash-style allocation | [code](ado/dea_nash.ado) · [help](ado/dea_nash.sthlp) |
| Exploratory game diagnostics | `dea_cooperative` | Shapley-style coalitions | [code](ado/dea_cooperative.ado) · [help](ado/dea_cooperative.sthlp) |
| Exploratory game diagnostics | `dea_stackelberg` | Stackelberg-style two-stage diagnostics | [code](ado/dea_stackelberg.ado) · [help](ado/dea_stackelberg.sthlp) |

## Versions and feedback

This is a development working copy. Individual command versions are recorded in
their `.ado` headers; help revisions have a separate dated marker. A help-only
revision does not change the computational command version. A future tagged
release will identify a tested combination of files.

For a problem report, include the command name, Git commit or release tag,
`which command_name` output, Stata version, a small reproducible example, and
expected versus actual results. Use [GitHub Issues](https://github.com/renegade68/DEA-Stata/issues).
Create synthetic example data instead of attaching private research data.

## Attribution and licensing

Preserve the author and reference information in the code and help files.
Cite the article above for the software and each command's references for the
method. The journal citation and article link will be updated when a DOI and
complete publication details become available.

The software is licensed under the [Academic Free License 3.0](LICENSE),
continuing the DEAS license used on SourceForge. See [the software notice](LICENSE_CODE.txt)
for scope and attribution. Article manuscripts and publisher PDFs are not included
and are not covered by this software license.
