# Command-line reference

Every `splicekit` invocation prints its version first, then runs the requested command. Run `splicekit -help` for the built-in summary; this page is the fuller reference.

## Processing

| Command | Description |
|---|---|
| `splicekit process` | Run all available analyses in order: setup, annotation, features, edgeR, juan, judge, motifs, promisc, clusterlogfc, june, report, jbrowse2. |

## Initialization

| Command | Description |
|---|---|
| `splicekit setup` | Initialize the project folder structure. |
| `splicekit annotation` | Load `samples.tab` and build comparisons. See [Sample annotation](samples.md). |
| `splicekit features` | Create junction/exon/anchor/gene count tables from BAM files. See [Features & count tables](features.md). |

## Splicing analyses

| Command | Description |
|---|---|
| `splicekit edgeR [junctions\|exons\|anchors\|genes]` | Run edgeR differential usage analysis. All feature types if no sub-command is given. See [Differential splicing (edgeR)](edgeR.md). |
| `splicekit dexseq [junctions\|exons\|anchors\|genes]` | Run DEXSeq as an alternative to edgeR. See [Additional analyses](additional-analyses.md#dexseq). |
| `splicekit juan` | Merge donor/acceptor anchor edgeR results into the junction results. |
| `splicekit judge` | Generate juDGE plots (junction logFC vs. gene logFC). See [juDGE plots](judge.md). |
| `splicekit promisc` | Promiscuity analysis: genes and their junction changes across conditions. |
| `splicekit clusterlogfc` | Cluster comparisons by logFC of significant (FDR < 0.05) changes, at junction/exon/gene level. |
| `splicekit june` | JUNE junction-event analysis (e.g. skipped/mutually exclusive exons). |
| `splicekit rmats` | Process the project's BAM files with rMATS. |
| `splicekit cassettes` | Cassette exon analysis. |
| `splicekit patterns` | Sequence pattern analysis. |

## Motif & scanRBP

| Command | Description |
|---|---|
| `splicekit motifs [dreme\|scanrbp]` | Run motif logos, DREME and scanRBP analysis. All three if no sub-command is given. See [Motif & RNA-binding analysis](motifs.md). |

## JBrowse2 & reporting

| Command | Description |
|---|---|
| `splicekit jbrowse2 [process\|start]` | Process (build) JBrowse2 files and/or start the JBrowse2 web server. Both steps if no sub-command is given. |
| `splicekit jbrowse` | Alias for `splicekit jbrowse2`. |
| `splicekit web` | Start the local web server serving both the HTML report and JBrowse2. See [Exploring results](jbrowse2.md). |
| `splicekit report` | Generate the HTML report under `report/`. |

## Other

| Command | Description |
|---|---|
| `splicekit version` | Print the installed splicekit version. |

## Global options

| Option | Description |
|---|---|
| `-help` | Print the built-in usage summary (or a sub-command's usage, e.g. `splicekit edgeR -help`). |
| `-version` | Print the installed version and exit. |
| `-verbose` | Print more detailed progress output. |
| `-force` | Force recomputation of steps that would otherwise be skipped if their output already exists (used by `process` and `jbrowse2`). |
