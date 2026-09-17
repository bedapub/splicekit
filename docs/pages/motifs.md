# Motif & RNA-binding analysis

`splicekit motifs` analyzes the sequences around regulated splicing events found by [edgeR](edgeR.md).

```bash
splicekit motifs           # run all: motif logos, DREME, scanRBP
splicekit motifs dreme     # only DREME
```

Motif analysis on donor site patterns (9nt sequences) runs on the top 100 hits of each comparison, producing motif logos and HTML reports under `results/motifs`. In addition to the logos, splicekit runs [DREME](https://meme-suite.org/meme/doc/dreme.html) on regulated sequences vs. control sequences to find enriched short motifs.

## scanRBP: RNA-protein binding enrichment

As part of the same motif analysis, splicekit identifies potential enrichment of RNA-protein binding at regulated sites (donor sites, acceptor sites and other regions), using [scanRBP](https://github.com/grexor/scanRBP).

Once sets of control and regulated sequences are identified, scanRBP computes the log-odds of the binding signal for a chosen protein from its PWM. Bootstrapping the sequence labels estimates the probability that binding at regulated sequences differs from binding at controls (a log-FC of the binding signal).

Configure which protein to scan in `splicekit.config`:

```python
scanRBP = True                    # run the scanRBP step? (True/False)
protein = "K562.TARDBP.0"         # PWM id, see: scanRBP search <term>
protein_label = "tdp43"           # short label used in file names and titles
```

See [Configuration](configuration.md#scanrbp-parameters) for the full parameter list, and the [scanRBP documentation](https://grexor.github.io/scanRBP/) for the standalone tool (`pip install scanRBP`) and its motif database.
