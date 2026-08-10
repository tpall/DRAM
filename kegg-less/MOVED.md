# Moved to `digestome`

This directory held the KEGG-less anaerobic-digestion marker panel and scorer. It
now lives in its own repository, named **`digestome`**.

**Nothing here depended on DRAM.** The Python was standard library only, the sole
external binaries were `hmmsearch`/`hmmfetch`/`hmmpress`, and the two mentions of
DRAM in the source were docstring text rather than calls. Keeping it inside a DRAM
fork implied a dependency a reader would go looking for and not find, and left
original work sitting under DRAM's GPL-3.

DRAM established the pattern the tool follows, distilling per-gene annotations into
pathway-level statements, and remains the more general tool. `digestome` is an
independent implementation for one domain.

## What moved

- `AD_methanogenesis_panel.tsv`: 86 markers over 20 modules, hydrolysis through
  all three methanogenesis branches
- `build_ad_hmm_db.sh`: builds the pressed HMM DB from NCBIfam and Pfam
- `panel_scored.py`, `aggregate_community.py`, `audit_symbol_matches.py`
- `secretion_modules.tsv`: extracellular-targeting modules for the export test
- `hpc/`: cluster entry points, `tests/`: the assertion suite
- `MANUSCRIPT.qmd`, `references.bib`

## History

All 48 commits that touched this directory were carried across with
`git subtree split`, so the reasoning behind the curated accessions, the
gene-symbol collision fixes and the pyrrolysine finding is still readable there
rather than only here.

## If you are looking for the code

Use the `digestome` repository. This directory is a tombstone and will not be
updated; anything found in the git history here is superseded.
