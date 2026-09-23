# Oral 16S release-shift study

This repository is independent of the submitted CREO-Research project. It contains an exploratory paired oral 16S reference-sequence comparison and a [40-reference author manuscript](manuscript.md).

**Main finding:** On the same 2,236 v16.03 sequences from historically represented HMTs, a full-length 16S nearest-reference classifier assigned the correct HMT on **90.03%** of queries, versus **80.23%** using in-silico V3–V4 fragments. The gain was **9.79 percentage points** (95% HMT-cluster-bootstrap interval **5.84–14.67 points**). The paired novelty AUROC difference had an interval spanning zero. This analysis was devised after an initial release inspection and requires external specimen or isolate confirmation.

## Reproduce

Use Python 3.11+ and install the version-pinned dependencies. The public FASTA data are fetched separately and excluded from Git:

~~~bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
python homd_data.py
python fetch_external_data.py
python homd_deep_benchmark.py
python paired_region_experiment.py
python homd_equal_budget.py
python homd_taxon_balance.py
python make_figures.py
~~~

The benchmark consumes public [HOMD reference files](https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/) and [16SGOSeq genome-derived sequences](https://zenodo.org/records/15209015). Check downloaded HOMD byte SHA-256 values against source_hashes.json. The 16SGOSeq downloader verifies its own SHA-256. Access to the same historical releases is needed for exact reproduction. The full-length reference comparison contains no patient data.

| File | Contents |
|:--|:--|
| manuscript.md, references.json | Exploratory paper and 40 verified DOI records |
| homd_data.py, fetch_external_data.py | Source downloads and taxonomy parser |
| paired_region_experiment.py, paired_region_experiment_results.json | Matched query study, paired counts and taxon bootstrap |
| homd_deep_benchmark.py, homd_deep_results.json | Temporal and same-release context analyses |
| homd_equal_budget.py, homd_equal_budget_results.json | Equal-query-size baseline splits |
| homd_taxon_balance.py, homd_taxon_balance_results.json | One-sequence-per-HMT sensitivity |
| source_hashes.json, requirements.txt | Exact inputs and environment |
| make_figures.py, figure_*.png | Plot source and regenerated figures |

**Limits:** HMTs are database taxa and do not universally correspond to named species. Both region representations come from curated reference sequences rather than laboratory reads. Fragment extraction and duplicate removal change the eligible query population. This repository offers a research draft, not a peer-reviewed or clinical result.
