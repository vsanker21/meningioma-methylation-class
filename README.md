# Meningioma methylation class

Analysis code for three measurements in meningioma: chromosomal loss on the methylation array, cell composition from matched RNA, and resemblance of the methylation profile to sorted normal cells. A further comparison uses the NF2 sequencing counts in the Nassiri series.

The public files read by the code are named in the data-availability statement of the manuscript.

- `purity_conumee.R` measures chromosomal loss from the array signal.
- `build_scrna_signature.R` and `rna_deconvolution.py` measure cell composition from matched RNA.
- `write_eligible_probes.py`, `methylation_reference.R` and `methylation_deconvolution.R` compare tumour methylation with sorted normal cells.
- `purity_scale_comparison.py` and `purity_scale_mechanism.py` compare the array loss with the RNA estimate.
- `scrna_chromosome_dosage.R` compares chromosome 22 expression in neoplastic and macrophage nuclei.
- `nf2_vaf_vs_intensity.py` compares NF2 sequencing counts with the array loss in the Nassiri series.
