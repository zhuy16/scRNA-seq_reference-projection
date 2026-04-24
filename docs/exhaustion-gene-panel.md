# Exhaustion Gene Panel

18 genes curated from published TIL exhaustion literature. Stored in `data/exhaustion_gene_panel.txt` and shared by both the CCA and scVI pipelines (NB03).

```
PDCD1  HAVCR2  LAG3   TIGIT  CTLA4  TOX    NR4A1  ENTPD1  CXCL13
PRDM1  EOMES   TBX21  GZMB   PRF1   IFNG   TNF    CD38    VCAM1
```

## Gene selection rationale

| Gene | Role |
|------|------|
| `PDCD1` (PD-1) | Canonical exhaustion checkpoint receptor |
| `HAVCR2` (TIM-3) | Late exhaustion; co-expressed with PD-1 in severe exhaustion |
| `LAG3` | Inhibitory co-receptor; upregulated in exhaustion |
| `TIGIT` | Exhaustion checkpoint; target of combination immunotherapy |
| `CTLA4` | Activation-induced inhibitory receptor |
| `TOX` | Master transcription factor driving exhaustion programme |
| `NR4A1` | Exhaustion-associated nuclear receptor; represses effector function |
| `ENTPD1` (CD39) | Marks tumour-antigen-experienced exhausted T cells |
| `CXCL13` | Strongly associated with antigen-specific exhausted TIL |
| `PRDM1` (BLIMP1) | Terminal differentiation TF; marks late exhaustion |
| `EOMES` | TF associated with exhausted/dysfunctional CD8 |
| `TBX21` (T-bet) | Effector/early-exhaustion marker; ratio with EOMES used as readout |
| `GZMB` | Cytotoxic effector; reduced in terminally exhausted cells |
| `PRF1` | Perforin; cytotoxic capacity marker |
| `IFNG` | Effector cytokine; attenuated in exhaustion |
| `TNF` | Effector cytokine; co-expressed with IFN-γ in functional T cells |
| `CD38` | Activation / exhaustion marker |
| `VCAM1` | Expressed on stem-like exhausted T cells (TCF1+ progenitors) |

## Source references

- Wherry EJ et al. (2007) *Molecular signature of CD8+ T cell exhaustion during chronic viral infection.* Immunity 27(4):670–684. https://doi.org/10.1016/j.immuni.2007.09.006
- Sade-Feldman M et al. (2018) *Defining T cell states associated with response to checkpoint immunotherapy in melanoma.* Cell 175(4):998–1013. https://doi.org/10.1016/j.cell.2018.10.038
- Oliveira G et al. (2021) *Phenotype, specificity and avidity of antitumour CD8+ T cells in melanoma.* Nature 596(7870):119–125. https://doi.org/10.1038/s41586-021-03704-y
