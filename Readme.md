```markdown
# 🚀 Usage

Run the script from the command line:

```bash
Rscript phyloPCoA.R [options]
```

## Example (Simulation Mode)

```bash
Rscript phyloPCoA.R \
  --bd 0.2,0.1,0.83,100 \
  -T 30 -B 15 \
  --outdir haha \
  --force \
  --exponent 1 \
  --pagel_lam_mode auto \
  --pagel_lam_sim 0.7 \
  --sim \
```

To draw simulated lambda values from a beta distribution:

```bash
Rscript phyloPCoA.R \
  --bd 0.2,0.1,0.83,100 \
  -T 30 -B 15 \
  --outdir haha \
  --force \
  --exponent 1 \
  --pagel_lam_mode hierarchical \
  --pagel_lam_sim beta:2,5 \
  --sim
```

## 🔍 Parameter Explanation

| Option | Description |
|---------|-------------|
| `--bd 0.2,0.1,0.83,100` | Birth–Death parameters: lambda=0.2, mu=0.1, sampling_ratio=0.83, root_age=100 |
| `-T 30` | Number of host species (tnum) |
| `-B 15` | Number of bacterial taxa/features (bnum) |
| `--exponent 1` | Controls feature correlation strength in the covariance matrix (Rho); higher → stronger correlation |
| `--outdir haha` | Output directory for all results |
| `--force` | Overwrite existing output directory |
| `--sim` | Enable simulation mode |
| `--species_exclude` | Drop listed species before analysis |
| `--sim_discrete_trait` | Simulate a binary host trait in simulation mode |
| `--pagel_lam_mode auto` | Pagel's lambda mode: `auto`, `global`, `per_feature`, `hierarchical`, or `none`; `auto` compares modes by AIC |
| `--pagel_lam_sim 0.7` | Simulate ground-truth Pagel lambda values for all features |
| `--pagel_lam_sim beta:2,5` | Simulate ground-truth per-feature Pagel lambda values drawn from `Beta(2,5)` |

## 🧬 Main Parameters

| Parameter | Description |
|-----------|-------------|
| `tnum` | Number of host species (tree tips) |
| `bnum` | Number of bacterial taxa/features |
| `exponent` | Controls how correlated the bacterial features are (higher = more correlated) |
| `BD` | Birth-Death parameters: speciation rate (λ), extinction rate (μ), sampling ratio (ρ), and root age |
| `pagel_lam_mode` | Controls how Pagel's lambda is estimated and applied |
| `pagel_lam_sim` | Controls simulated ground-truth Pagel lambda values when `--sim` is enabled |

## 🧩 Input and Output Files (In Simulation Mode)

All outputs are written to `--outdir`.

| File | Description |
|------|-------------|
| `sim.tree` | Simulated host phylogeny |
| `Sigma.tbl` | Covariance matrix of simulated microbial relative abundance |
| `Rho.tbl` | Correlation matrix of microbial relative abundance |
| `prop.tbl` | Relative abundance of microbes |
| `log_prop.tbl` | Log-transformed relative abundance data |
| `log_prop_geomean.tbl` | Log-proportion table used for phylogenetic transforms |
| `P.tbl` | Phylogenetically transformed abundance matrix |
| `prop2.tbl` | Back-transformed proportions from `P` |
| `pagel_lam_estimates.tbl` | Pagel's lambda estimates by feature; hierarchical mode reports posterior mean, MAP lambda, and MLE beta alpha/beta hyperparameters |
| `pagel_lam_model_selection.tbl` | AIC comparison table written when `--pagel_lam_mode auto` is used |
| `pcoa.pdf` | Two pages of PCoA visualizations, each with four PCoA plots (including Bray-Curtis on `prop2.tbl`) and one host tree |
| `adonis/` | PERMANOVA results for PCoA axes |
| `compare/` | Statistics and comparison results between PCoA and true covariance |

Example directory layout:

```
haha/
├── prop.tbl
├── log_prop.tbl
├── P.tbl
├── Sigma.tbl
├── Rho.tbl
├── sim.tree
├── pcoa.pdf
└── compare/
    ├── explained.tbl
    ├── determined_by_trait.tbl
    ├── determined_by_phylo.tbl
    └── compared_to_R_matrix.tbl
```

With `--pagel_lam_mode auto`, the output also includes `pagel_lam_model_selection.tbl`.

## 📊 Analysis Components

| Metric | Description |
|--------|-------------|
| **LDA Accuracy** | Classification accuracy for group assignments |
| **Fisher's Discriminant Ratio (FDR)** | Ratio of between-group variance to within-group variance |
| **Davies–Bouldin Index (DBI)** | Cluster separation index (lower = more separated) |
| **Correlation vs Rho** | Measures alignment between PCoA eigenstructure and true covariance |

## 🔬 Workflow Summary

1. Simulate host phylogeny using the BD model.
2. Generate microbial abundance matrix with covariance structure defined by `exponent`.
3. Transform abundance data via CLR and phylogenetically corrected CLR (pCLR).
4. Perform PCoA for multiple distance metrics.
5. Evaluate clustering and correspondence with the true covariance (Rho).
6. Save all outputs to the chosen directory.

## 🧠 Example Summary

Example command:

```bash
Rscript phyloPCoA.R --bd 0.2,0.1,0.83,100 -T 30 -B 15 --outdir haha --force --exponent 1 --sim
```

This creates a 30‑species host tree (`tnum=30`) with 15 bacterial taxa (`bnum=15`), simulates correlated abundance data (`exponent=1`), performs PCoA, and stores all outputs in `haha/`.

## 🖋️ Authors

- **Youhua Chen** — A first draft code.
- **Sishuo Wang** — Code development, simulation implementation, conceptualization, and major updates (2024–2026)

## License
© 2024-2026 Sishuo Wang. All rights reserved.
This software is provided for reference and demonstration purposes only.  
Reproduction or distribution of this code without written permission is prohibited.
