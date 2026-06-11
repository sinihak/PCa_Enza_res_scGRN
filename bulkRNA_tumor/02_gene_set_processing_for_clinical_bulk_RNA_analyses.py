
import utils as ut
import pandas as pd
import json
from pathlib import Path
datadir = Path.cwd()


gs= f"{datadir}/data/genesets_to_plot.xlsx"
gene_sets = ut.load_gene_sets_excel(gs)


keys_to_keep = ["ARS_stress_linked_sensitive","SOX4_linked_resistance","SOX4_regulon","Liu_PTEN_loss_UP","Liu_PTEN_loss_DOWN","hallmark_pi3k_akt_mtor_signaling","AR_activity","Basal", "Duct luminal","Luminal",
"NE_identity","Beltran_NEPC_UP","Beltran_NEPC_DOWN","hallmark_wnt_beta_catenin_signaling","AR_variant_ARV"]

gene_sets = {k: gene_sets[k] for k in keys_to_keep if k in gene_sets}

# ADD SOX4 to its regulon:
if "SOX4" not in gene_sets["SOX4_regulon"]:
    gene_sets["SOX4_regulon"].append("SOX4")

# ignore the canonical AR response genes in the sensitive sets
ar_genes = set(extra_gene_sets["hallmark_androgen_response"])
for key, genes in gene_sets.items():
    if "sensitive" in key:
        # remove androgen response genes
        gene_sets[key] = list(set(genes) - ar_genes)

# ignore the ribosomal genes
gene_sets["ARS_stress_linked_sensitive"] = (pd.Series(gene_sets["ARS_stress_linked_sensitive"])
      .loc[~pd.Series(gene_sets["ARS_stress_linked_sensitive"]).str.contains("RPL|RPS", na=False)]
      .tolist()
)

# Now remove overlaps for correlation analyses
gs_subset = ["Basal", "Duct luminal","Luminal","NE_identity","AR_activity","Beltran_NEPC_UP","Beltran_NEPC_DOWN","hallmark_wnt_beta_catenin_signaling","AR_variant_ARV","hallmark_epithelial_mesenchymal_transition"]
for gs in gs_subset:
    # remove SOX4 targets from the gene sets for correlation analysis
    gene_sets[gs] = list(set(gene_sets[gs]) - set(gene_sets["SOX4_regulon"]))
    # also, remove the gene set-associated genes from the larger "SOX4_associated_resistance"
    gene_sets["SOX4_associated_resistance"] = list(set(gene_sets["SOX4_associated_resistance"]) - set(gene_sets[gs]))
    gene_sets["SOX4_associated_resistance"] = list(set(gene_sets["SOX4_associated_resistance"]) - set(gene_sets["SOX4_regulon"])) # remove overlaps between the core regulon and trans-regulatory network
    gene_sets["ARS_stress_linked_sensitive"] = list(set(gene_sets["ARS_stress_linked_sensitive"]) - set(gene_sets[gs]))



with open(f"{datadir}/data/processed_gene_sets.json", "w") as fp:
    json.dump(gene_sets, fp)