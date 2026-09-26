"""The one prompt every condition uses (plain models in run_eval.py and the agents in run_agent.py)."""

# PBMC3k keeps its original wording so old and new runs stay comparable.
# The Hao prompt names the sample type, so answers like "neutrophil" count as ignoring context.
SAMPLE = {"pbmc": "human blood", "hao": "human peripheral blood mononuclear cells (PBMCs)"}


def make_prompt(dataset, genes):
    return (f"These are the top marker genes of a cluster from {SAMPLE[dataset]}. What cell type is it? "
            'Reply as JSON: {"cell_type": ..., "confidence": 0-1}\n\nGenes: ' + ", ".join(genes))
