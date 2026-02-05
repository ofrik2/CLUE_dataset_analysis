import os
import tempfile

from clue_pathway_enrichment.io.load_signature import load_signature
from clue_pathway_enrichment.io.load_pathways import load_pathways
from clue_pathway_enrichment.preprocessing.gene_id_normalize import normalize_gene_ids
from clue_pathway_enrichment.preprocessing.make_ranked_list import make_ranked_list
from clue_pathway_enrichment.pipeline.run_pipeline import run_pipeline


def test_end_to_end_alpha():
    # create temporary signature file and pathways file
    with tempfile.TemporaryDirectory() as td:
        sig_path = os.path.join(td, "sig.txt")
        pw_path = os.path.join(td, "pws.tsv")
        with open(sig_path, "w") as fh:
            fh.write("geneA\ngeneB\ngeneC\n")
        with open(pw_path, "w") as fh:
            fh.write("path1\tgeneA\npath1\tgeneX\npath2\tgeneB\n")

        sig = load_signature(sig_path)
        assert sig == ["geneA", "geneB", "geneC"]
        sig_norm = normalize_gene_ids(sig)
        assert sig_norm == ["GENEA", "GENEB", "GENEC"]

        pws = load_pathways(pw_path)
        assert "path1" in pws and "path2" in pws

        # create a ranked list (gene, score)
        ranked = make_ranked_list([(g, float(i)) for i, g in enumerate(sig_norm[::-1])])
        assert isinstance(ranked, list)

        scores = run_pipeline(ranked, pws, method="alpha")
        assert isinstance(scores, dict)

