from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Dict, Set, Tuple


def pathway_to_binary_vector(
    ranked_genes: pd.DataFrame,
    pathway_genes: Set[str],
) -> np.ndarray:
    """
    ranked_genes: DataFrame with column 'gene' in ranked order (0..N-1)
    pathway_genes: set of gene IDs in this pathway
    """
    genes = ranked_genes["gene"].astype(str).values
    v = np.isin(genes, list(pathway_genes)).astype(int)
    return v
