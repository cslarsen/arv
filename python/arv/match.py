"""
Contains matching functions.

Part of arv
Copyright 2014, 2016, 2017 Christian Stigen Larsen
Distributed under the GPL v3 or later. See COPYING.
"""

def assert_european(genome):
    """If ethnicity is set, make sure it's European."""
    if genome.ethnicity not in [None, "european"]:
        raise ValueError("Only applicable to Europeans")

def unphased_match(snp, phenotypes):
    """Match SNP with unphased genotypes and return phenotype.

    Disregards phasing when comparing genotypes, meaning that an input value of
    "AG" will be matched against both "AG" and "GA".

    Arguments:
        genotype: Genotype (str) or SNP (arv.SNP) to match.
        phenotypes: Dict mapping (unphased) genotype to phenotype.

    Example:
        unphased_match(genome.rs4988235, {
            "AA": "Likely lactose tolerant",
            "AG": "Likely lactose tolerant",
            "GG": "Likely lactose intolerant",
            None: "Unknown genotype"})

        The above example could return "Likely lactose tolerant", for example,
        or "Unknown genotype" if there was no match. Note that the key "AG"
        will match both "AG" and "GA" in the snp.

    Returns:
        Matching phenotype. If the `phenotypes` dict has a `None` key, it will
        be returned in case there is no match.
    """
    if isinstance(snp, str):
        genotype = snp
    else:
        genotype = str(snp)

    # Look for "IJ"
    if genotype in phenotypes:
        return phenotypes[genotype]

    # Look for "JI"
    genotype = "".join(reversed(str(snp)))
    if genotype in phenotypes:
        return phenotypes[genotype]

    # Use default value?
    if None in phenotypes:
        return phenotypes[None]
    else:
        raise KeyError(str(snp))
