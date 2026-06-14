from src.pathways import (
    gene_token_guard, extract_key_token,
    _normalize, fuzzy_canonical
)

def test_extract_key_token_gene_patterns():
    assert extract_key_token("Signaling by BRCA1 mutants") == "BRCA1"
    assert extract_key_token("Signaling by KRAS in Cancer") == "KRAS"
    assert extract_key_token("Loss of Function of TP53 in Cancer") == "TP53"
    assert extract_key_token("Defective FOXF1 causes ACDMPV") == "FOXF1"
    assert extract_key_token("EGFR mediated signaling") == "EGFR"
    assert extract_key_token("Drug resistance of ALK mutants") == "ALK"
    assert extract_key_token("Signaling by PI3K/AKT") == "PI3K/AKT"

def test_extract_key_token_non_templated():
    assert extract_key_token("Wnt signaling pathway") is None
    assert extract_key_token("JAK-STAT signaling pathway") is None
    assert extract_key_token("Apoptosis") is None
    assert extract_key_token("mRNA Translation") is None

def test_gene_token_guard_rejects_wrong_gene():
    assert gene_token_guard(
        "Signaling by BRCA1 mutants",
        "Signaling by AMER1 mutants") == False
    assert gene_token_guard(
        "Signaling by KRAS in Cancer",
        "Signaling by LTK in cancer") == False
    assert gene_token_guard(
        "Loss of Function of TP53 in Cancer",
        "Loss of Function of SMAD4 in Cancer") == False

def test_gene_token_guard_accepts_same_gene():
    assert gene_token_guard(
        "Signaling by BRCA1 mutants",
        "Signaling by BRCA1 mutants") == True
    assert gene_token_guard(
        "Signaling by EGFR in Cancer",
        "Signaling by EGFR") == True

def test_gene_token_guard_accepts_non_templated():
    assert gene_token_guard(
        "Wnt signaling pathway",
        "Signaling by WNT") == True
    assert gene_token_guard(
        "JAK-STAT signaling",
        "Signaling by JAK-STAT") == True

def test_normalize_rhsa_suffix():
    assert _normalize(
        "Oxidative Stress Induced Senescence (R-HSA-2559580)"
    ) == "oxidative stress induced senescence"

def test_normalize_whitespace():
    assert _normalize("  RAF/MAP kinase cascade  ") == \
        "raf/map kinase cascade"

def test_fuzzy_canonical_exact_match():
    ref = {"Signaling by WNT", "RAF/MAP kinase cascade",
           "Transcriptional Regulation by TP53"}
    synonyms = {}
    is_canon, name, method = fuzzy_canonical(
        "Signaling by WNT", ref, synonyms)
    assert is_canon == True
    assert name == "Signaling by WNT"
    assert method == "exact"

def test_fuzzy_canonical_rhsa_exact():
    ref = {"Signaling by WNT"}
    synonyms = {}
    is_canon, name, method = fuzzy_canonical(
        "Signaling by WNT (R-HSA-195721)", ref, synonyms)
    assert is_canon == True
    assert method == "exact"

def test_fuzzy_canonical_synonym():
    ref = {"Signaling by WNT"}
    synonyms = {"wnt signaling pathway": "Signaling by WNT"}
    is_canon, name, method = fuzzy_canonical(
        "Wnt signaling pathway", ref, synonyms)
    assert is_canon == True
    assert name == "Signaling by WNT"
    assert method == "synonym"

def test_fuzzy_canonical_non_canonical():
    ref = {"Signaling by WNT", "RAF/MAP kinase cascade"}
    synonyms = {}
    is_canon, name, method = fuzzy_canonical(
        "FOXF1-EZH2-DKK3 axis", ref, synonyms, threshold=85)
    assert is_canon == False
    assert method == "non_canonical"

def test_fuzzy_canonical_gene_guard_prevents_false_positive():
    ref = {"Signaling by AMER1 mutants", "Signaling by BRCA1 mutants"}
    synonyms = {}
    is_canon, name, method = fuzzy_canonical(
        "Signaling by BRCA1 mutants", ref, synonyms, threshold=85)
    assert is_canon == True
    assert name == "Signaling by BRCA1 mutants"
