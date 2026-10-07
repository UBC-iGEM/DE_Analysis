"""Map gene candidates to curated transcription units and promoter annotations.

Annotated promoter sequences are reference context, not cloning boundaries.
This stage deliberately requires operator/fragment review and wet-lab validation.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def read_product(path):
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")
    frame.columns = [re.sub(r"^\d+\)", "", str(c)).strip() for c in frame.columns]
    return frame.map(lambda value: value.strip())


def map_promoters(candidates, tus, promoters, release):
    required_tu = {"id", "tuGenes", "promoterId", "operonName"}
    required_promoter = {"id", "name", "strand", "posTSS", "sequence", "confidenceLevel"}
    if required_tu - set(tus) or required_promoter - set(promoters):
        raise ValueError("RegulonDB TUSet or PromoterSet schema is missing required fields")
    if promoters.id.duplicated().any():
        raise ValueError("Duplicate promoter IDs in reference annotations")
    promoter_lookup = promoters.set_index("id").to_dict("index")
    gene_tus = {}
    for tu in tus.to_dict("records"):
        for gene in tu["tuGenes"].split(";"):
            gene_tus.setdefault(gene.strip().lower(), []).append(tu)
    records = []
    for candidate in candidates.to_dict("records"):
        gene = str(candidate["gene"]).lower()
        matches = gene_tus.get(gene, [])
        if not matches:
            records.append({"gene": gene, "mapping_status": "no_curated_transcription_unit", "regulondb_release": release,
                            "construct_status": "requires_annotation_and_fragment_review"})
        for tu in matches:
            promoter = promoter_lookup.get(tu["promoterId"], {})
            records.append({"gene": gene, "group": candidate.get("group", ""),
                            "evidence_tier": candidate.get("evidence_tier", ""),
                            "tu_id": tu["id"], "operon_name": tu["operonName"], "tu_genes": tu["tuGenes"],
                            "promoter_id": tu["promoterId"], "promoter_name": promoter.get("name", ""),
                            "strand": promoter.get("strand", ""), "tss": promoter.get("posTSS", ""),
                            "sigma_factor": promoter.get("sigmaFactor", ""),
                            "annotated_sequence": promoter.get("sequence", ""),
                            "promoter_confidence": promoter.get("confidenceLevel", ""),
                            "promoter_evidence": promoter.get("promoterEvidence", ""),
                            "tu_confidence": tu.get("confidenceLevel", ""), "pmids": promoter.get("pmids", ""),
                            "regulondb_release": release,
                            "mapping_status": "curated_promoter_mapped" if promoter else "tu_without_promoter_annotation",
                            "construct_status": "requires_operator_and_fragment_review; not_ready_to_order"})
    mapped = pd.DataFrame(records) if records else pd.DataFrame(columns=["gene", "mapping_status", "promoter_id", "regulondb_release", "construct_status"])
    summaries = []
    if not mapped.empty and "promoter_id" in mapped:
        for promoter_id, group in mapped[mapped.promoter_id.fillna("").ne("")].groupby("promoter_id", sort=True):
            first = group.iloc[0].to_dict()
            first["candidate_genes"] = ", ".join(sorted(set(group.gene)))
            first["n_candidate_genes"] = group.gene.nunique()
            first["candidate_groups"] = ", ".join(sorted(set(group.group.fillna(""))))
            first["candidate_evidence_tiers"] = ", ".join(sorted(set(group.evidence_tier.fillna(""))))
            first["tu_ids"] = ", ".join(sorted(set(group.tu_id)))
            first.pop("gene", None)
            first.pop("tu_id", None)
            summaries.append(first)
    review = pd.DataFrame(summaries) if summaries else pd.DataFrame(columns=["promoter_id", "candidate_genes", "n_candidate_genes", "regulondb_release", "construct_status"])
    return mapped, review


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=root / "results/promoter_candidates/annotated_candidates.csv")
    parser.add_argument("--tu", type=Path, default=root / "data/references/regulondb/TUSet.tsv")
    parser.add_argument("--promoters", type=Path, default=root / "data/references/regulondb/PromoterSet.tsv")
    parser.add_argument("--lock", type=Path, default=root / "config/regulondb.lock.json")
    parser.add_argument("--out", type=Path, default=root / "results/promoter_review")
    args = parser.parse_args(argv)
    # Refuse to silently combine promoter annotations from a different release.
    from promoter_discovery.setup_data import sha256_file
    lock = json.loads(args.lock.read_text())
    releases = set()
    for path in (args.tu, args.promoters):
        record = next((r for r in lock["assets"].values() if (root / r["path"]).resolve() == path.resolve()), None)
        if not record or record["sha256"] != sha256_file(path):
            raise ValueError(f"Annotation is not verified by the reference lock: {path}")
        releases.add(record.get("regulondb_release"))
    if len(releases) != 1 or None in releases:
        raise ValueError("TU and promoter annotations must have the same verified RegulonDB release")
    release = releases.pop()
    candidates = pd.read_csv(args.candidates)
    if "regulondb_release" not in candidates or not candidates.regulondb_release.astype(str).eq(release).all():
        raise ValueError("Candidate regulatory evidence and promoter annotations must use the same release; rebuild and rescore the network")
    mapped, review = map_promoters(candidates, read_product(args.tu), read_product(args.promoters), release)
    args.out.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(args.out / "candidate_promoter_mapping.csv", index=False)
    review.to_csv(args.out / "promoter_review.csv", index=False)
    print(f"[promoters] {len(mapped)} gene/TU mappings; {len(review)} distinct annotated promoters; fragment review required")


if __name__ == "__main__":
    main()
