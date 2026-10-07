from promoter_discovery.build_database import build_database
from promoter_discovery.db_api import Database


def test_all_target_types_keep_evidence_and_missing_links(extended_database_root, tmp_path):
    database = Database(build_database(extended_database_root, tmp_path / "project.sqlite"))
    paths = database.reference_promoter_paths("RyhB")
    assert {row["interaction_id"] for row in paths} == {"I3", "I5", "I6"}
    assert all(row["target_kind"] == "gene" for row in paths)
    assert all(row["path_basis"] == "gene_tu_context" for row in paths)
    unmapped = next(row for row in paths if row["interaction_id"] == "I6")
    assert unmapped["gene_id"] == "unmappedgene"
    assert unmapped["promoter_id"] is None and unmapped["promoter_link"] == "missing"
    sigma = database.reference_promoter_paths("sigma70")
    assert {row["gene_id"] for row in sigma} == {"genea", "geneb"}
    assert all(row["target_kind"] == "tu" and row["promoter_link"] == "context" for row in sigma)
    crp = database.reference_promoter_paths("CRP")
    assert all(row["promoter_link"] == "direct" for row in crp)
    coverage = database.reference_coverage()
    assert coverage["regulatory_interactions"] == 7
    assert coverage["interactions_by_target"] == {"gene": 3, "promoter": 2, "tu": 2}
    assert coverage["interactions_without_promoter"] == 2
    database.close()


def test_candidate_graph_keeps_direct_targets_and_excludes_sibling_gene_effects(extended_database_root, tmp_path):
    database = Database(build_database(extended_database_root, tmp_path / "project.sqlite"))
    paths = database.get_candidate("GeneA")["regulatory_paths"]
    assert {row["interaction_id"] for row in paths} == {"I1", "I3", "I4"}
    graph = database.candidate_network("genea")
    direct = {edge["id"]: edge for edge in graph["edges"] if edge["kind"] == "regulation"}
    assert direct["I1"]["target"] == "promoter:P1"
    assert direct["I3"]["target"] == "gene:genea"
    assert direct["I4"]["target"] == "tu:T1"
    assert direct["I3"]["pmids"] == "3"
    assert {node["kind"] for node in graph["nodes"]} == {"actor", "promoter", "tu", "gene"}
    assert any(edge["kind"] == "context" for edge in graph["edges"])
    assert database.candidate_network("unknown") is None
    database.close()


def test_queries_default_to_one_completed_run(extended_database_root, tmp_path):
    database = Database(build_database(extended_database_root, tmp_path / "project.sqlite"))
    original = database.runs()[0]["run_id"]
    with database.connection:
        database.connection.execute("INSERT INTO analysis_runs(run_id,created_at,status) VALUES ('new','2099-01-01','complete')")
        database.connection.execute("INSERT INTO candidates(run_id,candidate_id,gene_id,antibiotic_class) VALUES ('new','genea','genea','aminoglycoside')")
        database.connection.execute("INSERT INTO regulatory_interactions(interaction_id,run_id,actor_id,interaction_type,target_kind,gene_id,effect,source_release,source_json) VALUES ('new-i','new','R2','sRNA-gene','gene','genea','repressor','14.5.0','{}')")
    assert [row["run_id"] for row in database.search_candidates()] == ["new"]
    assert database.get_candidate("genea")["run_id"] == "new"
    assert {row["run_id"] for row in database.reference_promoter_paths("RyhB")} == {"new"}
    assert {row["run_id"] for row in database.reference_interactions("RyhB")} == {"new"}
    assert {row["run_id"] for row in database.candidate_paths("genea")} == {"new"}
    assert database.panel() == []
    assert database.sample_qc() == []
    assert database.get_candidate("genea", original)["antibiotic_class"] == "beta_lactam"
    assert {row["run_id"] for row in database.reference_promoter_paths("RyhB", original)} == {original}
    database.close()
