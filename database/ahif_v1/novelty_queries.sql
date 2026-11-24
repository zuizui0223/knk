-- AHIF-DB v0.1. Run with: sqlite3 ahif_v1.sqlite < novelty_queries.sql
-- All summaries refer ONLY to the targeted, source-linked pilot.
-- Do not call any values a systematic literature prevalence or pooled effect.

-- 1. Independent source-study cluster counts, not plot counts.
SELECT COUNT(*) AS study_records,
       COUNT(DISTINCT publication_cluster) AS candidate_publication_clusters,
       COUNT(DISTINCT CASE WHEN raw_source_audited='yes' THEN publication_cluster END) AS audited_clusters
FROM studies;

-- 2. Where is identity-versus-rank history already measured?
SELECT s.study_id,s.citation,p.panel_id,p.independent_unit_type,p.response_grain
FROM panels p JOIN studies s ON s.study_id=p.study_id
WHERE p.named_structure_scored='yes' AND p.identity_free_rank_scored='yes'
ORDER BY s.study_id,p.panel_id;

-- 3. The strict causal identification gap in this pilot (expected 0).
SELECT COUNT(*) AS same_panel_mechanism_to_rank_storage_identified
FROM v_bridge_gate
WHERE evidence_status='identified_joint_mechanism_to_storage';

-- 4. Compare intervention pathways without pooling incomparable effects.
SELECT p.intervention_lane,p.putative_mechanism_channel,
       COUNT(DISTINCT s.publication_cluster) AS publication_clusters,
       COUNT(*) AS experimental_panels
FROM panels p JOIN studies s USING(study_id)
WHERE s.raw_source_audited='yes'
GROUP BY p.intervention_lane,p.putative_mechanism_channel
ORDER BY publication_clusters DESC;

-- 5. All experimental effect values are source-linked with original units.
SELECT study_id,publication_cluster,intervention_lane,panel_id,endpoint_class,
       estimand,contrast,estimate,ci95_low,ci95_high,p_value,
       experimental_grain,analysis_status,evidence_url
FROM v_effect_lineage
ORDER BY study_id,panel_id,endpoint_class,contrast;

-- 6. How many source groups have both structure and mechanism measurements?
-- This does not infer that the measurements share an experiment or identify mediation.
WITH scope AS (
 SELECT s.publication_cluster,
        MAX(CASE WHEN p.named_structure_scored='yes' AND p.identity_free_rank_scored='yes'
                 THEN 1 ELSE 0 END) AS has_rank_pair,
        MAX(CASE WHEN p.partner_or_absence_control='yes' OR
                          p.putative_mechanism_channel='conditioned_environment'
                 THEN 1 ELSE 0 END) AS has_process_or_partner_control
 FROM panels p JOIN studies s USING(study_id)
 WHERE s.raw_source_audited='yes'
 GROUP BY s.publication_cluster
)
SELECT SUM(has_rank_pair) AS structure_clusters,
       SUM(has_process_or_partner_control) AS process_clusters,
       SUM(has_rank_pair * has_process_or_partner_control) AS both_in_same_publication_clustering_not_same_panel
FROM scope;

-- 7. Unresolved study/laboratory dependence must remain visible.
SELECT lab_or_overlap_group,COUNT(*) AS candidate_study_rows,
       GROUP_CONCAT(study_id,', ') AS study_ids
FROM studies
GROUP BY lab_or_overlap_group HAVING COUNT(*)>1
ORDER BY candidate_study_rows DESC;

-- 8. Research questions to falsify, not claims already proved.
SELECT hypothesis_id,testable_hypothesis,required_joint_endpoints,
       required_independence_grain,falsifiable_prediction,interpretation_status
FROM hypotheses ORDER BY hypothesis_id;
