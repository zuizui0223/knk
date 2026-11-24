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

-- 9. Independent original 2026 field rank and spontaneous-species sensitivity
-- Same physical 40 plots repeated in four survey years; values not independent
-- study effects or independently identified competition coefficients.
SELECT e.contrast,e.endpoint_class,e.estimate,e.p_value,e.experimental_grain,
       s.publication_cluster,s.lab_or_overlap_group,e.evidence_url
FROM effects e
JOIN panels p USING(panel_id)
JOIN studies s USING(study_id)
WHERE p.panel_id='P_VANSTEIJNFIELD2026'
  AND e.endpoint_class IN ('structure','structure_trajectory',
      'structure_nontarget_roster','structure_sown_roster')
ORDER BY e.endpoint_class,e.contrast;

-- 10. Causally incomplete full community panels versus two-species geometry
SELECT p.panel_id,p.study_id,p.response_grain,
       p.same_panel_causal_interaction_to_structure_identified,
       s.raw_source_audited
FROM panels p JOIN studies s USING(study_id)
WHERE p.named_structure_scored='yes'
  AND p.identity_free_rank_scored='yes'
ORDER BY p.study_id;

-- 11. Methodological novelty threats never count as independently
-- measured effects (e.g. sPriority meta-analysis and LOVE 2024).
SELECT record_id,title,doi_or_unverified,source_verification_level,
       AHIF_novelty_overlap,not_independent_effect,source_url
FROM methodological_prior_art
ORDER BY record_id;

-- 12. Paper-level novelty audits, preserving original source-data status.
SELECT p.study_id,p.screening_scope,p.source_raw_audited,
       p.existing_novelty_claim_defeated,p.remaining_testable_gap
FROM prior_art_verified p
ORDER BY p.study_id;

-- 13. Field routing tests; different estimands must not be pooled.
SELECT site_scope,period,endpoint,response_channel,estimate,estimand,p_unadjusted,p_holm,family,publication_cluster FROM v_routing_tests ORDER BY family,site_scope,endpoint;
-- 14. Descriptive cage perturbations never count as treatment-replicated effects.
SELECT qc_lane,focal_species,control_cage,dropout_cage,n_control_bees,n_dropout_bees,strain_named_tv,strain_rank_tv,species_named_tv,species_rank_tv,median_nonfocal_copy_ratio,evidence_class FROM cage_descriptive_contrasts ORDER BY qc_lane,focal_species,dropout_cage;
