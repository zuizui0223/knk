PRAGMA foreign_keys=ON;
CREATE TABLE opportunity_sources (
 source TEXT PRIMARY KEY, publication_doi TEXT NOT NULL UNIQUE,
 source_locator TEXT NOT NULL, design TEXT NOT NULL,
 physical_dependency TEXT NOT NULL, observation_unit TEXT NOT NULL,
 inference_limit TEXT NOT NULL, meta_eligible TEXT NOT NULL CHECK(meta_eligible='no')
);
CREATE TABLE opportunity_lanes (
 lane TEXT PRIMARY KEY, source TEXT NOT NULL REFERENCES opportunity_sources(source),
 resolution TEXT NOT NULL, threshold REAL NOT NULL CHECK(threshold>=0),
 proxy TEXT NOT NULL, measurement TEXT NOT NULL, risk_definition TEXT NOT NULL,
 role TEXT NOT NULL CHECK(role IN ('primary','sensitivity')),
 UNIQUE(lane,source)
);
CREATE TABLE opportunity_taxa (
 source TEXT NOT NULL REFERENCES opportunity_sources(source), taxon TEXT NOT NULL,
 label TEXT NOT NULL, resolution TEXT NOT NULL,
 PRIMARY KEY(source,taxon)
);
CREATE TABLE opportunity_events (
 lane TEXT NOT NULL, source TEXT NOT NULL, event_id TEXT NOT NULL,
 context TEXT NOT NULL, dependency_group TEXT NOT NULL,
 physical_target TEXT NOT NULL, focal_profile TEXT NOT NULL, partner_profile TEXT NOT NULL,
 preceding_label TEXT NOT NULL, target_label TEXT NOT NULL,
 exposure_status TEXT NOT NULL CHECK(exposure_status IN ('no_transfer_placebo','assigned_transfer','assigned_coalescence')),
 orientation TEXT NOT NULL CHECK(orientation IN ('actual_recipient','analytical')),
 nominal_partner_fraction REAL CHECK(nominal_partner_fraction>0 AND nominal_partner_fraction<1),
 event_weight REAL NOT NULL CHECK(event_weight>0),
 PRIMARY KEY(lane,event_id), UNIQUE(lane,event_id,source),
 FOREIGN KEY(lane,source) REFERENCES opportunity_lanes(lane,source)
);
CREATE TABLE opportunity_risk (
 lane TEXT NOT NULL, source TEXT NOT NULL, event_id TEXT NOT NULL, taxon TEXT NOT NULL,
 pool_present INTEGER NOT NULL CHECK(pool_present>=0),
 pool_n INTEGER NOT NULL CHECK(pool_n>0 AND pool_present<=pool_n),
 partner_detected INTEGER NOT NULL CHECK(partner_detected IN (0,1)),
 partner_relative_abundance REAL NOT NULL CHECK(partner_relative_abundance>=0 AND partner_relative_abundance<=1.000000000001),
 target_detected INTEGER NOT NULL CHECK(target_detected IN (0,1)),
 unmixed_control_detected INTEGER CHECK(unmixed_control_detected IN (0,1)),
 PRIMARY KEY(lane,event_id,taxon),
 FOREIGN KEY(lane,event_id,source) REFERENCES opportunity_events(lane,event_id,source),
 FOREIGN KEY(source,taxon) REFERENCES opportunity_taxa(source,taxon)
);
CREATE TABLE opportunity_assays (
 assay_id TEXT PRIMARY KEY, source TEXT NOT NULL REFERENCES opportunity_sources(source),
 medium TEXT NOT NULL, parent_i TEXT NOT NULL, parent_j TEXT NOT NULL,
 pair_id TEXT NOT NULL, physical_sample TEXT NOT NULL, well TEXT NOT NULL,
 monoculture_i_sample TEXT NOT NULL, monoculture_j_sample TEXT NOT NULL,
 n_reads INTEGER NOT NULL CHECK(n_reads>=1000),
 unique_marker_i INTEGER NOT NULL CHECK(unique_marker_i>=0),
 unique_marker_j INTEGER NOT NULL CHECK(unique_marker_j>=0),
 signature_mass_i REAL NOT NULL CHECK(signature_mass_i BETWEEN 0 AND 1.000000000001),
 signature_mass_j REAL NOT NULL CHECK(signature_mass_j BETWEEN 0 AND 1.000000000001),
 unresolved_mass REAL NOT NULL CHECK(unresolved_mass BETWEEN -0.000000000001 AND 1.000000000001),
 signature_share_i REAL CHECK(signature_share_i BETWEEN 0 AND 1),
 monoculture_reference_tv REAL NOT NULL CHECK(monoculture_reference_tv BETWEEN 0 AND 1.000000000001),
 evidence_class TEXT NOT NULL CHECK(evidence_class='controlled_coculture_relative_read_endpoint'),
 UNIQUE(source,physical_sample)
);
CREATE TABLE opportunity_assay_links (
 assay_id TEXT PRIMARY KEY REFERENCES opportunity_assays(assay_id),
 source TEXT NOT NULL REFERENCES opportunity_sources(source),
 community_mixture_sample TEXT NOT NULL, community_control_i TEXT NOT NULL,community_control_j TEXT NOT NULL,
 community_reference_tv REAL NOT NULL CHECK(community_reference_tv BETWEEN 0 AND 1.000000000001),
 link_basis TEXT NOT NULL CHECK(link_basis='same_medium_origin_label_and_well_not_taxon_identity'),
 UNIQUE(source,community_mixture_sample)
);
CREATE TABLE opportunity_signed_interaction_effects (
 effect_id TEXT PRIMARY KEY, source TEXT NOT NULL REFERENCES opportunity_sources(source),
 evidence TEXT NOT NULL
);
CREATE VIEW opportunity_detection_summary AS
 SELECT r.lane,e.context,e.exposure_status,r.partner_detected,
 COUNT(*) AS dependent_risk_rows,COUNT(DISTINCT e.physical_target) AS physical_targets,
 AVG(r.target_detected) AS detection_fraction,
 AVG(1.0*r.pool_present/r.pool_n) AS mean_pool_frequency,
 AVG(r.partner_relative_abundance) AS mean_partner_relative_abundance
 FROM opportunity_risk r JOIN opportunity_events e USING(lane,source,event_id)
 GROUP BY r.lane,e.context,e.exposure_status,r.partner_detected;
CREATE VIEW opportunity_pairwise_assay_bridge AS
 SELECT a.source,a.medium,a.pair_id,COUNT(*) AS paired_wells,
 AVG(a.monoculture_reference_tv) AS isolate_reference_tv,
 AVG(l.community_reference_tv) AS community_reference_tv,
 AVG(a.unresolved_mass) AS unresolved_signature_mass
 FROM opportunity_assays a JOIN opportunity_assay_links l USING(assay_id,source)
 GROUP BY a.source,a.medium,a.pair_id;
CREATE VIEW opportunity_event_brier AS
 SELECT e.lane,e.source,e.event_id,e.context,e.dependency_group,e.physical_target,
 e.exposure_status,e.event_weight,COUNT(*) AS n_risk,
 AVG((1.0*r.pool_present/r.pool_n-r.target_detected)*(1.0*r.pool_present/r.pool_n-r.target_detected)) AS pool_loss,
 AVG((r.partner_detected-r.target_detected)*(r.partner_detected-r.target_detected)) AS donor_loss,
 AVG(((1.0*r.pool_present/r.pool_n+r.partner_detected)/2-r.target_detected)*((1.0*r.pool_present/r.pool_n+r.partner_detected)/2-r.target_detected)) AS hybrid_loss,
 AVG(((1.0*r.pool_present/r.pool_n+e.nominal_partner_fraction*r.partner_detected)/2-r.target_detected)*((1.0*r.pool_present/r.pool_n+e.nominal_partner_fraction*r.partner_detected)/2-r.target_detected)) AS dose_loss,
 AVG(((1.0*r.pool_present/r.pool_n+r.partner_detected)/2-r.unmixed_control_detected)*((1.0*r.pool_present/r.pool_n+r.partner_detected)/2-r.unmixed_control_detected)) AS control_hybrid_loss,
 AVG((1.0*r.pool_present/r.pool_n-r.unmixed_control_detected)*(1.0*r.pool_present/r.pool_n-r.unmixed_control_detected)) AS control_pool_loss
 FROM opportunity_events e JOIN opportunity_risk r USING(lane,source,event_id)
 GROUP BY e.lane,e.source,e.event_id;
