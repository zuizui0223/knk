# AHIF opportunity-to-outcome registry: query guide

The new SQLite tables are a derived evidence layer on top of the frozen mechanism v0.5 registry. They contain source-specific candidate taxon × exposure records and controlled co-culture measurements; they are not a global species interaction network. Build via `python -m benchmarks.build_ahif_opportunity_registry_v06 --out <output> --scratch <temporary source directory>`. Original source bytes remain temporary. Source-local marker keys are hashed to avoid accidental cross-study identity merges and raw sequence redistribution.

The analytical grain is `opportunity_risk(lane,event_id,taxon)`. Source and taxon ID jointly identify a local taxon. `opportunity_events.physical_target` identifies the actual target sample, while `orientation='analytical'` indicates symmetric focal/partner views of that same sample; those two records are not independent. Lanes are thresholds/proxy sensitivity variants, not independent studies. `opportunity_events.event_weight` reproduces the previous scoring weights within each comparison context. Goldman's nominal partner fraction is a volume fraction, not an observed cell count. The source pool deliberately excludes the actual partner and, for symmetric mixes, the focal parent. Villarreal uses all undetected author bins; Diaz/Goldman restrict to the regional pre-profile roster.

## Example queries
Coverage of source studies, lanes and physical samples:
```sql
SELECT l.source,l.lane,l.role,COUNT(DISTINCT e.physical_target) AS targets,
       COUNT(DISTINCT e.event_id) AS oriented_events,COUNT(*) AS dependent_risk_rows
FROM opportunity_risk r
JOIN opportunity_events e USING (lane,source,event_id)
JOIN opportunity_lanes l USING (lane,source)
GROUP BY l.source,l.lane,l.role;
```
Candidate detection frequencies by actual partner availability:
```sql
SELECT * FROM opportunity_detection_summary;
```
This view is cell-weighted descriptive, NOT the prior pair-balanced prediction metric. For valid comparison use the event-level Brier view with `event_weight`, equal source-appropriate group weighting and the frozen risk definitions. Never interpret number of risk rows as independent n.
```sql
SELECT lane,context,dependency_group,
       SUM(event_weight*(pool_loss-hybrid_loss))/SUM(event_weight) AS hybrid_gain
FROM opportunity_event_brier
GROUP BY lane,context,dependency_group;
```
Controlled isolate assays linked to whole-community endpoint:
```sql
SELECT * FROM opportunity_pairwise_assay_bridge;
```
This matches the source origin label, experimental medium and well position. It does not equate the isolated champion with the community dominant or a unique biological species. Marker signature mass is relative read abundance; the arithmetic-mean monoculture reference is a descriptive composition reference, not a null for absolute growth. The table `opportunity_signed_interaction_effects` must remain empty until the appropriate evidence is verified. Detection, contact assignment and co-culture read composition cannot establish competition or facilitation on their own.

## Next independent source
Venturelli et al. (2018), *Molecular Systems Biology* e8157, DOI 10.15252/msb.20178157, dataset EV1 (published supplementary XLSX): 12 strains, 66 pair combinations, monoculture and time-resolved OD600/relative-read profiles, species OD-equivalent abundance estimated by multiplication. Relevant because it could support an experimentally measured directed **net biomass response** against a verified matched monoculture, rather than simply read-share. First verify exact per-unit initial dose, dilution schedule, timing and source replicate/OD alignment; compare same-start conditions only, preserve detection limits and units. A finite net response is not a per-capita gLV interaction coefficient. The article already models pairwise interactions and community assembly. No source file acquisition or new score for this candidate in v0.6.
