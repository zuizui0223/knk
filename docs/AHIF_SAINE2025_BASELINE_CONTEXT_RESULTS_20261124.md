# Saine: initial community state predicts rank-response modestly, not causal moderation

Source-verified run https://github.com/zuizui0223/knk/actions/runs/45388722290, job 138426416400, commit 1a7be170c060b03cdcfcc018119895a7e2f43900; all tests, original RData/author-source digest checks and artifact upload passed. 446 complete logs from 460 originally allocated (347 treated, 99 control); same source experiment. No independent study count or formal meta-analysis eligibility increment.

| Training year | Endpoint | Nested gain over context-only | Context-only gain over zero | Sharp-null p | Within-family Holm | Improved held sites |
|---|---|---:|---:|---:|---:|---:|
| 2020 | named_TV | -0.2698% | 0.3721% | 0.772 | 0.772 | 0/5 |
| 2020 | rank_TV | 1.3348% | -0.1198% | 0.012 | 0.024 | 4/5 |
| 2021 | named_TV | -0.1606% | 1.0203% | 0.539 | 0.539 | 0/5 |
| 2021 | rank_TV | 1.8455% | 0.9179% | 0.013 | 0.026 | 3/5 |

Rank-TV predictive gain 1.335% across year/site, 1.846% same year across site. Named-TV gain negative in both lanes. Held-site gains not uniform; report all folds in the JSON. Both families use the same 2021 endpoints/logs, not independent confirmation. The two Holm families were prespecified separately; if adjusted across all four sharp-null tests posthoc, both rank results yield Holm .048. The reference shuffles binary labels under sharp no introduction effect, NOT a test of no effect modification allowing nonzero constant effects. Positive gains are modest exploratory prediction leads, not evidence of direct competitive coefficients, proof of rank buffering, or proof of causal moderators. Baseline read depth may be predictive; separate measurement from biological baseline variables. Natural/felled context already included in the baseline model. Site held out for fit and tuning; 2021 responses in other sites may tune alpha. Post-treatment complete-case selection and only five fixed sites remain limitations.

Next fixed falsification: compare depth-plus-log-context model with depth-plus-three-biological-summary model using same nested site splits; approximate restricted-model wild bootstrap under a nonzero introduction effect allowed to vary with depth/context, plus direct joint robust Wald test as a distinct diagnostic. Keep all results and avoid tuning after new scores. No change to 161-effect ledger or frozen manuscript. Source aggregate validation/ahif_v1/saine2025_baseline_context_v1.json.
