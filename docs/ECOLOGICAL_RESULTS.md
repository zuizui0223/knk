# 生態学的成果台帳 — Fungal and grassland assembly history

> 2026-11-22：EOGの元の数理・既知真値・予測性能とは区別。KNKには独立操作実験2系統と、同一データを再解析したv31–v37、STOPのv27/v29、および採点前の草地候補をまとめる。

## 核心

同じ実験的到着履歴が、現在の群集に**一様には保持されない**。何を測定対象（estimand）にするかで保持量が変わる。これは生態学的な現象であり、EOGの性能優位を主張するものではない。

## 独立操作実験（2系統のみ）

### Leopold / Populus 系・菌類（v28）

- n=233植物、12宿主遺伝子型×5到着順序（60組合せすべて含む）。
- 菌類5群組成：観測部分R²=0.3697、genotype層別置換帰無中央値0.2230、過剰 E=+0.1467、p=0.0001（9,999置換）。
- さび病病斑率：R²=0.2300、帰無中央値0.2147、E=+0.0153、p=0.3893。
- 結論：履歴の痕跡は菌類の組成に強く、病斑率へは同じ形では保持されない。さび病への履歴効果が存在しない、とは言えない。
- [原典の凍結結果](eog_original_idea_leopold_history_retention_v28_result.md)。

### 草地・地上部/地下部（v30）

- n=19 rhizobox、F/G/L-firstの3到着履歴、各replicate内9,999置換。
- 地上部forb/grass/legume組成：観測部分R²=0.9529、帰無中央値0.1495、E=+0.8034、p=0.0001。
- 6層根バイオマスの縦方向分布：R²=0.1853、帰無中央値0.1509、E=+0.0344、p=0.3765。
- 総地上部バイオマス：E=+0.0948、p=0.2562。独立の31,104ラベル完全列挙でも組成と根の解釈は変わらない。
- 結論：履歴は主に**何の機能群がバイオマスを担うか**に残り、総量には同じ強さで残らない。
- [原典の凍結結果](eog_grassland_above_below_history_retention_v30_result.md)。

## 同一実験の後続研究：独立反復ではない

- **v31:** 菌類のidentity-free rank abundance E=+0.1486, p=0.0001。草地では E=+0.1949, p=0.0769。[詳細](eog_identity_storage_rank_abundance_v31_result.md)
- **v32:** first-arriver role対応で見かけ上大きい差が出たが、v33の全対応placeboにより**その機構解釈は失効**。[詳細](eog_arrival_role_symmetry_v32_result.md)
- **v33:** 菌類の正しいrole mappingは120候補中53位（tail 0.4417）、草地は6候補中最下位。first-arriver role symmetryを特定機構とする主張を反証。[詳細](eog_role_mapping_specificity_v33_result.md)
- **v34:** Alternaria履歴を除くとrank-memory Eが+0.1486→+0.0875、差+0.0611。残った4履歴でもp=0.0086で履歴保持。[詳細](eog_history_level_leverage_v34_result.md)
- **v35:** Alternaria削除効果はrank-1 dominance +0.08647、下位順位+0.03957で前者が大きい。[詳細](eog_dominance_tail_history_storage_v35_result.md)
- **v36:** Alternaria対他4履歴でrank-1優占 E=+0.11812, p=0.0002、rank2–5 E=+0.04723, p=0.0316。[詳細](eog_alternaria_binary_components_v36_result.md)
- **v37:** Alternaria到着履歴によるrank-1増大は12遺伝子型中11型で正。平均+0.06099、p=0.0001。genotype-specific interaction p=0.7548で不支持。[詳細](eog_alternaria_genotype_generality_v37_result.md)

## STOPおよび未採点資料を肯定結果に含めない

- **v27:** Catano修復草地は再利用許諾とspecies/trait NMDSのplot ID欠如でSTOP。履歴効果は採点していない。[STOP文書](eog_original_idea_empirical_history_retention_v27_result.md)
- **v29:** 木材分解菌の元CSVの正規バイトを取得できずSTOP。組成・分解への効果量は存在しない。[STOP文書](eog_wood_decomposer_history_retention_v29_result.md)
- **Czech dry-grassland:** データ源・スキーマ・幾何の審査段階で、生態アウトカムは未採点。[契約](czech_grassland_geometry_contract.md)

## 未実証の機構・一般化

- 優占を作る生理的・分子機構、niche preemptionとniche modificationの識別。
- 菌類・草地2系統を超えた一般則。
- 組成変化が病徴・地下部へ媒介する因果パス。
- 名前なし順位解析による首位の菌種の同定。
- v32のrole interchangeability機構はv33で否定されているため復活させない。

独立の EOG-WF log-loss 検証・BAM known-truth 768条件は[方法論 EOG](https://github.com/zuizui0223/eog)の成果であり、この生態成果台帳には算入しません。
