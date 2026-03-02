### Supplementary Table 1: Complete parameter list of the NED model

| Parameter | Module | Species / Process | Value in notebook | Typical range | Status | Basis | Section / Code | Rationale / interpretation |
|-----------|--------|-------------------|-------------------|---------------|--------|-------|----------------|-----------------------------|
| **$NICD$ module** |||||||||
| `Vmni` | NICD_reg | $NICD$ max production (Hp) | 25 | — | fixed | phenomenological | §1.1 | With `Kdni=0.25`, the effective $NICD$ ceiling is `Vmni/Kdni ≈ 100`; `Vmni` is fixed while `Kdni` modulates timescale. |
| `hMni` | NICD_reg | $NICD$ activation K₀.₅ | 2000 | — | fixed | phenomenological | §1.1 | Sensitivity to total productive input (`tmDN + cmDN`, typically 0–5000). |
| `Kdni` | NICD_reg | $NICD$ degradation rate | 0.25 | 0.1–0.5 | scanned | experimentally constrained | §1.1 | Sets $NICD$ relaxation time while preserving ~0–100 dynamic range set by `Vmni/Kdni`. |
| `tmDN` | CC3D → $NICD$ | Trans $DLL4$–$NOTCH1$ interactions | computed each MCS | 0–5000 | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | Contact-area–weighted productive trans interactions from CC3D. |
| `cmDN` | cis_interactions | Productive cis $DLL4$–$NOTCH1$ | computed | 0–5000 | dynamic | structural (coupling) | §2.4 | Productive cis contribution to $NICD$ within the SBML layer. |
| **Notch / ligand expression** |||||||||
| `bN1` | ligands_receptors_reg | Basal $NOTCH1$ production | 250 | — | fixed | phenomenological | §2.1 | Sets $N1$ production (lumped transcription/translation); steady state ~5,000 supports the 10³–10⁴ regime of Boareto et al. PNAS 2015 |
| `Kdn1` | ligands_receptors_reg | $NOTCH1$ degradation | 0.05 | — | fixed | phenomenological/literature constrained | §2.1 | τ≈20 min, t½≈14 min consistent with rapid receptor turnover in activated context. |
| `hMn1` | ligands_receptors_reg | $NOTCH1$ inducibility K₀.₅ | 50 | — | fixed | phenomenological | §2.1 | Half-max $HES1$ for optional $N1$ inducibility. |
| $N1$ | ligands_receptors_reg | $NOTCH1$ level | dynamic | 0–5000 | dynamic | structural (state) | §2.1 | Receptor abundance updated by Antimony + CC3D terms. |
| `bD4` | ligands_receptors_reg | Basal $DLL4$ production | 250 | — | fixed | phenomenological | §2.2 | Fully derepressed $DLL4$ flux; steady state ~5,000 in derepressed state. |
| `hMd4` | ligands_receptors_reg | $DLL4$ repression K₀.₅ | 50 | — | fixed | phenomenological | §2.2 | Half-max $HES1$ for negative regulation. |
| `Kdd4` | ligands_receptors_reg | $DLL4$ degradation | 0.05 | — | fixed | phenomenological/literature constrained | §2.2 | τ≈20 min, t½≈14 min. |
| `dD4` | CC3D → Antimony | Engaged $DLL4$ loss term | computed each MCS | 0–$D4$-dependent | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | Accumulates trans-engaged ligand per CC3D step. |
| $D4$ | ligands_receptors_reg | $DLL4$ level | dynamic | 0–5000 | dynamic | structural (state) | §2.2 | Ligand abundance after production/repression and losses. |
| `VmJ1` | ligands_receptors_reg | Max $Jagged1$ induction | 250 | — | fixed | phenomenological | §2.3 | Max inducible $Jagged1$ flux; steady state ~5,000 when fully induced. |
| `hMj1` | ligands_receptors_reg | $Jagged1$ induction K₀.₅ | `ps_ja` (≈50) | 25–100 | scanned | phenomenological | §2.3 | Half-max $HES1$ for lateral induction. |
| `bJ1` | ligands_receptors_reg | Basal $Jagged1$ production | 0 | — | fixed | phenomenological | §2.3 | Default off; can be enabled for basal $Jagged1$. |
| `Kdj1` | ligands_receptors_reg | $Jagged1$ degradation | 0.05 | — | fixed | phenomenological/literature constrained | §2.3 | τ≈20 min, t½≈14 min. |
| `dJ1` | CC3D → Antimony | Engaged $Jagged1$ loss term | computed each MCS | 0–$J1$-dependent | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | Accumulates trans-engaged $Jagged1$ per CC3D step. |
| $J1$ | ligands_receptors_reg | $Jagged1$ level | dynamic | 0–5000 | dynamic | structural (state) | §2.3 | $Jagged1$ abundance after induction and losses. |
| **Trans interactions** |||||||||
| `tmJN` | CC3D → $NICD$ | Trans $Jagged1$–$NOTCH1$ | computed each MCS | 0–5000 | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | Analogous to `tmDN`; used where $Jagged1$ trans engagement is modeled. |
| `Rd4` | CC3D | Reacted $DLL4$ | computed | 0–5000 | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | $Dll4$ consumed in trans interactions. |
| `Rj1` | CC3D | Reacted $Jagged1$ | computed | 0–5000 | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | $Jagged1$ consumed in trans interactions. |
| **Cis interactions** |||||||||
| `cJN` | cis_interactions | Non-productive $Jagged1$–Notch cis | mass-action | n/a | fixed law | structural | §2.4 | Canonical cis-inhibition via $Jagged1$. |
| `cdn` | cis_interactions | Productive $DLL4$–Notch cis | mass-action | n/a | fixed law | structural | §2.4 | Allows $NICD$ production from cis complexes. |
| `cadn` | cis_interactions | Active cis $DLL4$–Notch | derived | 0–5000 | dynamic | structural (state) | §2.4 | Intermediate contributing to $NICD$ production. |
| **Competition parameters** |||||||||
| `kpDJ` | CC3D | $DLL4$/$Jagged1$ partitioning | computed each MCS | 0–1 | dynamic | phenomenological control | EC_Connect_core.py:203-235 | Allocates receptor usage between $DLL4$ and $Jagged1$. |
| `KJC` | ligands_receptors_reg | $Jagged1$ affinity modifier | 1.0 | 0.5–2 | fixed | phenomenological control | EC_Connect_core.py:203-235 | Equal-affinity baseline; can modulate $Jagged1$ competitiveness. |
| **$HES1$ – no autoregulation** |||||||||
| `VmH1` | HES_no_auto | Max $HES1$ production | 10 | — | fixed | phenomenological | §3.1 | Sets $HES1$ scale in reduced model (0–~10 corresponds to $NICD$ 0–100). |
| `hMh1` | HES_no_auto | $NICD$→$HES1$ K₀.₅ | 35 | — | fixed | phenomenological | §3.1 | Matches $NICD$ range to $HES1$ activation threshold. |
| `Kdh1` | HES_no_auto | $HES1$ degradation | 0.09 | — | fixed | experimentally constrained | §3.1 | Places $HES1$ on hour-scale. |
| `bH1` | HES_no_auto | Basal $HES1$ production | 0.0 | — | fixed | phenomenological | §3.1 | Default zero; enables pure $NICD$ dependence. |
| **$HES1$ – autoregulation** |||||||||
| `vm_prna` | HES_auto | Max pre-$mRNA$ transcription | 0.1 | 0.05–0.2 | scanned | phenomenological | §3.2.2 | $NICD$-driven transcription strength under feedback. |
| `K` | HES_auto | Shared K₀.₅ (Hill/MM) | 0.5 | 0.3–0.7 | scanned | phenomenological | §3.2.2 | Common K₀.₅ / K_m for Hill and MM terms. |
| `hh1` | HES_auto | $HES1$ autorepression Hill | 8 | 6–10 | fixed | structural required | §3.2.2 | Cooperativity needed for oscillations. |
| `kmrp` | HES_auto | $mRNA$ processing rate | 0.1 | 0.05–0.2 | scanned | structural required | §3.2.2 | Part of distributed transcriptional delay. |
| `kmrd` | HES_auto | $mRNA$ degradation | 0.05 | 0.02–0.1 | scanned | structural required | §3.2.2 | Controls transcript lifetime and oscillation stability. |
| `kphp` | HES_auto | Translation rate | 0.1 | 0.05–0.2 | scanned | structural required | §3.2.2 | Delay component between $mRNA$ and protein. |
| `khp` | HES_auto | Protein maturation | 0.1 | 0.05–0.2 | scanned | structural required | §3.2.2 | Nuclear/active conversion delay. |
| `khd` | HES_auto | Active $HES1$ degradation | 0.075 | 0.05–0.15 | scanned | experimentally constrained | §3.2.2 | Rapid turnover required for sustained oscillations. |
| $NICD$ | HES_auto | $NICD$ input (scaled) | 0–1 | 0–1 | controlled | experimental driver | §3.2.3 | Bifurcation parameter for oscillatory regime (scaled from observed $NICD$ range). |
| **Utility / kinetics** |||||||||
| `Hp()` | utility | Positive Hill | — | — | structural | — | §0 (utility) | Canonical sigmoidal activation. |
| `Hc()` | utility | Competitive Hill | — | — | structural | — | §0 (utility) | Implements $NICD$–$HES1$ competition. |
| `Ma()` | utility | Mass action | — | — | structural | — | §0 (utility) | Linear degradation / interaction law. |
| `MM()` | utility | Michaelis–Menten | — | — | structural | — | §0 (utility) | Saturable processing and degradation. |

Notes:
- $NICD$ amplitude is calculated by $Vmni/Kdni$ (~100 with 25/0.25); scanning $Kdni$ changes the timescale while preserving this range.
- $N1$, $D4$, and $J1$ basals/decays place steady states in the 10³–10⁴ regime (lumped units), consistent with Boareto et al PNAS 2015.
- Time unit in the Antimony layer is 1 min; characteristic times (e.g., 1/$Kdn1$, 1/$Kdd4$, 1/$Kdj1$) are in real minutes.
- Typical ranges are provided where parameters are scanned or vary during simulations; fixed entries list “—”.
