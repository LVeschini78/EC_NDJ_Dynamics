### Supplementary Table 1: Complete parameter list of the NED model

| Parameter | Module | Species / Process | Value in notebook | Typical range | Status | Basis | Section / Code | Rationale / interpretation |
|-----------|--------|-------------------|-------------------|---------------|--------|-------|----------------|-----------------------------|
| **Model control** |||||||||
| `h1_auto` | model_control | Pathway selection switch | 0 | 0 or 1 | fixed per run | structural | §3 | As implemented in equations: `h1_auto=1` enables the auto-regulatory $HES1$ production path (Eqs 11-17), while `h1_auto=0` enables the standard model without feedback (eqs 9-10). |
| **PS scanned controls (declared scan knobs)** |||||||||
| `ps_bd4` | scan_controls | Basal $DLL4$ production scaler | 0.0 | 50-150 (search scripts) | scanned (workflow) | phenomenological | NED_model_string.py, EC_Connect_search.py | Declared as scan control; in this table's model string `bD4` is currently fixed to 250, so this control is not directly wired there. |
| `ps_bn1` | scan_controls | Basal $NOTCH1$ production scaler | 0.0 | 50-150 (search scripts) | scanned (workflow) | phenomenological | NED_model_string.py, EC_Connect_search.py | Declared as scan control; in this table's model string `bN1` is currently fixed to 250, so this control is not directly wired there. |
| `ps_bJ1` | scan_controls | Basal $Jagged1$ production scaler | 0.0 | 0-150 (search scripts) | scanned (workflow) | phenomenological | NED_model_string.py, EC_Connect_search.py | Declared as scan control; in this table's model string `bJ1` is fixed to 0. Legacy code may use lowercase alias `ps_bj1`. |
| `ps_J1ind` (`ps_ind` alias) | scan_controls | Inducible $Jagged1$ amplitude control | 0.0 | context dependent | scanned | phenomenological | NED_model_string.py, EC_Connect_core.py | This is the key "PS_ind" knob: in this model string `VmJ1 := amp * ps_J1ind`; legacy CC3D strings often use alias `ps_ind`. |
| `ps_ja` | scan_controls | $Jagged1$ induction threshold control | 0.0 | typically 50-5000 | scanned | phenomenological | NED_model_string.py, EC_Connect_search.py | Wired directly via `hMj1 := ps_ja` (half-max $HES1$ for Jagged induction). |
| `ps_ni` | scan_controls | $NICD$ activation threshold control | 2000 | typically 2000-200000 | scanned | phenomenological | NED_model_string.py, EC_Connect_search.py | Wired directly via `hMni := ps_ni`. |
| `ps_Kdni` | scan_controls | $NICD$ degradation control | 0.0 | typically 0-100 (search scripts) | scanned | phenomenological | NED_model_string.py, EC_Connect_search.py | Wired as `Kdni := ps_Kdni/4` in this model string (note the scaling by 1/4). |
| `ps_kconv` | scan_controls | Productive cis conversion control | 0.0 | typically 0-15 (search scripts) | scanned | phenomenological | NED_model_string.py, EC_Connect_search.py | Wired directly via `kconv := ps_kconv`; controls cDN to cmDN conversion strength. |
| **$NICD$ module** |||||||||
| `Vmni` | NICD_reg | $NICD$ max production (Hp) | 25 | — | fixed | phenomenological | §1.1 | With `Kdni=0.25`, the effective $NICD$ ceiling is `Vmni/Kdni ≈ 100`; `Vmni` is fixed while `Kdni` modulates timescale. |
| `hMni` | NICD_reg | $NICD$ activation K0.5 | `ps_ni` (default 2000) | context dependent | scanned | phenomenological | §1.1 | Sensitivity to total productive input (`tmDN + cmDN`, typically 0-5000). |
| `Kdni` | NICD_reg | $NICD$ degradation rate | `ps_Kdni/4` (default 0) | context dependent | scanned | experimentally constrained | §1.1 | Directly controlled by scan knob `ps_Kdni`; effective relaxation and amplitude depend on scanned value. |
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
| `VmJ1` | ligands_receptors_reg | Max $Jagged1$ induction | `amp * ps_J1ind` | context dependent | scanned | phenomenological | §2.3 | Inducible $Jagged1$ amplitude set by PS_ind control (`ps_J1ind`; alias `ps_ind` in some code paths). |
| `hMj1` | ligands_receptors_reg | $Jagged1$ induction K0.5 | `ps_ja` | context dependent | scanned | phenomenological | §2.3 | Half-max $HES1$ for lateral induction. |
| `bJ1` | ligands_receptors_reg | Basal $Jagged1$ production | 0 | — | fixed | phenomenological | §2.3 | Default off; can be enabled for basal $Jagged1$. |
| `Kdj1` | ligands_receptors_reg | $Jagged1$ degradation | 0.05 | — | fixed | phenomenological/literature constrained | §2.3 | τ≈20 min, t½≈14 min. |
| `dJ1` | CC3D → Antimony | Engaged $Jagged1$ loss term | computed each MCS | 0–$J1$-dependent | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | Accumulates trans-engaged $Jagged1$ per CC3D step. |
| $J1$ | ligands_receptors_reg | $Jagged1$ level | dynamic | 0–5000 | dynamic | structural (state) | §2.3 | $Jagged1$ abundance after induction and losses. |
| **Trans interactions** |||||||||
| `tmJN` | CC3D → $NICD$ | Trans $Jagged1$–$NOTCH1$ | computed each MCS | 0–5000 | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | Analogous to `tmDN`; used where $Jagged1$ trans engagement is modeled. |
| `Rd4` | CC3D | Reacted $DLL4$ | computed | 0–5000 | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | $Dll4$ consumed in trans interactions. |
| `Rj1` | CC3D | Reacted $Jagged1$ | computed | 0–5000 | dynamic | structural (coupling) | EC_Connect_core.py:203-235 | $Jagged1$ consumed in trans interactions. |
| **Cis interactions** |||||||||
| `cr` | cis_interactions | Cis interaction rate scaling | 0 | 0-1 | scanned | phenomenological | §2.4 | Scales the flux of cis $DLL4$-$N1$ and $Jagged1$-$N1$ complex formation; `cr=0` disables cis interactions entirely. |
| `crd` | cis_interactions | Effective cis rate (derived) | `cr*(1-tmDN/(300+tmDN))` | 0-`cr` | dynamic | structural (coupling) | §2.4 | Trans interactions suppress cis formation competitively; `crd` decreases as `tmDN` rises, implementing the trans-cis competition. |
| `kcd` | cis_interactions | Cis complex degradation | 0.10 | 0.05–0.2 | fixed | phenomenological | §2.4 | Sets turnover of unproductive cis complexes (`cDN`, `cJN`); τ ≈ 10 min. |
| `kconv` (`ps_kconv`) | cis_interactions | cDN → cmDN conversion rate | 0 | 0–0.5 | scanned | phenomenological | §2.4 | Rate at which the cis $DLL4$–$N1$ complex is converted to the transcriptionally active form `cmDN` that contributes to $NICD$ production (Eq 18). |
| `cDN` | cis_interactions | Cis $DLL4$–$N1$ complex | dynamic | 0–5000 | dynamic | structural (state) | §2.4 | Formed by mass-action (`Ma(cDNi,crd)`); depleted by conversion to `cmDN` and degradation. |
| `cmDN` | cis_interactions | Active cis $DLL4$–$N1$ | derived | 0–5000 | dynamic | structural (coupling) | §2.4 | Productive cis intermediate; feeds into $NICD$ production alongside `tmDN` (Eq 1). |
| `cJN` | cis_interactions | Cis $Jagged1$–$N1$ complex | dynamic | 0–5000 | dynamic | structural (state) | §2.4 | Non-productive cis complex that sequesters receptor; formed by `Ma(cJNi,crd)` and degraded at rate `kcd`. |
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
| `Hc()` | utility | Competitive Hill | — | — | structural | — | §0 (utility) | Implements $NICD$–$HES1$ competition in autoregulatory `pmRNA` production. |
| `Hn()` | utility | Negative Hill | — | — | structural | — | §0 (utility) | Sigmoidal repression: `bp/(1+(S/hM)^h)`; used in $DLL4$ repression by $HES1$ (Eq 5). |
| `Ma()` | utility | Mass action | — | — | structural | — | §0 (utility) | Linear degradation / interaction law. |
| `MM()` | utility | Michaelis–Menten | — | — | structural | — | §0 (utility) | Saturable processing and degradation. |

Notes:
- In this model string, Parameter Scan knobs are explicitly declared (`ps_*`).
- $NICD$ amplitude follows $Vmni/Kdni$ when `Kdni>0`; because `Kdni := ps_Kdni/4`, amplitude/timescale depend on the chosen scan setting.
- $N1$, $D4$, and $J1$ basals/decays place steady states in the 10³–10⁴ regime (lumped units), consistent with Boareto et al PNAS 2015.
- Time unit in the Antimony layer is 1 min; characteristic times (e.g., 1/$Kdn1$, 1/$Kdd4$, 1/$Kdj1$) are in real minutes.
- Typical ranges are provided where parameters are scanned or vary during simulations; fixed entries list “—”.
