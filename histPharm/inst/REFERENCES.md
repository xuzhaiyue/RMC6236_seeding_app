# Methodological reference points

## d-chain / sequential dose-response benchmark

Koplev S, Longden J, Ferkinghoff-Borg J, Blicher P, Linding R. Dynamic rearrangement of cell states detected by systematic screening of sequential anticancer treatments. *Cell Reports*. 2017. DOI: 10.1016/j.celrep.2017.08.095.

Public code: `skoplev/d-chain`. Its public example data use the columns `Experiment, CellLine, Run, Plate, Pretreatment, Compound, Concentration, RelCount`; the C++ source describes AB, A and A0 experiment classes and a Bliss-like interval model.

## histPharm positioning

histPharm should be benchmarked against established sequential-response methods rather than described as the first sequential therapy framework. Its intended novelty is the reciprocal history-conditioned estimand, explicit interaction-versus-efficacy decomposition, effect-space mapping, measured-state held-out prediction, and delayed-fate layer.
