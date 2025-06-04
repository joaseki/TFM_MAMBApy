#!/usr/bin/env python3

import os
import re
import pandas as pd
from cobra.io import read_sbml_model
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties

#Substitute values depending on dataset and model
# ── CONFIG ──────────────────────────────────────────────────────────────────────
EXP_FILE       = "expression.csv"
OUTPUT_DIR     = "gimme_outputs"
WO_INTAKES_DIR = os.path.join(OUTPUT_DIR, "wo_intakes")
W_INT_DIR      = os.path.join(OUTPUT_DIR, "with_intakes")
MODELS_DIR     = "models"
NO_INTAKE_SBML = os.path.join(MODELS_DIR, "iJO1366_gerosa_gimme.xml")
TTG_RATIO      = 9999
OBJ_FRAC       = 0.8
FLUX_THRESH    = 0.8
BIOMASS_RXN    = "BIOMASS_Ec_iJO1366_core_53p95M"
GPR_CLEANUP    = re.compile(r"__COBAMPGPRDOT__[0-9]+")

# ── HELPERS ────────────────────────────────────────────────────────────────────
def clean_gene(x: str) -> str:
    return GPR_CLEANUP.sub("", x)

def zero_none(mapping):
    return {rid: (0 if v is None else v) for rid, v in mapping.items()}

# ── PREP ───────────────────────────────────────────────────────────────────────
os.makedirs(WO_INTAKES_DIR, exist_ok=True)
os.makedirs(W_INT_DIR, exist_ok=True)

expr_df = pd.read_csv(EXP_FILE, index_col=0)
samples = expr_df.columns.tolist()

print("Starting pipeline…")
for sample in samples:
    print(f"\n▶ Sample: {sample}")
    # build per‐sample OmicsContainer
    samp = expr_df[[sample]].T
    samp.index = [sample]
    container = TabularReader(
        path_or_df=samp,
        nomenclature="entrez_id",
        omics_type="transcriptomics"
    ).to_containers()[0]

    for case, sbml, outdir in [
        ("no_intake", NO_INTAKE_SBML, WO_INTAKES_DIR),
        ("intake",    os.path.join(MODELS_DIR, f"iJO1366_gerosa_gimme_{sample}.xml"), W_INT_DIR)
    ]:
        print(f"  • Case: {case}")
        if not os.path.isfile(sbml):
            print(f"    ⚠️ Missing model: {os.path.basename(sbml)}")
            continue

        # load and wrap
        model = read_sbml_model(sbml)
        wrap  = ReconstructionWrapper(
            model=model.copy(),
            ttg_ratio=TTG_RATIO,
            gpr_gene_parse_function=clean_gene
        )
        r_ids   = wrap.model_reader.r_ids
        obj_idx = r_ids.index(BIOMASS_RXN)

        # map expression → reaction scores
        data_map   = container.get_integrated_data_map(wrap.model_reader, and_func=min, or_func=sum)
        integrator = ContinuousScoreIntegrationStrategy(score_apply=zero_none)
        scores     = integrator.integrate(data_map=data_map)

        # build GIMME props
        props = GIMMEProperties(
            exp_vector=list(scores.values()),
            obj_frac=OBJ_FRAC,
            objectives=[{obj_idx: 1}],
            preprocess=True,
            flux_threshold=FLUX_THRESH,
            solver="GUROBI",
            reaction_ids=r_ids,
            metabolite_ids=wrap.model_reader.m_ids
        )

        # run GIMME
        gm = GIMME(S=wrap.S, lb=wrap.lb, ub=wrap.ub, properties=props)
        active_idx = gm.run()
        if active_idx is None:
            print(f"  GIMME infeasible")
            continue

        # prune inactive
        active_set = {r_ids[i] for i in active_idx}
        for rx in model.reactions:
            if rx.id not in active_set:
                rx.knock_out()

        # FBA
        sol    = model.optimize()
        fluxes = sol.fluxes

        # save as CSV
        flux_df = fluxes.reset_index()
        flux_df.columns = ["reaction", "flux"]
        csv_out = os.path.join(outdir, f"{sample}_gimme_solution.csv")
        flux_df.to_csv(csv_out, index=False)
        print(f" Saved  {os.path.basename(csv_out)}")

print("\n Pipeline complete.")

