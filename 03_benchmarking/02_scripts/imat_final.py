#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.imat import IMAT, IMATProperties

#change according to the dataset analysed and its corresponding model
# ── CONFIG ──────────────────────────────────────────────────────────────────────
BASE_MODEL        = "models/iJO1366_gerosa_iMAT.xml"
INTAKE_MODEL_TPL  = "models/iJO1366_gerosa_iMAT_{sample}.xml"
EXP_FILE          = "expression.csv"
OUTPUT_DIR        = "imat_fluxes_gerosa"
WO_DIR            = os.path.join(OUTPUT_DIR, "wo_intakes")
WI_DIR            = os.path.join(OUTPUT_DIR, "with_intakes")
TTG_RATIO         = 9999
DISC_SIGMA_FACTOR = 0.5          # ±0.5σ discretization
IMAT_THRESHOLDS   = (-0.99, 0.99)
GPR_CLEANER       = re.compile(r"__COBAMPGPRDOT__[0-9]+")

# ── HELPERS ─────────────────────────────────────────────────────────────────────
def discretize(s: pd.Series) -> pd.Series:
    μ, σ = s.mean(), s.std()
    up, lo = μ + DISC_SIGMA_FACTOR*σ, μ - DISC_SIGMA_FACTOR*σ
    return s.apply(lambda v: 1 if v > up else (-1 if v < lo else 0))

def clean_gpr(g: str) -> str:
    return GPR_CLEANER.sub("", g)

def zero_none(m: dict) -> dict:
    return {rid: (0 if (v is None or np.isnan(v)) else v)
            for rid, v in m.items()}

# ── PREP ────────────────────────────────────────────────────────────────────────
for d in (WO_DIR, WI_DIR):
    os.makedirs(d, exist_ok=True)

print("Loading expression matrix…")
expr_df = pd.read_csv(EXP_FILE, index_col=0)

# ── MAIN ────────────────────────────────────────────────────────────────────────
for sample in expr_df.columns:
    print(f"\n▶ Sample: {sample}")
    # discretize one column
    disc = discretize(expr_df[sample])
    disc_df = pd.DataFrame(disc).T
    disc_df.index = [sample]
    container = TabularReader(
        path_or_df=disc_df,
        nomenclature="entrez_id",
        omics_type="transcriptomics"
    ).to_containers()[0]

    for scenario, model_path, out_dir in (
        ("wo_intakes",    BASE_MODEL,        WO_DIR),
        ("with_intakes",  INTAKE_MODEL_TPL.format(sample=sample),  WI_DIR)
    ):
        if not os.path.isfile(model_path):
            print(f"   Model not found for {scenario}: {model_path}")
            continue

        print(f"  • Running IMAT ({scenario})")
        # wrap model
        m = read_sbml_model(model_path)
        wrap = ReconstructionWrapper(
            model=m,
            ttg_ratio=TTG_RATIO,
            gpr_gene_parse_function=clean_gpr
        )

        # reaction scoring
        data_map   = container.get_integrated_data_map(
            model_reader=wrap.model_reader,
            and_func=min,
            or_func=sum
        )
        integrator = ContinuousScoreIntegrationStrategy(score_apply=zero_none)
        scores     = integrator.integrate(data_map=data_map)

        # run IMAT
        props = IMATProperties(
            exp_vector=list(scores.values()),
            exp_thresholds=IMAT_THRESHOLDS
        )
        imat = IMAT(S=wrap.S, lb=wrap.lb, ub=wrap.ub, properties=props)
        sol  = imat.run_imat()
        if sol is None:
            print(f"     No IMAT solution for {scenario}")
            continue

        # extract + pad flux vector
        flux_vec     = sol.x()
        reaction_ids = wrap.model_reader.r_ids
        full_flux = [
            flux_vec[i] if i < len(flux_vec) else 0.0
            for i in range(len(reaction_ids))
        ]

        # save CSV
        out_csv = os.path.join(out_dir, f"{sample}_imat_fluxes.csv")
        pd.DataFrame({
            "reaction": reaction_ids,
            "flux":     full_flux
        }).to_csv(out_csv, index=False)
        print(f"     Saved  {os.path.basename(out_csv)}")



