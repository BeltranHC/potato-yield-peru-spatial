"""MGWR con inferencia exacta (Oshan et al. 2019, paquete mgwr).

Uso: python py/mgwr_inferencia.py <entrada.csv> <carpeta_salida>
Entrada: columnas y, X_km, Y_km y las covariables (ya estandarizadas en R).
Salidas (CSV): mgwr_local.csv (betas, EE, t, número de condición local por distrito),
               mgwr_global.csv (bw, IC del bw, ENP_j, t crítico ajustado por covariable),
               mgwr_diag.csv (AICc, R2, ENP total de GWR y MGWR).
"""
import sys
import numpy as np
import pandas as pd
from mgwr.sel_bw import Sel_BW
from mgwr.gwr import GWR, MGWR

entrada, salida = sys.argv[1], sys.argv[2]
d = pd.read_csv(entrada)
xs = [c for c in d.columns if c not in ("ubigeo", "y", "X_km", "Y_km")]
coords = d[["X_km", "Y_km"]].values
y = d[["y"]].values
X = d[xs].values

# GWR (un solo ancho de banda), kernel bicuadrado adaptativo, AICc
sel_g = Sel_BW(coords, y, X, kernel="bisquare", fixed=False)
bw_g = sel_g.search(criterion="AICc")
gwr = GWR(coords, y, X, bw_g, kernel="bisquare", fixed=False).fit()

# MGWR (un ancho de banda por covariable), backfitting iniciado en el bw de GWR
sel_m = Sel_BW(coords, y, X, multi=True, kernel="bisquare", fixed=False)
sel_m.search(criterion="AICc", multi_bw_min=[30], init_multi=bw_g, verbose=False)
mgwr = MGWR(coords, y, X, sel_m, kernel="bisquare", fixed=False).fit()

nombres = ["Intercept"] + xs
tcrit = mgwr.critical_tval()                     # ajustado por ENP_j (da Silva & Fotheringham 2016)
try:
    intervalos = mgwr.get_bws_intervals(sel_m)
except Exception:                                # noqa: BLE001
    intervalos = [(np.nan, np.nan)] * len(nombres)
cn, vdp = mgwr.local_collinearity()

pd.DataFrame({
    "variable": nombres, "bw": sel_m.bw[0], "bw_ic_inf": [i[0] for i in intervalos],
    "bw_ic_sup": [i[1] for i in intervalos], "ENP_j": mgwr.ENP_j, "t_critico": tcrit,
}).to_csv(f"{salida}/mgwr_global.csv", index=False)

loc = pd.DataFrame({"ubigeo": d["ubigeo"]})
for j, v in enumerate(nombres):
    loc[v] = mgwr.params[:, j]
    loc[f"{v}_SE"] = mgwr.bse[:, j]
    loc[f"{v}_TV"] = mgwr.tvalues[:, j]
loc["cond_local"] = np.asarray(cn).ravel()
loc["residual"] = mgwr.resid_response.ravel()
loc.to_csv(f"{salida}/mgwr_local.csv", index=False)

pd.DataFrame({
    "modelo": ["GWR", "MGWR"], "bw_unico": [bw_g, np.nan],
    "AICc": [gwr.aicc, mgwr.aicc], "R2": [gwr.R2, mgwr.R2], "adjR2": [gwr.adj_R2, mgwr.adj_R2],
    "ENP": [gwr.ENP, mgwr.ENP],
    "cond_local_max": [np.nan, float(np.max(cn))], "pct_cond_mayor_30": [np.nan, float(np.mean(np.asarray(cn) > 30) * 100)],
}).to_csv(f"{salida}/mgwr_diag.csv", index=False)
print("bw MGWR:", dict(zip(nombres, sel_m.bw[0])))
print("ENP_j:", np.round(mgwr.ENP_j, 2))
print("AICc GWR/MGWR:", round(gwr.aicc, 1), round(mgwr.aicc, 1), "| R2:", round(gwr.R2, 3), round(mgwr.R2, 3))
print("Número de condición local: máx", round(float(np.max(cn)), 1))
