"""medir.py -- protocolo experimental de la Parte III.

Corre ./montecarlo REPS veces por cada p, toma la mediana de T(p) y genera
resultados.csv, tabla.md y speedup.png. Verifica ademas que los aciertos
sean identicos en todas las corridas (reproducibilidad).
Uso: python3 medir.py   (compilar antes montecarlo.c)
"""
import csv
import os
import re
import statistics
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 100_000_000
HILOS = [1, 2, 4, 5, 8, 10, 15]
REPS = 5


def correr(p):
    env = dict(os.environ, OMP_NUM_THREADS=str(p))
    salida = subprocess.run(["./montecarlo", str(N)], env=env, capture_output=True,
                            text=True, check=True).stdout
    aciertos = int(re.search(r"Aciertos: (\d+)", salida).group(1))
    tiempo = float(re.search(r"Tiempo: ([\d.]+)", salida).group(1))
    return aciertos, tiempo


correr(1)  # calentamiento: carga del binario y de la libreria, no se mide

tiempos = {}
aciertos_ref = None
for p in HILOS:
    corridas = [correr(p) for _ in range(REPS)]
    for aciertos, _ in corridas:
        aciertos_ref = aciertos_ref or aciertos
        assert aciertos == aciertos_ref, f"p={p}: {aciertos} != {aciertos_ref}"
    tiempos[p] = [t for _, t in corridas]
    print(f"p={p:2d}  mediana={statistics.median(tiempos[p]):.4f} s  {tiempos[p]}")

T1 = statistics.median(tiempos[1])
filas = []
for p in HILOS:
    Tp = statistics.median(tiempos[p])
    S = T1 / Tp
    kf = (1 / S - 1 / p) / (1 - 1 / p) if p > 1 else None  # Karp-Flatt
    filas.append((p, Tp, S, S / p, kf, tiempos[p]))

S4 = next(S for p, _, S, *_ in filas if p == 4)
f = (1 - 1 / S4) / (1 - 1 / 4)  # fraccion paralelizable despejada de Amdahl con p=4

with open("resultados.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["p", "mediana_T_s", "S", "E", "karp_flatt"] + [f"T{i+1}" for i in range(REPS)])
    for p, Tp, S, E, kf, ts in filas:
        w.writerow([p, f"{Tp:.4f}", f"{S:.3f}", f"{E:.3f}", "" if kf is None else f"{kf:.4f}"]
                   + [f"{t:.4f}" for t in ts])

with open("tabla.md", "w") as fh:
    fh.write(f"N = {N:.0e}, {REPS} corridas por p, mediana. f estimado con S(4): {f:.4f}\n\n")
    fh.write("| p | Mediana T(p) [s] | S(p) | E(p) | Karp-Flatt e |\n|---|---|---|---|---|\n")
    for p, Tp, S, E, kf, _ in filas:
        fh.write(f"| {p} | {Tp:.4f} | {S:.3f} | {E:.3f} | {'-' if kf is None else f'{kf:.4f}'} |\n")
    fh.write(f"\nAciertos en todas las corridas: {aciertos_ref} (pi = {4 * aciertos_ref / N:.6f})\n")

ps = [r[0] for r in filas]
plt.figure(figsize=(6, 4.5))
plt.plot(ps, ps, "--", color="gray", label="Ideal S(p) = p")
plt.plot(ps, [1 / ((1 - f) + f / p) for p in ps], ":", color="tab:orange",
         label=f"Amdahl, f = {f:.4f}")
plt.plot(ps, [r[2] for r in filas], "o-", color="tab:blue", label="Medido (mediana)")
plt.xlabel("p (hilos)")
plt.ylabel("Speedup S(p)")
plt.title(f"Monte Carlo para pi, N = {N:.0e}")
plt.xticks(ps)
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("speedup.png", dpi=150)
print(open("tabla.md").read())
