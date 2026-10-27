"""medir_schedule.py -- II.5: compara planificaciones con el protocolo de la Parte III.

Para cada planificacion y cada p corre ./montecarlo REPS veces y reporta la
mediana de T(p). Escribe schedule.md.
Uso: python3 medir_schedule.py   (compilar antes montecarlo.c)
"""
import os
import re
import statistics
import subprocess

HILOS = [1, 2, 4, 5, 8, 10, 15]
PLANES = ["static", "dynamic,1048576", "dynamic,1"]
REPS = 5


def tiempo(p, plan):
    env = dict(os.environ, OMP_NUM_THREADS=str(p), OMP_SCHEDULE=plan)
    salida = subprocess.run(["./montecarlo"], env=env, capture_output=True,
                            text=True, check=True).stdout
    return float(re.search(r"Tiempo: ([\d.]+)", salida).group(1))


tiempo(1, "static")  # calentamiento
med = {}
for plan in PLANES:
    for p in HILOS:
        ts = [tiempo(p, plan) for _ in range(REPS)]
        med[plan, p] = statistics.median(ts)
        print(f"{plan:16s} p={p:2d} mediana={med[plan, p]:.4f}  {ts}")

with open("schedule.md", "w") as fh:
    fh.write(f"Mediana de T(p) [s], N = 1e8, {REPS} corridas\n\n")
    fh.write("| p | " + " | ".join(PLANES) + " |\n|---|" + "---|" * len(PLANES) + "\n")
    for p in HILOS:
        fh.write(f"| {p} | " + " | ".join(f"{med[pl, p]:.4f}" for pl in PLANES) + " |\n")
print(open("schedule.md").read())
