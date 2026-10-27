/* montecarlo.c -- version paralela con OpenMP.
   Uso: ./montecarlo <numero de muestras>
   Hilos: OMP_NUM_THREADS=p.  Planificacion: OMP_SCHEDULE=static | "dynamic,1048576".

   Generador: LCG de 64 bits propio (no rand_r) para poder saltar k valores
   en O(log k). La muestra i usa siempre los valores 2i y 2i+1 de una unica
   secuencia que parte de la semilla, sin importar que hilo la procese: cada
   hilo salta al inicio de su bloque de iteraciones (>= N/p muestras con
   schedule static), asi los flujos no se solapan y el resultado es identico
   para cualquier p y cualquier planificacion. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>

#define SEMILLA 2223102u
#define LCG_A 6364136223846793005ULL   /* constantes de Knuth (MMIX) */
#define LCG_C 1442695040888963407ULL

/* Siguiente valor uniforme en [0,1): usa los 53 bits altos del estado. */
static inline double aleatorio(uint64_t *estado) {
    *estado = LCG_A * *estado + LCG_C;
    return (*estado >> 11) * 0x1.0p-53;
}

/* Estado tras avanzar k pasos desde s (Brown, 1994): compone la
   transformacion x -> a*x + c consigo misma por cuadrados sucesivos. */
static uint64_t saltar(uint64_t s, uint64_t k) {
    uint64_t a = LCG_A, c = LCG_C, a_tot = 1, c_tot = 0;
    while (k) {
        if (k & 1) {
            a_tot *= a;
            c_tot = c_tot * a + c;
        }
        c *= a + 1;
        a *= a;
        k >>= 1;
    }
    return a_tot * s + c_tot;
}

int main(int argc, char *argv[]) {
    long N = (argc > 1) ? atol(argv[1]) : 100000000L;
    long dentro = 0;
    uint64_t estado = SEMILLA;  /* copia por hilo: flujo del hilo */
    long siguiente = -1;        /* proxima i esperada por el hilo; -1 obliga a saltar */

    /* Sin OMP_SCHEDULE, libgomp usa dynamic,1: un salto por muestra. */
    if (!getenv("OMP_SCHEDULE"))
        omp_set_schedule(omp_sched_static, 0);

    double t0 = omp_get_wtime();
    #pragma omp parallel for reduction(+:dentro) firstprivate(estado, siguiente) schedule(runtime)
    for (long i = 0; i < N; i++) {
        if (i != siguiente)   /* inicio de un bloque nuevo */
            estado = saltar(SEMILLA, 2 * (uint64_t)i);
        double x = aleatorio(&estado);
        double y = aleatorio(&estado);
        if (x*x + y*y <= 1.0)
            dentro++;
        siguiente = i + 1;
    }
    double t = omp_get_wtime() - t0;

    printf("Estimacion de pi: %.6f\n", 4.0 * dentro / (double)N);
    printf("Aciertos: %ld de %ld\n", dentro, N);
    printf("Hilos: %d  Tiempo: %.4f s\n", omp_get_max_threads(), t);
    return 0;
}
