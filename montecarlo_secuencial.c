/* montecarlo_secuencial.c -- version secuencial de referencia (II.1).
   Igual al codigo del enunciado, con la semilla asignada y medicion de tiempo.
   Uso: ./montecarlo_secuencial <numero de muestras> */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[]) {
    long N = (argc > 1) ? atol(argv[1]) : 100000000L;
    unsigned int estado = 2223102u;   /* semilla del generador */
    long dentro = 0;

    double t0 = omp_get_wtime();
    for (long i = 0; i < N; i++) {
        double x = rand_r(&estado) / (double)RAND_MAX;
        double y = rand_r(&estado) / (double)RAND_MAX;
        if (x*x + y*y <= 1.0)
            dentro++;
    }
    double t = omp_get_wtime() - t0;

    printf("Estimacion de pi: %.6f\n", 4.0 * dentro / (double)N);
    printf("Tiempo: %.4f s\n", t);
    return 0;
}
