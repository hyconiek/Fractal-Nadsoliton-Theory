/* Exact dense modular row rank for a raw row-major uint32 matrix.
 *
 * Usage: fin_modrank_u32 matrix.bin rows cols prime [echelon.bin]
 * The implementation performs row-echelon elimination over GF(prime).
 * It is intentionally dependency-free and parallelizes independent row
 * eliminations with OpenMP when compiled with -fopenmp.
 */

#include <errno.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static uint32_t mod_pow(uint32_t a, uint32_t e, uint32_t p) {
    uint64_t x = a, out = 1;
    while (e) {
        if (e & 1U) out = (out * x) % p;
        x = (x * x) % p;
        e >>= 1U;
    }
    return (uint32_t)out;
}

static void swap_rows(uint32_t *a, size_t n, size_t i, size_t j) {
    if (i == j) return;
    for (size_t k = 0; k < n; ++k) {
        uint32_t t = a[i*n+k];
        a[i*n+k] = a[j*n+k];
        a[j*n+k] = t;
    }
}

int main(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        fprintf(stderr, "usage: %s matrix.bin rows cols prime [echelon.bin]\n", argv[0]);
        return 2;
    }
    const char *path = argv[1];
    size_t m = (size_t)strtoull(argv[2], NULL, 10);
    size_t n = (size_t)strtoull(argv[3], NULL, 10);
    uint32_t p = (uint32_t)strtoul(argv[4], NULL, 10);
    if (!m || !n || p < 3) return 3;

    size_t count = m*n;
    uint32_t *a = (uint32_t *)malloc(count*sizeof(uint32_t));
    if (!a) {
        fprintf(stderr, "allocation failed for %zu entries\n", count);
        return 4;
    }
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "open failed: %s\n", path);
        free(a);
        return 5;
    }
    size_t got = fread(a, sizeof(uint32_t), count, f);
    fclose(f);
    if (got != count) {
        fprintf(stderr, "short read: got %zu expected %zu\n", got, count);
        free(a);
        return 6;
    }

    size_t *pivot_cols = (size_t *)malloc(m*sizeof(size_t));
    if (!pivot_cols) { free(a); return 7; }
    size_t rank = 0;
    for (size_t col = 0; col < n && rank < m; ++col) {
        size_t pivot = rank;
        while (pivot < m && a[pivot*n+col] == 0U) ++pivot;
        if (pivot == m) continue;
        swap_rows(a, n, rank, pivot);

        uint32_t inv = mod_pow(a[rank*n+col], p-2U, p);
        for (size_t j = col; j < n; ++j)
            a[rank*n+j] = (uint32_t)(((uint64_t)a[rank*n+j]*inv) % p);

        #pragma omp parallel for schedule(static)
        for (ptrdiff_t ii = (ptrdiff_t)rank+1; ii < (ptrdiff_t)m; ++ii) {
            size_t i = (size_t)ii;
            uint32_t factor = a[i*n+col];
            if (!factor) continue;
            a[i*n+col] = 0U;
            for (size_t j = col+1; j < n; ++j) {
                uint32_t sub = (uint32_t)(((uint64_t)factor*a[rank*n+j]) % p);
                uint32_t value = a[i*n+j];
                a[i*n+j] = value >= sub ? value-sub : value+p-sub;
            }
        }
        pivot_cols[rank] = col;
        ++rank;
        if (rank % 100 == 0 || rank == m)
            fprintf(stderr, "rank_progress=%zu col=%zu\n", rank, col);
    }
    printf("rank=%zu\n", rank);
    printf("pivots=");
    for (size_t i = 0; i < rank; ++i)
        printf("%s%zu", i ? "," : "", pivot_cols[i]);
    printf("\n");
    if (argc == 6) {
        FILE *out = fopen(argv[5], "wb");
        if (!out || fwrite(a, sizeof(uint32_t), count, out) != count) {
            fprintf(stderr, "failed to write echelon output\n");
            if (out) fclose(out);
            free(pivot_cols); free(a); return 8;
        }
        fclose(out);
    }
    free(pivot_cols);
    free(a);
    return 0;
}
