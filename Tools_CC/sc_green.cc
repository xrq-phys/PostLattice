#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <vector>
#include "sc_classes.hh"

using namespace std;

inline int idx_2d(int w, int x, int y)
{ return w * y + x; }

void r_2d(int *r, int w, int i)
{
    r[0] = i % w;
    r[1] = i / w;
}

void nn_2d(int *nn, int w, int l, int xi, int yi)
{
    nn[0] = idx_2d(w, xi, (yi + 1)     % l);
    nn[1] = idx_2d(w, xi, (yi - 1 + l) % l);
    nn[2] = idx_2d(w, (xi + 1)     % w, yi);
    nn[3] = idx_2d(w, (xi - 1 + w) % w, yi);
}

int main(const int argc, const char *argv[])
{
    assert(argc > 2);
    int l = atoi(argv[1]),
        w = atoi(argv[2]);
    int ntot = 0;
    int n = w * l;
    int ti[4], tj[4], ta[4];

    assert(w == l);
    if (argc == 3)
        for (int yi = 0; yi < l; yi++)
            for (int xi = 0; xi < w; xi++)
                for (int yj = 0; yj < l; yj++)
                    for (int xj = 0; xj < w; xj++) { // Forall a, i
                        ntot += 32;
                        nn_2d(ti, w, l, xi, yi);
                        nn_2d(tj, w, l, xj, yj);
                        for (int i = 0; i < 4; i++)
                            for (int j = 0; j < 4; j++) { // Forall b, j
                                printf("%d %d %d %d %d %d %d %d\n",
                                    idx_2d(w, xi, yi), 1, idx_2d(w, xj, yj), 1, ti[i], 0, tj[j], 0);
                                printf("%d %d %d %d %d %d %d %d\n",
                                    idx_2d(w, xi, yi), 0, idx_2d(w, xj, yj), 0, ti[i], 1, tj[j], 1);
                            }
                    }
    else if (argv[3][0] == '-') {
        int ya = 0, xa = 0;
        for (int yj = 0; yj < l; yj++)
            for (int xj = 0; xj < w; xj++) {
                ntot += 32;
                nn_2d(ta, w, l, xa, ya);
                nn_2d(tj, w, l, xj, yj);
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++) {
                        printf("%d %d %d %d %d %d %d %d\n", ta[i], 1, idx_2d(w, xj, yj), 1, idx_2d(w, xa, ya), 0, tj[j], 0);
                        printf("%d %d %d %d %d %d %d %d\n", ta[i], 0, idx_2d(w, xj, yj), 0, idx_2d(w, xa, ya), 1, tj[j], 1);
                    }
            }
    } else if (argv[3][0] == '+') {
        lattice physics(l);
        sc_corr quantity(physics, 8);
        for (int yi = 0; yi < l; yi++) for (int xi = 0; xi < l; xi++)
            for (int yj = 0; yj < l; yj++) for (int xj = 0; xj < l; xj++)
                for (int ya = 0; ya < l; ya++) for (int xa = 0; xa < l; xa++)
                    for (int yb = 0; yb < l; yb++) for (int xb = 0; xb < l; xb++)
                        for (int si = 0; si < 2; si++) for (int sj = 0; sj < 2; sj++)
                            for (int sa = 0; sa < 2; sa++) for (int sb = 0; sb < 2; sb++) {
                                if (quantity.validate(physics.idx_2d(xi, yi), si, 
                                                      physics.idx_2d(xj, yj), sj,
                                                      physics.idx_2d(xa, ya), sa,
                                                      physics.idx_2d(xb, yb), sb)) {
                                    printf("%d %d %d %d %d %d %d %d\n", physics.idx_2d(xi, yi), si, 
                                                                        physics.idx_2d(xj, yj), sj,
                                                                        physics.idx_2d(xa, ya), sa,
                                                                        physics.idx_2d(xb, yb), sb);
                                    ntot++;
                                }
                            }
    }
    printf(" NCisAjsCktAltSC %d ", ntot);

    return 0;
}
