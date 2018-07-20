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
    if (argc < 3) abort();
    int l = atoi(argv[1]),
        w = atoi(argv[2]);
    int ntot = 0;
    int n = w * l;
    int ti[4], tj[4];

    assert(w == l);
    if (argc == 3)
        for (auto yi = 0; yi < l; yi++)
            for (auto xi = 0; xi < w; xi++)
                for (auto yj = 0; yj < l; yj++)
                    for (auto xj = 0; xj < w; xj++) { // Forall a, i
                        ntot += 32;
                        nn_2d(ti, w, l, xi, yi);
                        nn_2d(tj, w, l, xj, yj);
                        for (auto i = 0; i < 4; i++)
                            for (auto j = 0; j < 4; j++) { // Forall b, j
                                printf("%d %d %d %d %d %d %d %d\n",
                                    idx_2d(w, xi, yi), 1, idx_2d(w, xj, yj), 1, ti[i], 0, tj[j], 0);
                                printf("%d %d %d %d %d %d %d %d\n",
                                    idx_2d(w, xi, yi), 0, idx_2d(w, xj, yj), 0, ti[i], 1, tj[j], 1);
                            }
                    }
    else {
        lattice physics(l);
        sc_corr quantity(physics, 8);
        for (auto yi = 0; yi < l; yi++) for (auto xi = 0; xi < l; xi++)
            for (auto yj = 0; yj < l; yj++) for (auto xj = 0; xj < l; xj++)
                for (auto ya = 0; ya < l; ya++) for (auto xa = 0; xa < l; xa++)
                    for (auto yb = 0; yb < l; yb++) for (auto xb = 0; xb < l; xb++)
                        for (auto si = 0; si < 2; si++) for (auto sj = 0; sj < 2; sj++)
                            for (auto sa = 0; sa < 2; sa++) for (auto sb = 0; sb < 2; sb++) {
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
