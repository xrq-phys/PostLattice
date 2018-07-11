#include <cstdio>
#include <cstdlib>
#include <vector>

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
    int n = w * l;
    int ti[4], tj[4];
    for (auto yi = 0; yi < l; yi++)
        for (auto xi = 0; xi < l; xi++)
            for (auto yj = 0; yj < l; yj++)
                for (auto xj = 0; xj < l; xj++) {
                    nn_2d(ti, w, l, xi, yi);
                    nn_2d(tj, w, l, xj, yj);
                    for (auto i = 0; i < 4; i++)
                        for (auto j = 0; j < 4; j++) {
                            printf("%d %d %d %d %d %d %d %d\n",
                                   idx_2d(w, xi, yi), 0, idx_2d(w, xj, yj), 1, ti[i], 0, tj[j], 1);
                            printf("%d %d %d %d %d %d %d %d\n",
                                   idx_2d(w, xi, yi), 1, idx_2d(w, xj, yj), 0, ti[i], 1, tj[j], 0);
                            printf("%d %d %d %d %d %d %d %d\n",
                                   idx_2d(w, xi, yi), 0, idx_2d(w, xj, yj), 1, ti[i], 1, tj[j], 0);
                            printf("%d %d %d %d %d %d %d %d\n",
                                   idx_2d(w, xi, yi), 1, idx_2d(w, xj, yj), 0, ti[i], 0, tj[j], 1);
                        }
                }

    return 0;
}