#include "operators.hh"

void operators::sc_corr::measure(int a, int sa, int b, int sb, 
                                 int i, int si, int j, int sj, double x)
{
    int sign;
    int rmin;
    int rb[2], 
        ri[2];
    
    if (validate(a, sa, b, sb, i, si, j, sj)) {
        sign = 1;

        if (sb != si) 
            sign = -sign;

        // Waveform factor only affects this part.
        if (form == 's') { /* Do nothing. */ }
        else if (form == 'd')
            if (system.nn[a][b] != system.nn[i][j] && 
                system.nn[a][b] +  system.nn[i][j] == 1)
                sign = -sign;

        rmin = system.calc_rmin(b, i);
        for (int ii = 0; ii < rc_count; ii++)
            if (rmin == system.r_c[ii])
#pragma omp critical
            values[ii] += x * sign / (2. * system.n);
    }
}