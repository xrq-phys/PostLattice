PostLattice
===========

A simple helper program for many-body calculations.

Compiling
---------

Use CMake, but **remember to recursively clone all submodules**.

Usage
-----

Master branch of this program supports 3 modes: `Gen` (Generate Green's function entries required for calculating spin and uperconducting correlations), `Measure` (Measure physical quantities from Green's function output of mVMC/HPhi) and `Plot` (To be described in another branch).

### Input

PostLattice is called in this way:

```
Post post.ini [G1F] [G2F]
```
Here `Post` is the executable file, `post.ini` serves as input configuration and `G1F G2F` are optional parameters that overrides some entries in `post.ini`.

A full `post.ini` file should look like this:

```ini
[Control]
GreenOne = Green1.def
GreenTwo = Green2.def
GreenOneOut = CisAjs.dat
GreenTwoOut = CisAjsCktAltDC.dat
Mode     = Gen
Verbose  = 0
Green_Verbose = 0
Green_Legacy = 0

[Physics]
System = Square
W      = 10
L      = 10

[Operator]
Spin_Structure = 1
AF_N_X         = 10
AF_N_Y         = 10
AF_Out         = sstruct.4119.txt
SC             = d
SC_NumR        = 90
SC_Out         = sc.d.4119.txt
SC_Simple      = 1
SC_Stat        = s
```

As is seen from the example, configuration file consists of 3 sections:

- Control section: controls file and other IO behaviours;
 - `Mode`: Sets calculation mode described above.
 - `GreenOne` and `GreenTwo`: 1/2-body Green's function index files. Output filename when `Mode=Gen`. Input if `Mode=Measure` and `Green_Verbose=0`. Not used otherwise. See option `Green_Legacy` for more details.
 - `GreenOneOut` and `GreenTwoOut`: Measured value of 1/2-body Green's functions. See `Green_Verbose` for the files' format.
 - `Green_Verbose`: Specifies format of Green's function components calculated by mVMC/HPhi. If you are using legacy version of mVMC or TNVMC which writes GF values consecutively in one line without repeating indicies, set this value to `0` and `1` otherwise.
 - `Green_Legacy`: As when `Green_Verbose=0`, PostLattice needs re-reading the supplied `GreenOne` and `GreenTwo` (`GreenOne/Two` is not used otherwise when measuring because new mVMC reprints site indices when writing output), please switch on this option if the indices file you're using is for legacy version of mVMC or TNVMC. But be aware that GF indices generated by PostLattice is always for new version.
 - `Verbose`: Adds some debug output.
- `Physics`: Specifies lattice system.
 - `System`: Choose between `Square`, `Triangular` and `Honeycomb`, where `Honeycomb` contains 2 sites per unit cell.
 - `a0W` or `W`: Width of the lattice. Note that currently only horizontal $a_0$ is supported.
 - `a1L` and `a1W`: secondary lattice vector. See mVMC document for detailed explanation.
- `Measure`: Sets parameter for physical quantity measurements.
 - `Spin_Structure`: Whether spin structure will be calculated.
 - `AF_N_X` and `AF_N_Y`: Number of Fourier transformation points in $x$ and $y$ direction.
 - `AF_Out`: Output file for spin structure factor measurement.
 - `SC`: Waveform for superconducting correlation. This can be any of the supported waveforms: `s`, `1`, `2` (common), `3`, `f` (for triangular) or `d` (for square). Custom waveform is also supported: We can assign a sign to each of the neighbour site like `SC=+-+-` (equivalent to `SC=d`). Set to `-` if you do not want to measure superconducting correlation.
 - `SC_NumR`: Number of SC to measure. Specify a large enough number and the program will automatically cut it off.
 - `SC_Simple`: Whether to skip one level of averaging the SC correlation. Please set to `1` in all cases.
 - `SC_Stat`: In measuring SC correlation, there are always different points that have the same distance to the origin. Set `s` to average them and `m` to pick maximum value among them. (`m` recommended)
 
### Output
 
Spin structure output format is like:

```
0.000 0.000 1.000 0.000
0.000 1.040 0.000 0.000
...
0.000 5.413 0.000 0.000

1.040 0.000 0.000 0.000
...
```
The first 2 columns are $k$-point coordinates and the last 2 are real and imaginary parts of the measured spin structure.

SC output format is like:

```
0.000 2.000 1
1.000 1.400 4
...
```
The first column is distance and column 2 / column 3 is the SC correlation at this distance.
 




