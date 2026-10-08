# Scripts for the WES brief communication (wes-2026-138)

Python studies, results and figures of the revision of the brief communication
on the radial induction of non-planar rotors (Wind Energy Science, wes-2026-138).
All runs use KSS-Blade V6 and the modified NREL 5 MW case (11 m/s, 12.1 rpm).

| file | content |
|---|---|
| `driver.py` | writes `Machine.in` / `BlaDes.in`, runs the binary, reads `Bem.out` and `VC.out` |
| `study3.py` | reference case, local distributions, grid study |
| `study4.py` | parametric study (45 cases) |
| `study5.py` | local loads, curved tip, drag sensitivity |
| `lw_check.py` | Limacher-Wood relation |
| `vcfull.py` | axial + radial velocity of semi-infinite vortex cylinders |
| `decomp.py` | power decomposition at fixed circulation |
| `figs2.py` | Fig. 1 and Fig. 2 |

`*.json` and `*.log` are the results of these runs, `*.png` the figures.

## How to run

The scripts expect one folder that holds the compiled binary `kss` together with
the polars (`*.aer`), `ProThick.in` and `ThickDis.in`. Its default location is
`../KSS-Blade-V6` relative to `scripts/`; any other folder can be given by the
environment variable `KSS_SRC`.

    mkdir KSS-Blade-V6
    cp ../SourceCode/*.f ../NREL-5MW/* KSS-Blade-V6/
    cd KSS-Blade-V6
    gfortran -std=legacy -fno-automatic -O1 -o kss mem.f KSS.f Sub1.f Sub2.f Sub3.f SubNum.f
    cd ../scripts
    python study3.py

Python needs `numpy` (and `matplotlib` for `figs2.py`).
