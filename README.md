# Two Floors: Housing Benefit, the Public Pension and the Shape of Income among Older People Living Alone

Samir Orujov, School of Business, ADA University, Baku

The disposable income of Danes aged 65 and over who live alone is bimodal. In register data from the Luxembourg Income Study for 2004 to 2019, a test of unimodality for binned incomes rejects for this group from 2013 onwards, and for those who live on public transfers it rejects in every wave. The two modes are the public pension without and with the pensioners' housing benefit: the gap between the modes is close to the median benefit, the recipients are renters and most non-recipients are owners, and removing the benefit from income removes the second mode. With LIS waves pooled, the same structure appears in the Netherlands, the United Kingdom and France; in Norway, Sweden, Germany, the United States and Ireland no second mode is due to the housing allowance. A statutory floor shapes the income distribution when a large share of a group lives on it, and a second transfer that only some of them receive splits it in two.

JEL: D31, H55, I38, R21. Keywords: income distribution, bimodality, minimum pension, housing benefit, one-person households, Luxembourg Income Study.

Working paper: [SSRN 6445278](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6445278). The test of unimodality used here is developed in a companion paper, Orujov (2026), *Identification and inference for the shape of a distribution from rounded reports*, [SSRN 6454722](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6454722); its code is at [sorujov/rounded-reports](https://github.com/sorujov/rounded-reports).

## Contents

| Folder | What it holds |
|---|---|
| `manuscript/` | LaTeX source of the paper (`main.tex`; `main_anon.tex` is the anonymised version) and the script for Figure 1 |
| `lis/` | The R jobs submitted to LISSY, the LIS remote-execution system |
| `lis/output/` | The aggregate output returned by LISSY for each job |
| `simulation/` | Size and power of the cell test, simulated outside LIS |

## Data and access

LIS microdata cannot be downloaded. The jobs in `lis/` run on LISSY (<https://www.lisdatacenter.org>), which any registered user can access; they read the LIS household and person files through the `lissyrtools` package and print only test statistics, shares and other aggregates. The files in `lis/output/` are those aggregates as returned. Reproducing them requires a LIS account.

## The LIS jobs

| Job | Content |
|---|---|
| `job6_dk.R` | Denmark: decomposition by income floor, one-person households |
| `job7_dk.R` | Denmark: decomposition by pension source |
| `job8_no.R` | Norway: wave by wave decomposition |
| `job9_dk.R` | Denmark: the two modes among transfer-reliant older people |
| `job10_sefi.R` | Sweden and Finland: wave by wave decomposition |
| `job11b_dk.R` | Denmark: income without the housing benefit (counterfactual) |
| `job12_dk.R` | Denmark: other definitions of transfer reliance, benefit maximum |
| `job13_tenure.R` | Denmark and Norway: housing tenure of transfer-reliant older people |
| `job17_fig.R` | Denmark: binned distribution for Figure 1 (cells of 5,000 kroner; cells with fewer than 30 observations suppressed) |
| `job18_light.R` | Denmark and Norway: power, grid shift, reweighting, low-income lines |
| `job20_dk.R` | Denmark: refined transfer-reliant group, poverty lines, cells of 10,000 kroner |
| `job22_other.R` | Eight further countries: the two-floor test wave by wave |
| `job23_base.R`, `job23a_pool.R`, `job23b_pool.R`, `job23c_pool.R` | Pooled waves, incomes scaled to the pension-only median (cross-country section; Denmark as control) |
| `job23d_c2.R` | Denmark 2019: the test when the share of benefit recipients is reduced |
| `job24_ie_uk.R` | Ireland and the United Kingdom: which transfer makes the second mode |
| `job25_checks.R` | Wave by wave checks, alternative scale, early French waves, parts of the refined Danish group |

## The simulation

`simulation/cell_test.R` is the base-R implementation of the cell test run on LISSY (it is the opening block of every job). `size_sim.R` draws from four unimodal distributions (normal, lognormal, flat-topped, sharp peak with a long right tail) and reports the rejection rate; `size_sim2.R` stores the p-values for the flat-topped case. The files in `simulation/output/` are the results reported in Section 3 and the appendix. Usage:

```
cd simulation
Rscript size_sim.R flattop 7000 100
Rscript size_sim2.R 20000 100 1
```

## Building the paper

```
cd manuscript
python3 make_fig.py          # Figure 1 from ../lis/output/job17_fig.txt
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

Requirements: R 4.2 or later (base R only for the simulation; `lissyrtools` and `purrr` on LISSY), Python 3 with NumPy and Matplotlib, and a LaTeX distribution.

## Acknowledgement

The empirical results use the Luxembourg Income Study Database. I am grateful to the LIS Cross-National Data Center for access.

## Citation

Orujov, S. (2026). Two Floors: Housing Benefit, the Public Pension and the Shape of Income among Older People Living Alone. Working paper, ADA University. SSRN 6445278.

## Licence

Code: MIT licence (see `LICENSE`). The text of the paper remains the author's; please cite it rather than reuse it.
