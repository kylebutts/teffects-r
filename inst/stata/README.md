# Stata reference results

`run_teffects_examples.do` runs the Stata examples used by the R test suite and
writes their target estimates and standard errors to
`teffects_reference.csv`. 

Run it with `inst/stata` as the working directory:

```sh
cd inst/stata
stata-mp -b do run_teffects_examples.do
```

The script executes all 40 displayed `teffects` commands extracted from the
dedicated manual entries. 
It also runs 29 compact reference cases used to exercise target statistics and matching-inference branches directly. 
Exact repetitions are retained because they belong to distinct documented examples.

The manual census uses `cattaneo2.dta`, `cattaneo3.dta`, `bdsianesi5.dta`, and `pollution.dta` from the `data` directory. Every manual command posts every coefficient and standard error in `e(b)` and `e(V)`.

The committed `teffects_reference.csv` is provisionally populated from the numbers printed in the manual. 
It contains 169 rows total: the 35 commands with displayed results plus the existing default reference cases. 
Blank standard errors occur for four coefficients shown only in a `coeflegend` table. 
A real Stata run should replace these printed-precision values with full-precision results.

Return `teffects_reference.csv`. 
If the batch does not finish, also return `run_teffects_examples.log` so the failing manual case can be identified.
