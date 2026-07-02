# GmQuant - Interest Rate Swap Curve Construction & Risk Engine

A single-file C++17 quantitative finance engine that builds discount curves from
market data, prices interest rate swaps against them, and computes **analytical**
(closed-form) risk sensitivities. It also ships with a small, extensible
object-oriented framework for adding new instruments and interpolation methods.

The whole solution lives in [`solution.cpp`](solution.cpp) with no external
dependencies it only uses the C++ standard library.

---

## What the program does

Given a set of quoted market instruments (cash deposits and par swaps), the engine:

1. **Builds discount curves** by *bootstrapping* discount factors (DFs) from the quotes.
2. **Interpolates** those curves using two methods: log-linear and an
   averaged-quadratic (Lagrange) scheme.
3. **Prices** a user-specified vanilla swap against every curve/method combination.
4. **Computes analytical risk**: the sensitivity of the swap's PV to each input
   market rate: using calculus (chain rule) rather than slow bump-and-revalue.

The work is organised into four parts (Q1–Q3), all driven from `main()`:

| Part  | What it produces |
| :---- | :--------------- |
| **Q1**   | Discount factor at a query maturity, for each curve × interpolation method |
| **Q2.1** | Present value **and** par rate of a new vanilla swap |
| **Q2.2** | Analytical risk vector: ∂PV/∂(each input rate) |
| **Q3**   | A generic OO framework (Strategy / Template Method / Factory) that the Q1-Q2 curve builders now use internally |

---

## Financial concepts implemented

### Curve bootstrapping
- **Cash curve**: each deposit's DF is solved directly:
  `DF(Tᵢ) = 1 / (1 + cᵢ · Tᵢ/360)`. Nodes are independent of one another.
- **Swap curve**: a par swap has PV = 0 at inception. Using the *telescoping
  identity* for the floating leg (`float PV = N·(1 − DF(T))`), each node reduces to
  `DF(Tₖ) = (1 − pₖ·Aₖ) / (1 + pₖ·δₖ)`, where `Aₖ` is the annuity of intermediate
  payments and `δₖ` the final-period day-count fraction.
- **Circularity for maturities > 1Y** (e.g. an 18-month node needs DFs that don't
  exist yet) is resolved with a **fixed-point iteration**: seed a guess via
  log-linear extrapolation, re-solve using the partially-built curve, and repeat
  until successive DFs differ by less than `1e-12`.

### Interpolation (on log-DF)
- **Linear**: straight line in log-DF space.
- **Averaged Quadratic (AQ)**: blends two overlapping Lagrange quadratics for a
  smoother curve without kinks, with fallbacks to linear / single-quadratic near
  the boundaries where fewer than 3–4 surrounding nodes exist.

### Pricing
- Fixed leg valued payment-by-payment with interpolated DFs; floating leg valued
  via the telescoping identity. Par rate = `(1 − DF(T)) / annuity`.

### Analytical risk (the hard part)
Instead of bumping each input rate and re-pricing, sensitivities are derived
in closed form:
- **Cash risk**: one term per node (`dDF/dc = −DF²·T/360`), no cascade.
- **Swap risk**: a change in one swap rate propagates through *all* later
  bootstrapped nodes, so the code performs a **forward-substitution cascade**
  combining the direct term with the interpolation sensitivities of every
  downstream node.

The full derivations (telescoping, bootstrap, cash/swap direct terms, and both
interpolation sensitivities) are written up in
[`Implementation_approach.docx`](Implementation_approach.docx).

---

## Design (Q3 framework)

The Q3 additions refactor the curve machinery behind clean OO abstractions so new
instruments and interpolators can be added **without touching existing code**
(Open/Closed Principle):

- **`Interpolator`** (Strategy): `LinearInterpolator`, `AQInterpolator`;
  a `CubicSplineInterpolator` can be dropped in by subclassing.
- **`CalibrationInstrument`** (Template Method): `CashDeposit`, `ParSwap`;
  an FRA/OIS/bond is added by implementing `maturity()`, `marketRate()`, `solveDF()`.
- **`CurveBootstrapper`** (Factory): calibrates a curve from any list of
  instruments with any interpolator, agnostic to concrete types.
- **`PricingInstrument`** / **`VanillaSwap`**: separates pricing from calibration.

---

## Input / Output format

The program reads **`input.csv`** from the working directory and writes
**`output.csv`**. (A sample result is included as [`Output.csv`](Output.csv).)

### `input.csv`
```
<numInstruments>
<maturity>,<cashRate%>,<swapRate%>     # one row per instrument
...
<queryMaturityInDays>                  # Q1 query point
<fixedRate%>,<maturity>,<fixedFreq>,<floatFreq>   # the swap to price (Q2)
```
Maturities accept unit suffixes: `D` (days), `W` (weeks ×7), `M` (months ×30),
`Y` (years ×360). Frequencies are like `6m`. Rates are entered in percent.

### `output.csv`
Each row has 4 columns: one per curve/method combination
(*cash-linear, cash-aq, swap-linear, swap-aq*):
1. Q1: discount factor at the query maturity
2. Q2.1: swap PV
3. Q2.1: par swap rate
4. …onward : Q2.2 analytical risk, one row per input instrument

---

## Build & run

Requires a C++17 compiler.

```bash
# Compile
g++ -std=c++17 -O2 -o gmquant solution.cpp

# Provide an input.csv in the same directory, then run
./gmquant
# → writes output.csv
```

---

## Files

| File | Description |
| :--- | :--- |
| `solution.cpp` | The complete engine (parsing, bootstrapping, pricing, risk, Q3 framework) |
| `Implementation_approach.docx` | Write-up: approach, full mathematical derivations, assumptions, complexity, design notes |
| `Output.csv` | Sample output for reference |
| `README.md` | This file |
