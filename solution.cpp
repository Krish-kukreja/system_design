// ============================================================
//  - Strategy    : Interpolator hierarchy (swap method at runtime)
//  - Template Method : CalibrationInstrument defines the bootstrap contract
//  - Factory     : CurveBootstrapper creates calibrated curves generically
//  - Open/Closed : New instruments and interpolators added without
//                  modifying any existing code
// ============================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <memory>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <utility>

// ────────────────────────────────────────────────────────────
//  UTILITY FUNCTIONS
// ────────────────────────────────────────────────────────────

// Trim whitespace for ribustness
static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    std::string result = s.substr(a, b - a + 1);
    if (result.size() >= 3 &&
        (unsigned char)result[0] == 0xEF &&
        (unsigned char)result[1] == 0xBB &&
        (unsigned char)result[2] == 0xBF)
        result = result.substr(3);
    return result;
}

// Day count fraction: (t2 - t1) / 360
static double dcf(double t2, double t1) { return (t2 - t1) / 360.0; }

// Parse maturity string → days  ("nD"→n, "nW"→n*7, "nM"→n*30, "nY"→n*360)
static double parseMaturity(const std::string& raw) {
    std::string s = trim(raw);
    if (s.empty()) throw std::invalid_argument("Empty maturity string");
    char   unit = s.back();
    double v    = std::stod(s.substr(0, s.size() - 1));
    switch (unit) {
        case 'D': case 'd': return v;
        case 'W': case 'w': return v * 7.0;
        case 'M': case 'm': return v * 30.0;
        case 'Y': case 'y': return v * 360.0;
        default: throw std::invalid_argument("Unknown maturity unit in: " + s);
    }
}

// Parse frequency string → days  ("1m"→30, "3m"→90, "6m"→180, "12m"→360)
static double parseFrequency(const std::string& raw) {
    std::string s = trim(raw);
    if (!s.empty() && (s.back() == 'm' || s.back() == 'M')) s.pop_back();
    return std::stod(s) * 30.0;
}

// Split CSV line into trimmed tokens
static std::vector<std::string> splitCSV(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ','))
        tokens.push_back(trim(tok));
    return tokens;
}

// ────────────────────────────────────────────────────────────
//  MARKET DATA — raw input row (maturity + two quoted rates)
// ────────────────────────────────────────────────────────────

struct Instrument {
    double maturity;   // days
    double cashRate;   // decimal
    double swapRate;   // decimal
};

// ────────────────────────────────────────────────────────────
//  LAGRANGE HELPERS
// ────────────────────────────────────────────────────────────

// Evaluates unique quadratic through (t0,y0),(t1,y1),(t2,y2) at t
static double lagrangeQuad(double t,
                            double t0, double y0,
                            double t1, double y1,
                            double t2, double y2) {
    double l0 = (t-t1)*(t-t2)/((t0-t1)*(t0-t2));
    double l1 = (t-t0)*(t-t2)/((t1-t0)*(t1-t2));
    double l2 = (t-t0)*(t-t1)/((t2-t0)*(t2-t1));
    return y0*l0 + y1*l1 + y2*l2;
}

// Returns {L0,L1,L2} — Lagrange basis polynomials at t for nodes (t0,t1,t2)
static std::array<double,3> lagrangeBases(double t,
                                           double t0, double t1, double t2) {
    double l0 = (t-t1)*(t-t2)/((t0-t1)*(t0-t2));
    double l1 = (t-t0)*(t-t2)/((t1-t0)*(t1-t2));
    double l2 = (t-t0)*(t-t1)/((t2-t0)*(t2-t1));
    return {l0, l1, l2};
}

// ────────────────────────────────────────────────────────────
//  CURVE
//  Stores (time, DF) node pairs. Interpolation methods on
//  log(DF): linear and averaged-quadratic (AQ).
// ────────────────────────────────────────────────────────────

class Curve {
public:
    std::vector<double> times;
    std::vector<double> dfs;

    void addNode(double t, double df) { times.push_back(t); dfs.push_back(df); }
    int  size() const { return static_cast<int>(times.size()); }

    double linearInterp(double t) const {
        int n = size();
        if (n == 0) throw std::runtime_error("Curve is empty");
        if (t <= times[0])     return dfs[0];
        if (t >= times[n-1])   return dfs[n-1];
        int i = static_cast<int>(
            std::upper_bound(times.begin(), times.end(), t) - times.begin()) - 1;
        double w = (t - times[i]) / (times[i+1] - times[i]);
        return std::exp(std::log(dfs[i]) + w*(std::log(dfs[i+1])-std::log(dfs[i])));
    }

    double aqInterp(double t) const {
        int n = size();
        if (n == 0) throw std::runtime_error("Curve is empty");
        if (t <= times[0])     return dfs[0];
        if (t >= times[n-1])   return dfs[n-1];
        if (n < 3)             return linearInterp(t);
        int i = static_cast<int>(
            std::upper_bound(times.begin(), times.end(), t) - times.begin()) - 1;
        if (i == 0)            return linearInterp(t);
        double wR = (t - times[i])     / (times[i+1] - times[i]);
        double wL = (times[i+1] - t)   / (times[i+1] - times[i]);
        double qi = lagrangeQuad(t, times[i-1], std::log(dfs[i-1]),
                                    times[i],   std::log(dfs[i]),
                                    times[i+1], std::log(dfs[i+1]));
        if (i + 2 >= n) return std::exp(qi);
        double qi1 = lagrangeQuad(t, times[i],   std::log(dfs[i]),
                                     times[i+1], std::log(dfs[i+1]),
                                     times[i+2], std::log(dfs[i+2]));
        return std::exp(wL*qi + wR*qi1);
    }

    double interpolate(double t, const std::string& method) const {
        if (method == "linear") return linearInterp(t);
        if (method == "aq")     return aqInterp(t);
        throw std::invalid_argument("Unknown interpolation method: " + method);
    }
};

// ────────────────────────────────────────────────────────────
//  FORWARD DECLARATIONS
//  Required so Interpolator and PricingInstrument can reference
//  sensitivity / pricing functions defined later in the file.
// ────────────────────────────────────────────────────────────

static std::vector<double> computeInterpSensitivities(
    double t, const Curve& curve, const std::string& method);

static std::vector<double> computePVNodeSensitivities(
    const Curve& curve, const std::string& method,
    double fixedRate, double maturityDays, double fixedFreqDays,
    double notional = 100.0);

static std::vector<double> computeCashRisk(
    const Curve& cashCurve,
    const std::vector<Instrument>& insts,
    const std::vector<double>& pvNodeSens);

static std::vector<double> computeSwapRisk(
    const Curve& swapCurve,
    const std::vector<Instrument>& insts,
    const std::string& method,
    const std::vector<double>& pvNodeSens);

// ════════════════════════════════════════════════════════════
//  Q3 — GENERIC CURVE CONSTRUCTION FRAMEWORK
// ════════════════════════════════════════════════════════════

// ────────────────────────────────────────────────────────────
//  ABSTRACT INTERPOLATOR  (Strategy Pattern)
//
//  Encapsulates an interpolation algorithm. Any new method
//  (e.g. cubic spline) is added by:
//    1. Subclassing Interpolator
//    2. Implementing getDF() and sensitivity()
//  No changes to bootstrapper, pricing, or risk code.
// ────────────────────────────────────────────────────────────

class Interpolator {
public:
    virtual ~Interpolator() = default;

    // Discount factor at t, interpolated from curve nodes
    virtual double getDF(double t, const Curve& curve) const = 0;

    // ∂DF(t)/∂DF(T_j) for every node j — used in analytical risk
    virtual std::vector<double> sensitivity(double t, const Curve& curve) const = 0;

    // String tag — used to bridge into existing string-based dispatch
    virtual std::string name() const = 0;
};

// Concrete: linear interpolation on log(DF)
class LinearInterpolator : public Interpolator {
public:
    double getDF(double t, const Curve& curve) const override {
        return curve.linearInterp(t);
    }
    std::vector<double> sensitivity(double t, const Curve& curve) const override {
        return computeInterpSensitivities(t, curve, "linear");
    }
    std::string name() const override { return "linear"; }
};

// Concrete: averaged-quadratic interpolation on log(DF)
class AQInterpolator : public Interpolator {
public:
    double getDF(double t, const Curve& curve) const override {
        return curve.aqInterp(t);
    }
    std::vector<double> sensitivity(double t, const Curve& curve) const override {
        return computeInterpSensitivities(t, curve, "aq");
    }
    std::string name() const override { return "aq"; }
};

// ── Extension point: cubic spline would be added here ────────
//
//  class CubicSplineInterpolator : public Interpolator {
//  public:
//      double getDF(double t, const Curve& curve) const override {
//          // fit natural cubic spline on log(DF) and evaluate at t
//      }
//      std::vector<double> sensitivity(double t, const Curve& curve) const override {
//          // cubic spline sensitivity — spreads across more nodes than linear
//      }
//      std::string name() const override { return "cubic_spline"; }
//  };
//
//  Usage: CurveBootstrapper::calibrate(instruments, CubicSplineInterpolator{})
//  Zero changes to any other part of the code.

// ────────────────────────────────────────────────────────────
//  BOOTSTRAP SWAP NODE (implementation detail used by ParSwap)
// ────────────────────────────────────────────────────────────

static double bootstrapSwapNode(double swapRate, double T,
                                 Curve& curve, const std::string& method) {
    if (T <= 180.0)
        return 1.0 / (1.0 + swapRate * T / 360.0);

    std::vector<double> dates;
    for (double t = 180.0; t <= T + 1e-9; t += 180.0)
        dates.push_back(std::round(t));
    dates.back() = T;
    int n = static_cast<int>(dates.size());

    double guess;
    {
        int sz = curve.size();
        if (sz >= 2) {
            double t1 = curve.times[sz-2], t2 = curve.times[sz-1];
            double y1 = std::log(curve.dfs[sz-2]), y2 = std::log(curve.dfs[sz-1]);
            guess = std::exp(y2 + (y2-y1)/(t2-t1) * (T-t2));
        } else {
            guess = 1.0 / (1.0 + swapRate * T / 360.0);
        }
        guess = std::max(1e-12, std::min(1.0, guess));
    }

    for (int iter = 0; iter < 100; iter++) {
        curve.addNode(T, guess);
        double annuity = 0.0;
        for (int j = 0; j < n - 1; j++) {
            double t_prev = (j == 0) ? 0.0 : dates[j-1];
            annuity += curve.interpolate(dates[j], method) * dcf(dates[j], t_prev);
        }
        curve.times.pop_back();
        curve.dfs.pop_back();
        double new_df = (1.0 - swapRate * annuity) /
                        (1.0 + swapRate * dcf(T, dates[n-2]));
        new_df = std::max(1e-12, std::min(1.0, new_df));
        if (std::abs(new_df - guess) < 1e-12) return new_df;
        guess = new_df;
    }
    return guess;
}

// ────────────────────────────────────────────────────────────
//  ABSTRACT CALIBRATION INSTRUMENT  (Template Method Pattern)
//
//  Defines the contract for any instrument that can calibrate
//  a discount curve node. The CurveBootstrapper calls solveDF()
//  without knowing which instrument type it is.
//
//  New instrument types (FRA, OIS, bond) are added by:
//    1. Subclassing CalibrationInstrument
//    2. Implementing maturity(), marketRate(), solveDF()
// ────────────────────────────────────────────────────────────

class CalibrationInstrument {
public:
    virtual ~CalibrationInstrument() = default;

    // Maturity of this instrument in days
    virtual double maturity() const = 0;

    // Quoted market rate for this instrument
    virtual double marketRate() const = 0;

    // Returns the discount factor at maturity such that PV(instrument) = 0.
    // The partially-built curve and interpolator are available for
    // instruments that need earlier nodes (e.g. swaps).
    virtual double solveDF(Curve& curve, const Interpolator& interp) const = 0;
};

// Concrete: cash deposit — DF solved directly, no curve needed
class CashDeposit : public CalibrationInstrument {
    double mat_, rate_;
public:
    CashDeposit(double maturity, double rate) : mat_(maturity), rate_(rate) {}
    double maturity()    const override { return mat_; }
    double marketRate()  const override { return rate_; }
    double solveDF(Curve& /*curve*/, const Interpolator& /*interp*/) const override {
        return 1.0 / (1.0 + rate_ * mat_ / 360.0);
    }
};

// Concrete: par interest rate swap — DF solved via bootstrap
class ParSwap : public CalibrationInstrument {
    double mat_, rate_;
public:
    ParSwap(double maturity, double rate) : mat_(maturity), rate_(rate) {}
    double maturity()   const override { return mat_; }
    double marketRate() const override { return rate_; }
    double solveDF(Curve& curve, const Interpolator& interp) const override {
        return bootstrapSwapNode(rate_, mat_, curve, interp.name());
    }
};

// ── Extension point: FRA would be added here ─────────────────
//
//  class ForwardRateAgreement : public CalibrationInstrument {
//  public:
//      ForwardRateAgreement(double startDate, double endDate, double fraRate)
//          : start_(startDate), end_(endDate), rate_(fraRate) {}
//      double maturity()   const override { return end_; }
//      double marketRate() const override { return rate_; }
//      double solveDF(Curve& curve, const Interpolator& interp) const override {
//          // PV(FRA) = 0 → DF(end) = DF(start) / (1 + rate * DCF(end, start))
//          double dfStart = interp.getDF(start_, curve);
//          return dfStart / (1.0 + rate_ * dcf(end_, start_));
//      }
//  private:
//      double start_, end_, rate_;
//  };

// ────────────────────────────────────────────────────────────
//  CURVE BOOTSTRAPPER  (Factory Method Pattern)
//
//  Calibrates a discount curve from any vector of
//  CalibrationInstruments using any Interpolator.
//  Completely agnostic to instrument type and interpolation.
// ────────────────────────────────────────────────────────────

class CurveBootstrapper {
public:
    static Curve calibrate(
        const std::vector<std::unique_ptr<CalibrationInstrument>>& instruments,
        const Interpolator& interp)
    {
        Curve curve;
        for (const auto& inst : instruments) {
            double df = inst->solveDF(curve, interp);
            curve.addNode(inst->maturity(), df);
        }
        return curve;
    }
};

// ────────────────────────────────────────────────────────────
//  ABSTRACT PRICING INSTRUMENT
//
//  Any instrument that can be priced against a calibrated curve.
//  Separates pricing concerns from calibration concerns
//  (Interface Segregation Principle).
// ────────────────────────────────────────────────────────────

class PricingInstrument {
public:
    virtual ~PricingInstrument() = default;

    // Present value from fixed-rate payer perspective
    virtual double pv(const Curve& curve,
                      const Interpolator& interp) const = 0;

    // Fixed rate that makes PV = 0
    virtual double parRate(const Curve& curve,
                           const Interpolator& interp) const = 0;

    // ∂PV/∂DF(T_j) for every curve node j — building block for risk
    virtual std::vector<double> pvNodeSensitivities(
        const Curve& curve,
        const Interpolator& interp) const = 0;
};

// Concrete: vanilla fixed-for-floating interest rate swap
class VanillaSwap : public PricingInstrument {
    double fixedRate_, maturity_, fixedFreq_, floatFreq_, notional_;

    // Shared helper: annuity sum and DF(T)
    std::pair<double,double> annuityAndDFT(
        const Curve& curve, const Interpolator& interp) const
    {
        std::vector<double> fixedDates;
        for (double t = fixedFreq_; t <= maturity_ + 1e-9; t += fixedFreq_)
            fixedDates.push_back(std::round(t));
        double annuity = 0.0;
        for (int i = 0; i < static_cast<int>(fixedDates.size()); i++) {
            double t_prev = (i == 0) ? 0.0 : fixedDates[i-1];
            annuity += interp.getDF(fixedDates[i], curve)
                       * dcf(fixedDates[i], t_prev);
        }
        return {annuity, interp.getDF(maturity_, curve)};
    }

public:
    VanillaSwap(double fixedRate, double maturity,
                double fixedFreq, double floatFreq,
                double notional = 100.0)
        : fixedRate_(fixedRate), maturity_(maturity),
          fixedFreq_(fixedFreq), floatFreq_(floatFreq),
          notional_(notional) {}

    double pv(const Curve& curve, const Interpolator& interp) const override {
        auto [annuity, dfT] = annuityAndDFT(curve, interp);
        return notional_ * (1.0 - dfT) - notional_ * fixedRate_ * annuity;
    }

    double parRate(const Curve& curve, const Interpolator& interp) const override {
        auto [annuity, dfT] = annuityAndDFT(curve, interp);
        return (1.0 - dfT) / annuity;
    }

    std::vector<double> pvNodeSensitivities(
        const Curve& curve, const Interpolator& interp) const override
    {
        return computePVNodeSensitivities(
            curve, interp.name(), fixedRate_, maturity_, fixedFreq_, notional_);
    }
};

// ════════════════════════════════════════════════════════════
//  CURVE BUILDERS
//  Now implemented using the generic Q3 framework above.
//  The external interface is unchanged — Q1/Q2/Q2.2 code
//  in main() requires zero modification.
// ════════════════════════════════════════════════════════════

static Curve buildCashCurve(const std::vector<Instrument>& insts) {
    std::vector<std::unique_ptr<CalibrationInstrument>> calibInsts;
    for (const auto& inst : insts)
        calibInsts.push_back(
            std::make_unique<CashDeposit>(inst.maturity, inst.cashRate));
    // CashDeposit::solveDF ignores the interpolator; LinearInterpolator used
    // as a placeholder to satisfy the generic interface
    return CurveBootstrapper::calibrate(calibInsts, LinearInterpolator{});
}

static Curve buildSwapCurve(const std::vector<Instrument>& insts,
                             const std::string& method) {
    std::vector<std::unique_ptr<CalibrationInstrument>> calibInsts;
    for (const auto& inst : insts)
        calibInsts.push_back(
            std::make_unique<ParSwap>(inst.maturity, inst.swapRate));
    if (method == "linear")
        return CurveBootstrapper::calibrate(calibInsts, LinearInterpolator{});
    if (method == "aq")
        return CurveBootstrapper::calibrate(calibInsts, AQInterpolator{});
    throw std::invalid_argument("Unknown interpolation method: " + method);
}

// ────────────────────────────────────────────────────────────
//  PRICE SWAP  (Q2.1)
// ────────────────────────────────────────────────────────────

static std::pair<double, double> priceSwap(const Curve& curve,
                                            const std::string& method,
                                            double fixedRate,
                                            double maturityDays,
                                            double fixedFreqDays,
                                            double floatFreqDays,
                                            double notional = 100.0) {
    (void)floatFreqDays; // floating leg uses telescoping identity
    std::vector<double> fixedDates;
    for (double t = fixedFreqDays; t <= maturityDays + 1e-9; t += fixedFreqDays)
        fixedDates.push_back(std::round(t));

    double annuity = 0.0;
    for (int i = 0; i < static_cast<int>(fixedDates.size()); i++) {
        double t_prev = (i == 0) ? 0.0 : fixedDates[i-1];
        annuity += curve.interpolate(fixedDates[i], method) * dcf(fixedDates[i], t_prev);
    }
    double df_T       = curve.interpolate(maturityDays, method);
    double pv         = notional * (1.0 - df_T) - notional * fixedRate * annuity;
    double parSwapRate = (1.0 - df_T) / annuity;
    return {pv, parSwapRate};
}

// ────────────────────────────────────────────────────────────
//  INTERPOLATION SENSITIVITIES  (Q2.2 core building block)
//  ∂DF(t)/∂DF(T_j) for every node j
// ────────────────────────────────────────────────────────────

static std::vector<double> computeInterpSensitivities(
    double t, const Curve& curve, const std::string& method)
{
    int n = curve.size();
    std::vector<double> result(n, 0.0);
    if (n == 0) return result;
    if (t <= curve.times[0])     { result[0]     = 1.0; return result; }
    if (t >= curve.times[n - 1]) { result[n - 1] = 1.0; return result; }

    int i = static_cast<int>(
        std::upper_bound(curve.times.begin(), curve.times.end(), t)
        - curve.times.begin()) - 1;

    if (std::abs(t - curve.times[i])     < 1e-9) { result[i]     = 1.0; return result; }
    if (std::abs(t - curve.times[i + 1]) < 1e-9) { result[i + 1] = 1.0; return result; }

    double df_t = curve.interpolate(t, method);

    if (method == "linear" || n < 3 || i == 0) {
        double w = (t - curve.times[i]) / (curve.times[i+1] - curve.times[i]);
        result[i]     = (1.0 - w) * df_t / curve.dfs[i];
        result[i + 1] = w         * df_t / curve.dfs[i+1];
        return result;
    }

    double wR = (t - curve.times[i])      / (curve.times[i+1] - curve.times[i]);
    double wL = (curve.times[i+1] - t)    / (curve.times[i+1] - curve.times[i]);

    if (i + 2 >= n) {
        auto L = lagrangeBases(t, curve.times[i-1], curve.times[i], curve.times[i+1]);
        result[i-1] = df_t / curve.dfs[i-1] * L[0];
        result[i]   = df_t / curve.dfs[i]   * L[1];
        result[i+1] = df_t / curve.dfs[i+1] * L[2];
        return result;
    }

    auto Li  = lagrangeBases(t, curve.times[i-1], curve.times[i],   curve.times[i+1]);
    auto Li1 = lagrangeBases(t, curve.times[i],   curve.times[i+1], curve.times[i+2]);

    result[i-1] = df_t / curve.dfs[i-1] *  wL * Li[0];
    result[i]   = df_t / curve.dfs[i]   * (wL * Li[1] + wR * Li1[0]);
    result[i+1] = df_t / curve.dfs[i+1] * (wL * Li[2] + wR * Li1[1]);
    result[i+2] = df_t / curve.dfs[i+2] *  wR * Li1[2];
    return result;
}

// ────────────────────────────────────────────────────────────
//  PV NODE SENSITIVITIES  (Q2.2)
//  ∂PV/∂DF(T_j) for every node j
// ────────────────────────────────────────────────────────────

static std::vector<double> computePVNodeSensitivities(
    const Curve& curve, const std::string& method,
    double fixedRate, double maturityDays, double fixedFreqDays,
    double notional)
{
    int n = curve.size();
    std::vector<double> result(n, 0.0);

    std::vector<double> fixedDates;
    for (double t = fixedFreqDays; t <= maturityDays + 1e-9; t += fixedFreqDays)
        fixedDates.push_back(std::round(t));

    for (int k = 0; k < static_cast<int>(fixedDates.size()); k++) {
        double t_prev = (k == 0) ? 0.0 : fixedDates[k-1];
        double weight = notional * fixedRate * dcf(fixedDates[k], t_prev);
        auto   sens   = computeInterpSensitivities(fixedDates[k], curve, method);
        for (int j = 0; j < n; j++) result[j] -= weight * sens[j];
    }

    auto matSens = computeInterpSensitivities(maturityDays, curve, method);
    for (int j = 0; j < n; j++) result[j] -= notional * matSens[j];
    return result;
}

// ────────────────────────────────────────────────────────────
//  CASH RISK VECTOR  (Q2.2)
//  ∂PV/∂c_i — cash nodes are independent, no cascade
// ────────────────────────────────────────────────────────────

static std::vector<double> computeCashRisk(
    const Curve& cashCurve,
    const std::vector<Instrument>& insts,
    const std::vector<double>& pvNodeSens)
{
    int n = static_cast<int>(insts.size());
    std::vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        double dDF_dci = -cashCurve.dfs[i] * cashCurve.dfs[i]
                         * (insts[i].maturity / 360.0);
        result[i] = pvNodeSens[i] * dDF_dci;
    }
    return result;
}

// Helper: semi-annual payment dates for a calibration swap
static std::vector<double> calibSwapDates(double T) {
    std::vector<double> dates;
    if (T <= 180.0) return dates;
    for (double t = 180.0; t <= T + 1e-9; t += 180.0)
        dates.push_back(std::round(t));
    if (!dates.empty()) dates.back() = T;
    return dates;
}

// ────────────────────────────────────────────────────────────
//  SWAP RISK VECTOR  (Q2.2)
//  ∂PV/∂p_k — forward substitution cascade through bootstrap
// ────────────────────────────────────────────────────────────

static std::vector<double> computeSwapRisk(
    const Curve& swapCurve,
    const std::vector<Instrument>& insts,
    const std::string& method,
    const std::vector<double>& pvNodeSens)
{
    int n = static_cast<int>(insts.size());
    std::vector<double> result(n, 0.0);

    for (int k = 0; k < n; k++) {
        double T_k = insts[k].maturity;
        double p_k = insts[k].swapRate;
        std::vector<double> dDF(n, 0.0);

        // Direct term: ∂DF(T_k)/∂p_k
        if (T_k <= 180.0) {
            double denom = 1.0 + p_k * T_k / 360.0;
            dDF[k] = -(T_k / 360.0) / (denom * denom);
        } else {
            auto   dates_k = calibSwapDates(T_k);
            int    nk      = static_cast<int>(dates_k.size());
            double A_k     = 0.0;
            for (int j = 0; j < nk - 1; j++) {
                double t_prev = (j == 0) ? 0.0 : dates_k[j-1];
                A_k += swapCurve.interpolate(dates_k[j], method) * dcf(dates_k[j], t_prev);
            }
            double delta_k = dcf(T_k, dates_k[nk-2]);
            double denom   = 1.0 + p_k * delta_k;
            dDF[k] = -(A_k + delta_k) / (denom * denom);
        }

        // Cascade: ∂DF(T_m)/∂p_k for m > k
        for (int m = k + 1; m < n; m++) {
            double T_m = insts[m].maturity;
            double p_m = insts[m].swapRate;
            if (T_m <= 180.0) { dDF[m] = 0.0; continue; }
            auto   dates_m = calibSwapDates(T_m);
            int    nm      = static_cast<int>(dates_m.size());
            double dA_m    = 0.0;
            for (int j = 0; j < nm - 1; j++) {
                double t_prev = (j == 0) ? 0.0 : dates_m[j-1];
                auto   sens   = computeInterpSensitivities(dates_m[j], swapCurve, method);
                double dot    = 0.0;
                for (int nn = 0; nn < n; nn++) dot += sens[nn] * dDF[nn];
                dA_m += dcf(dates_m[j], t_prev) * dot;
            }
            double delta_m = dcf(T_m, dates_m[nm-2]);
            dDF[m] = (-p_m / (1.0 + p_m * delta_m)) * dA_m;
        }

        for (int m = 0; m < n; m++) result[k] += pvNodeSens[m] * dDF[m];
    }
    return result;
}

// ════════════════════════════════════════════════════════════
//  MAIN
// ════════════════════════════════════════════════════════════

int main() {
    // ── Parse input.csv ───────────────────────────────────
    std::ifstream fin("input.csv");
    if (!fin.is_open()) { std::cerr << "Error: cannot open input.csv\n"; return 1; }

    auto readLine = [&]() -> std::string {
        std::string line;
        while (std::getline(fin, line)) { line = trim(line); if (!line.empty()) return line; }
        return "";
    };

    int numInstruments = std::stoi(splitCSV(readLine())[0]);
    std::vector<Instrument> instruments(numInstruments);
    for (int i = 0; i < numInstruments; i++) {
        auto toks = splitCSV(readLine());
        if (toks.size() < 3)
            throw std::runtime_error("Malformed instrument row " + std::to_string(i));
        instruments[i].maturity = parseMaturity(toks[0]);
        instruments[i].cashRate = std::stod(toks[1]) / 100.0;
        instruments[i].swapRate = std::stod(toks[2]) / 100.0;
    }
    double queryT     = std::stod(splitCSV(readLine())[0]);
    auto   swapToks   = splitCSV(readLine());
    if (swapToks.size() < 4) throw std::runtime_error("Malformed swap row in input");
    double newFixedRate = std::stod(swapToks[0]) / 100.0;
    double newMaturity  = parseMaturity(swapToks[1]);
    double newFixedFreq = parseFrequency(swapToks[2]);
    double newFloatFreq = parseFrequency(swapToks[3]);
    fin.close();

    // ── Build curves (via Q3 generic framework internally) ─
    Curve cashCurve    = buildCashCurve(instruments);
    Curve swapCurveLin = buildSwapCurve(instruments, "linear");
    Curve swapCurveAQ  = buildSwapCurve(instruments, "aq");

    // ── Q1: discount factor at queryT ─────────────────────
    double q1a = cashCurve.interpolate(queryT, "linear");
    double q1b = cashCurve.interpolate(queryT, "aq");
    double q1c = swapCurveLin.interpolate(queryT, "linear");
    double q1d = swapCurveAQ.interpolate(queryT, "aq");

    // ── Q2.1: price new swap (also via Q3 VanillaSwap) ────
    //  Note: priceSwap() and VanillaSwap::pv() produce identical results.
    //  Both paths shown here; output uses the original priceSwap() for
    //  consistency.
    auto [pv_a, par_a] = priceSwap(cashCurve,    "linear", newFixedRate, newMaturity, newFixedFreq, newFloatFreq);
    auto [pv_b, par_b] = priceSwap(cashCurve,    "aq",     newFixedRate, newMaturity, newFixedFreq, newFloatFreq);
    auto [pv_c, par_c] = priceSwap(swapCurveLin, "linear", newFixedRate, newMaturity, newFixedFreq, newFloatFreq);
    auto [pv_d, par_d] = priceSwap(swapCurveAQ,  "aq",     newFixedRate, newMaturity, newFixedFreq, newFloatFreq);

    // ── Q2.2: analytical risk vectors ─────────────────────
    auto pvns_a = computePVNodeSensitivities(cashCurve,    "linear", newFixedRate, newMaturity, newFixedFreq);
    auto pvns_b = computePVNodeSensitivities(cashCurve,    "aq",     newFixedRate, newMaturity, newFixedFreq);
    auto pvns_c = computePVNodeSensitivities(swapCurveLin, "linear", newFixedRate, newMaturity, newFixedFreq);
    auto pvns_d = computePVNodeSensitivities(swapCurveAQ,  "aq",     newFixedRate, newMaturity, newFixedFreq);

    auto risk_a = computeCashRisk(cashCurve,    instruments, pvns_a);
    auto risk_b = computeCashRisk(cashCurve,    instruments, pvns_b);
    auto risk_c = computeSwapRisk(swapCurveLin, instruments, "linear", pvns_c);
    auto risk_d = computeSwapRisk(swapCurveAQ,  instruments, "aq",     pvns_d);

    // ── Write output.csv ───────────────────────────────────
    std::ofstream fout("output.csv");
    if (!fout.is_open()) { std::cerr << "Error: cannot open output.csv\n"; return 1; }
    fout << std::fixed << std::setprecision(10);

    fout << q1a   << "," << q1b   << "," << q1c   << "," << q1d   << "\n";
    fout << pv_a  << "," << pv_b  << "," << pv_c  << "," << pv_d  << "\n";
    fout << par_a << "," << par_b << "," << par_c << "," << par_d << "\n";
    for (int i = 0; i < numInstruments; i++)
        fout << risk_a[i] << "," << risk_b[i] << ","
             << risk_c[i] << "," << risk_d[i] << "\n";

    fout.close();
    std::cout << "Done. Results written to output.csv\n";
    return 0;
}