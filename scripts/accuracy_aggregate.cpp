// Aggregate accessibility errors of one method against Raccess Exact over a
// bpRNA sequence list, by Table 1 length bin and window length.
//
//   accuracy_aggregate <seq_list.tsv> <exact_dir> <method_dir> <format> <w,w,...>
//
// format: linracc  (<dir>/<id>.txt, "start\tw,prob;..." lines, '#' metadata)
//         raccess  (<dir>/<id>.txt, "start\tw,energy;...", P = exp(-E/RT))
//         lunp     (<dir>/<id>_lunp, RNAplfold: row = 1-based end i,
//                   column l, start = i - l)
// The exact reference is always Raccess output (<exact_dir>/<id>.txt).
// Prints JSON: for each bin ("100-300", ..., "700-1000", "pooled") and each
// w (and "all"): n_cells, median/p95/max absolute error, and for rare cells
// (exact P in [1e-8, 1e-5]) n and median/p90 relative error; opening-energy
// error |dG_a - dG_e| with dG = -RT ln P (kcal/mol) median/p95 over cells
// where both P > 0, and n_inf = cells with exact P > 0 but approximate P <= 0;
// and n_seq_fail = sequences with some cell of absolute error > 0.5.
// Sequence-list order is the loop order; missing files are counted and
// skipped.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kRT = 0.61633008;
constexpr double kRareLow = 1e-8, kRareHigh = 1e-5;
const int kBinLo[] = {100, 300, 500, 700};
const int kBinHi[] = {300, 500, 700, 1001};
const char* kBinLabel[] = {"100-300", "300-500", "500-700", "700-1000"};

// table[w_index][start] = probability; NaN when absent.
using Table = std::vector<std::vector<double>>;

bool read_start_payload(const std::string& path, bool energy, const std::vector<int>& ws,
                        int n, Table& t) {
    std::ifstream in(path);
    if (!in) return false;
    t.assign(ws.size(), {});
    for (size_t k = 0; k < ws.size(); k++) t[k].assign(n - ws[k] + 1, NAN);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '>' || line[0] == '#') continue;
        size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        int start = std::atoi(line.c_str());
        size_t p = tab + 1;
        while (p < line.size()) {
            size_t semi = line.find(';', p);
            if (semi == std::string::npos) semi = line.size();
            size_t comma = line.find(',', p);
            if (comma != std::string::npos && comma < semi) {
                int w = std::atoi(line.c_str() + p);
                double v = std::strtod(line.c_str() + comma + 1, nullptr);
                for (size_t k = 0; k < ws.size(); k++) {
                    if (ws[k] == w && start >= 0 && start < (int)t[k].size()) {
                        t[k][start] = energy ? (std::isfinite(v) ? std::exp(-v / kRT) : 0.0) : v;
                    }
                }
            }
            p = semi + 1;
        }
    }
    return true;
}

bool read_lunp(const std::string& path, const std::vector<int>& ws, int n, Table& t) {
    std::ifstream in(path);
    if (!in) return false;
    t.assign(ws.size(), {});
    for (size_t k = 0; k < ws.size(); k++) t[k].assign(n - ws[k] + 1, NAN);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == ' ') continue;
        std::vector<std::string> cols;
        std::stringstream ss(line);
        std::string c;
        while (std::getline(ss, c, '\t')) cols.push_back(c);
        int i = std::atoi(cols[0].c_str());
        for (size_t k = 0; k < ws.size(); k++) {
            int l = ws[k];
            if (l >= (int)cols.size() || cols[l] == "NA") continue;
            int start = i - l;
            if (start >= 0 && start < (int)t[k].size()) t[k][start] = std::strtod(cols[l].c_str(), nullptr);
        }
    }
    return true;
}

struct Acc {
    std::vector<double> abs, rel, de;
    size_t n_inf = 0;
};

double quantile(std::vector<double>& v, double q) {
    if (v.empty()) return NAN;
    size_t idx = std::min(v.size() - 1, (size_t)std::ceil(q * v.size()) - (q > 0 ? 1 : 0));
    std::nth_element(v.begin(), v.begin() + idx, v.end());
    return v[idx];
}

double median(std::vector<double>& v) {
    if (v.empty()) return NAN;
    size_t m = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + m, v.end());
    double hi = v[m];
    if (v.size() % 2) return hi;
    double lo = *std::max_element(v.begin(), v.begin() + m);
    return 0.5 * (lo + hi);
}

void num(double x) {
    if (std::isfinite(x)) std::printf("%.6g", x); else std::printf("null");
}

void emit(Acc& a) {
    std::printf("{\"n_cells\":%zu,\"median_abs\":", a.abs.size());
    num(median(a.abs));
    std::printf(",\"p95_abs\":");
    num(quantile(a.abs, 0.95));
    std::printf(",\"max_abs\":");
    num(a.abs.empty() ? NAN : *std::max_element(a.abs.begin(), a.abs.end()));
    std::printf(",\"n_rare\":%zu,\"median_rel_rare\":", a.rel.size());
    num(median(a.rel));
    std::printf(",\"p90_rel_rare\":");
    num(quantile(a.rel, 0.90));
    std::printf(",\"median_dG_kcal\":");
    num(median(a.de));
    std::printf(",\"p95_dG_kcal\":");
    num(quantile(a.de, 0.95));
    std::printf(",\"n_inf\":%zu}", a.n_inf);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::fprintf(stderr, "usage: %s seq_list exact_dir method_dir linracc|raccess|lunp w,w,...\n", argv[0]);
        return 2;
    }
    const std::string list = argv[1], exact_dir = argv[2], method_dir = argv[3], fmt = argv[4];
    std::vector<int> ws;
    {
        std::stringstream ss(argv[5]);
        std::string tok;
        while (std::getline(ss, tok, ',')) ws.push_back(std::atoi(tok.c_str()));
    }
    const int nb = 4, nw = ws.size();
    // acc[bin][w] with bin nb = pooled and w nw = all windows.
    std::vector<std::vector<Acc>> acc(nb + 1, std::vector<Acc>(nw + 1));
    std::vector<int> nseq(nb + 1, 0), nfail(nb + 1, 0);
    int missing = 0, missing_cells = 0;

    std::ifstream lf(list);
    std::string id;
    int n;
    Table ex, ap;
    while (lf >> id >> n) {
        int b = -1;
        for (int k = 0; k < nb; k++) if (n >= kBinLo[k] && n < kBinHi[k]) b = k;
        if (b < 0) continue;
        bool ok = read_start_payload(exact_dir + "/" + id + ".txt", true, ws, n, ex);
        if (fmt == "lunp") ok = ok && read_lunp(method_dir + "/" + id + "_lunp", ws, n, ap);
        else ok = ok && read_start_payload(method_dir + "/" + id + ".txt", fmt == "raccess", ws, n, ap);
        if (!ok) { missing++; continue; }
        nseq[b]++; nseq[nb]++;
        double worst = 0;
        for (int k = 0; k < nw; k++) {
            for (size_t s = 0; s < ex[k].size(); s++) {
                double e = ex[k][s], a = ap[k][s];
                if (!std::isfinite(e) || !std::isfinite(a)) { missing_cells++; continue; }
                double d = std::fabs(a - e);
                worst = std::max(worst, d);
                const bool e_pos = e > 0, a_pos = a > 0;
                const double dg = (e_pos && a_pos) ? kRT * std::fabs(std::log(a) - std::log(e)) : NAN;
                for (int bb : {b, nb}) for (int kk : {k, nw}) {
                    acc[bb][kk].abs.push_back(d);
                    if (e >= kRareLow && e <= kRareHigh) acc[bb][kk].rel.push_back(d / e);
                    if (std::isfinite(dg)) acc[bb][kk].de.push_back(dg);
                    else if (e_pos && !a_pos) acc[bb][kk].n_inf++;
                }
            }
        }
        if (worst > 0.5) { nfail[b]++; nfail[nb]++; }
    }
    std::printf("{\"n_missing_sequences\":%d,\"n_missing_cells\":%d,\"bins\":{", missing, missing_cells);
    for (int bb = 0; bb <= nb; bb++) {
        std::printf("%s\"%s\":{\"n_sequences\":%d,\"n_seq_fail\":%d", bb ? "," : "", bb < nb ? kBinLabel[bb] : "pooled", nseq[bb], nfail[bb]);
        for (int kk = 0; kk <= nw; kk++) {
            std::printf(",\"%s\":", kk < nw ? std::to_string(ws[kk]).c_str() : "all");
            emit(acc[bb][kk]);
        }
        std::printf("}");
    }
    std::printf("}}\n");
    return 0;
}
