// Per-loop-type window accessibility from sampled structures (LinearSampling
// dot-bracket lines on stdin), the sampling counterpart of LinRacc -byloop.
//
//   window_freq_counter_byloop <seq_len> <w,w,...>
//
// Each maximal unpaired run [s,e] of a structure gets the loop type of its
// innermost enclosing base pair (i,j), following LinearRaccess's contexts:
//   E  exterior  (no enclosing pair)
//   H  hairpin   (i,j) encloses no pair
//   B  bulge     (i,j) encloses exactly one pair (p,q) and exactly one of
//                p-i-1, j-q-1 is > 0
//   I  internal  (i,j) encloses exactly one pair and both gaps are > 0
//   M  multiloop (i,j) encloses two or more pairs
// and every window of length w inside the run counts once for that type
// (difference arrays, O(n) per structure).  Output, with frequencies over
// all samples: lines "<type> <w> <start> <freq>" for nonzero cells; the
// number of samples goes to stderr as "n_samples=<k>".
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s seq_len w,w,...\n", argv[0]);
        return 2;
    }
    const int n = std::atoi(argv[1]);
    std::vector<int> ws;
    {
        std::stringstream ss(argv[2]);
        std::string t;
        while (std::getline(ss, t, ',')) ws.push_back(std::atoi(t.c_str()));
    }
    const char kTypes[] = {'E', 'H', 'B', 'I', 'M'};
    const int nt = 5, nw = static_cast<int>(ws.size());
    // diff[type][w] over window starts 0..n-w (+1 sentinel)
    std::vector<std::vector<std::vector<long long>>> diff(
        nt, std::vector<std::vector<long long>>(nw));
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < nw; k++) diff[t][k].assign(n - ws[k] + 2, 0);

    std::vector<int> partner(n), stack, children(n), first_child(n), enclosing(n);
    std::string line;
    long long samples = 0;
    while (std::getline(std::cin, line)) {
        if ((int)line.size() < n) continue;
        samples++;
        stack.clear();
        for (int i = 0; i < n; i++) {
            children[i] = 0;
            first_child[i] = -1;
            const char c = line[i];
            enclosing[i] = stack.empty() ? -1 : stack.back();
            if (c == '(') {
                if (!stack.empty()) {
                    if (children[stack.back()]++ == 0) first_child[stack.back()] = i;
                }
                stack.push_back(i);
            } else if (c == ')') {
                const int o = stack.back();
                stack.pop_back();
                partner[o] = i;
                partner[i] = o;
            }
        }
        auto loop_type = [&](int o) {  // o: opening index of the enclosing pair
            if (o < 0) return 0;
            if (children[o] == 0) return 1;
            if (children[o] >= 2) return 4;
            const int p = first_child[o], q = partner[p];
            const bool left = p - o - 1 > 0, right = partner[o] - q - 1 > 0;
            return (left && right) ? 3 : 2;
        };
        for (int s = 0; s < n;) {
            if (line[s] != '.') { s++; continue; }
            int e = s;
            while (e + 1 < n && line[e + 1] == '.') e++;
            const int t = loop_type(enclosing[s]);
            const int len = e - s + 1;
            for (int k = 0; k < nw; k++) {
                if (len < ws[k]) continue;
                diff[t][k][s]++;
                diff[t][k][e - ws[k] + 2]--;
            }
            s = e + 1;
        }
    }
    for (int t = 0; t < nt; t++) {
        for (int k = 0; k < nw; k++) {
            long long run = 0;
            for (int i = 0; i <= n - ws[k]; i++) {
                run += diff[t][k][i];
                if (run) std::printf("%c %d %d %.10g\n", kTypes[t], ws[k], i, (double)run / samples);
            }
        }
    }
    std::fprintf(stderr, "n_samples=%lld\n", samples);
    return 0;
}
