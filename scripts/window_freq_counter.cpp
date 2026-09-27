// Fast window-unpaired-frequency counter for LinearSampling's dot-bracket
// output stream. Reads dot-bracket structure lines from stdin (the caller is
// responsible for skipping LinearSampling's leading echoed-sequence line,
// e.g. via `tail -n +2`), accumulates a difference array per requested
// window length w (same idea as this paper's own Algorithm 2: each maximal
// unpaired run of length L contributes to every w<=L via one O(1) update),
// then prefix-sums and prints one line per (w, position): "w pos freq".
//
// Usage:
//   linearsampling -k K -b B <<< "$SEQ" | tail -n +2 | \
//     ./window_freq_counter <length> <w1,w2,...>
//
// This exists because a pure-Python equivalent (even using re.finditer for
// the run-scan, offloading that one part to C) was still ~6x slower than
// LinearSampling's own raw execution time at K=5,000,000 -- Python's
// per-line overhead (file iteration, regex match-object creation, dict
// iteration) dominates at this sample count, not the run-scan itself.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s <length> <w1,w2,...>\n", argv[0]);
        return 1;
    }
    int length = atoi(argv[1]);

    std::vector<int> ws;
    {
        char* tok = strtok(argv[2], ",");
        while (tok) {
            ws.push_back(atoi(tok));
            tok = strtok(nullptr, ",");
        }
    }

    std::vector<std::vector<long long>> diffs(ws.size());
    for (size_t wi = 0; wi < ws.size(); ++wi) {
        int w = ws[wi];
        int n = length - w + 1;
        if (n > 0) diffs[wi].assign(n + 1, 0);
    }

    char* line = nullptr;
    size_t cap = 0;
    ssize_t len;
    long long n_samples = 0;

    while ((len = getline(&line, &cap, stdin)) != -1) {
        // strip trailing newline
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[--len] = '\0';
        }
        if (len == 0) continue;
        n_samples++;

        int i = 0;
        while (i < length) {
            if (line[i] == '.') {
                int j = i;
                while (j < length && line[j] == '.') j++;
                int run_len = j - i;
                for (size_t wi = 0; wi < ws.size(); ++wi) {
                    int w = ws[wi];
                    if (diffs[wi].empty() || run_len < w) continue;
                    int n_positions = run_len - w + 1;
                    diffs[wi][i] += 1;
                    diffs[wi][i + n_positions] -= 1;
                }
                i = j;
            } else {
                i++;
            }
        }
    }
    free(line);

    for (size_t wi = 0; wi < ws.size(); ++wi) {
        int w = ws[wi];
        int n = length - w + 1;
        if (n <= 0) continue;
        long long running = 0;
        for (int pos = 0; pos < n; ++pos) {
            running += diffs[wi][pos];
            printf("%d %d %.10g\n", w, pos, n_samples > 0 ? (double)running / (double)n_samples : 0.0);
        }
    }
    fprintf(stderr, "n_samples=%lld\n", n_samples);
    return 0;
}
