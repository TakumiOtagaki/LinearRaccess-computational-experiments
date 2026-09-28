#!/usr/bin/env python3
"""Run DeepRaccess (hmdlab/DeepRaccess, published FCN_structured.pth, window
length 35) and write Raccess-format output, one file per sequence:
    <outdir>/<id>.txt   ">id" then "start\\t35,energy;" lines
so benchmark/scripts/accuracy_aggregate (format "raccess") can score it.

Uses DeepRaccess's own make_dl / predict / windowed (predict.py), so the
preprocessing (sequence reversal, 440-nt chunks with stride 330 and 55-nt
trimming) is exactly the published one.  Index convention, verified against
the teacher CSVs shipped in sample_data/ (see docs/pi_feedback_20260925.md):
output column j of a length-n sequence is the 35-nt window starting at
forward position n - 35 - j, for j = 0..n-35; later columns are padding.

--mode batch   all sequences, grouped by length so padding stays small;
               records total wall time (incl. model load) and per-group time.
--mode single  one sequence per forward pass (batch size 1), per-sequence
               wall time, for the runtime comparison; --threads 1 for CPU.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEEPRACCESS = REPO_ROOT / "vendor/DeepRaccess"
sys.path.insert(0, str(DEEPRACCESS))

import torch  # noqa: E402
import mymodel  # noqa: E402
import predict as dr  # noqa: E402
import process  # noqa: E402

W = 35


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name is not None:
            seqs[name] += line.strip().upper().replace("T", "U")
    return seqs


def load_model(device, pretrain: Path):
    # Same wrapping as DeepRaccess's process.model_device (DataParallel, so the
    # checkpoint's "module." keys load; deterministic cuDNN), but over the
    # GPUs actually present -- model_device hard-codes device_ids=[0,1,2,3].
    model = getattr(mymodel, "FCN")().to(device)
    ids = list(range(torch.cuda.device_count())) if device.type == "cuda" else None
    model = torch.nn.DataParallel(model, device_ids=ids)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    model.load_state_dict(torch.load(pretrain, map_location=device)["model_state_dict"])
    return model.module.to(device) if hasattr(model, "module") else model.to(device)


def run_group(model, device, names, seqs, batch, tmpdir):
    fa = tmpdir / "group.fa"
    fa.write_text("".join(f">{n}\n{seqs[n]}\n" for n in names))
    dl, flag, division = dr.make_dl(str(fa), batch)
    _, out = dr.predict(device, model, dl)
    return dr.windowed(out, flag, division)


def write_output(outdir: Path, name: str, seq: str, row) -> None:
    n = len(seq)
    lines = [f">{name}"]
    for j in range(n - W + 1):
        lines.append(f"{n - W - j}\t{W},{float(row[j]):.6g};")
    (outdir / f"{name.replace('/', '_')}.txt").write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--ids", type=Path, help="optional list of sequence IDs (first column)")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mode", choices=("batch", "single"), default="batch")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--pretrain", type=Path, default=DEEPRACCESS / "path/FCN_structured.pth")
    args = ap.parse_args()

    if args.threads:
        torch.set_num_threads(args.threads)
    device = torch.device(args.device)
    seqs = read_fasta(args.fasta)
    if args.ids:
        wanted = [l.split()[0] for l in args.ids.read_text().splitlines() if l.strip()]
        seqs = {n: seqs[n] for n in wanted}
    seqs = {n: s for n, s in seqs.items() if len(s) >= W}
    args.outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    model = load_model(device, args.pretrain)
    if device.type == "cuda":
        torch.cuda.synchronize()
    t_load = time.perf_counter() - t0
    record = {"mode": args.mode, "device": str(device), "threads": torch.get_num_threads(),
              "torch": torch.__version__, "n_sequences": len(seqs), "model_load_s": t_load,
              "gpu": torch.cuda.get_device_name(0) if device.type == "cuda" else None}

    with tempfile.TemporaryDirectory(dir=args.outdir) as td:
        tmp = Path(td)
        if args.mode == "batch":
            # Groups share one padded length: <=440 nt, then one group per
            # extra 330-nt chunk, as in DeepRaccess's own windowing.
            groups: dict[int, list[str]] = {}
            for n, s in seqs.items():
                groups.setdefault(max(0, (len(s) - 111) // 330), []).append(n)
            record["groups"] = []
            t1 = time.perf_counter()
            for key in sorted(groups):
                names = groups[key]
                g0 = time.perf_counter()
                out = run_group(model, device, names, seqs, args.batch, tmp)
                if device.type == "cuda":
                    torch.cuda.synchronize()
                record["groups"].append({"key": key, "n": len(names),
                                         "max_len": max(len(seqs[n]) for n in names),
                                         "wall_s": time.perf_counter() - g0})
                for i, n in enumerate(names):
                    write_output(args.outdir, n, seqs[n], out[i])
            record["inference_wall_s"] = time.perf_counter() - t1
        else:
            record["per_sequence"] = []
            # One untimed warm-up pass (CUDA context, kernel selection).
            first = next(iter(seqs))
            run_group(model, device, [first], seqs, 1, tmp)
            for n, s in seqs.items():
                g0 = time.perf_counter()
                out = run_group(model, device, [n], seqs, 1, tmp)
                if device.type == "cuda":
                    torch.cuda.synchronize()
                record["per_sequence"].append({"id": n, "length": len(s),
                                               "wall_s": time.perf_counter() - g0})
                write_output(args.outdir, n, s, out[0])
        if device.type == "cuda":
            record["max_gpu_mem_mb"] = torch.cuda.max_memory_allocated() / 2**20
    record["total_wall_s"] = time.perf_counter() - t0
    (args.outdir / f"run_{args.mode}_{device.type}.json").write_text(json.dumps(record, indent=1))
    print(json.dumps({k: v for k, v in record.items() if k not in ("groups", "per_sequence")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
