"""Compiles the generated C code with gcc and compares it with the Python reference.

* float/double: against scipy.signal.sosfilt
* fixed16/fixed32: bit-exact against fdcore.fixed_filter

usage: python test_codegen.py      (needs gcc on PATH and a built ../fdesign)
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import signal

import codegen
import fdcore

SPECS = [
    fdcore.Spec("lowpass", "butterworth", 8000, [400], [800], 3, 40),
    fdcore.Spec("highpass", "chebyshev", 48000, [3000], [2000], 0.5, 60),
    fdcore.Spec("bandpass", "elliptic", 8000, [800, 1200], [400, 1800], 1, 50),
    fdcore.Spec("bandstop", "elliptic", 44100, [900, 1300], [1000, 1150], 0.5, 40),
    fdcore.Spec("lowpass", "elliptic", 16000, [1000], [1300], 0.1, 70),
]
N = 600


def run_c(tmp: Path, hname, header, cname, source, ident, ctype, x) -> np.ndarray:
    (tmp / hname).write_text(header)
    (tmp / cname).write_text(source)
    is_float = ctype in ("float", "double")
    fmt = "%.17g" if is_float else "%lld"
    cast = "(double)" if is_float else "(long long)"
    data = ", ".join(repr(float(v)) if is_float else str(int(v)) for v in x)
    main = (f'#include <stdio.h>\n#include <stdint.h>\n#include "{hname}"\n'
            f"static const {ctype} in[{len(x)}] = {{{data}}};\n"
            f"int main(void) {{ {ident}_state_t st; size_t i; {ident}_init(&st);\n"
            f"  for (i = 0; i < {len(x)}; i++) printf(\"{fmt}\\n\", {cast}{ident}_process(&st, in[i]));\n"
            f"  return 0; }}\n")
    (tmp / "main.c").write_text(main)
    exe = tmp / "t.exe"
    subprocess.run(["gcc", "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror", "-o", str(exe),
                    str(tmp / "main.c"), str(tmp / cname)], check=True)
    out = subprocess.run([str(exe)], capture_output=True, text=True, check=True).stdout.split()
    return np.array([float(v) if is_float else int(v) for v in out])


def main() -> int:
    rng = np.random.default_rng(0)
    failures = 0
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        for spec in SPECS:
            d = fdcore.design(spec)
            label = f"{spec.characteristic} {spec.type} (order {d.digital_order})"
            for arith in codegen.ARITHMETICS:
                for scaled in (False, True):
                    sos = fdcore.scale_sections(d.sos, spec.fs) if scaled else d.sos
                    fixed = None
                    if arith.startswith("fixed"):
                        word = int(arith[5:])
                        fixed = fdcore.quantize(sos, word, fdcore.auto_frac_bits(sos, word))
                    files = codegen.generate(d, "test_filter", arith, sos, fixed)
                    if fixed is None:
                        x = rng.uniform(-1, 1, N)
                        x[0] = 1.0
                        y = run_c(tmp, *files, "test_filter", arith, x)
                        ref = signal.sosfilt(sos, x)
                        err = np.max(np.abs(y - ref)) / max(np.max(np.abs(ref)), 1e-12)
                        ok = err < (1e-3 if arith == "float" else 1e-9)
                        detail = f"rel. error {err:.2e}"
                    else:
                        full = (1 << (fixed.word - 1)) - 1
                        x = np.round(rng.uniform(-0.5, 0.5, N) * full).astype(np.int64)
                        x[0] = full // 2
                        y = run_c(tmp, *files, "test_filter", f"int{fixed.word}_t", x)
                        ref, sat = fdcore.fixed_filter(fixed, x)
                        ok = np.array_equal(y, ref)
                        detail = f"Q{fixed.frac}, {'bit-exact' if ok else 'MISMATCH'}, {sat} saturations"
                    failures += not ok
                    print(f"{'ok  ' if ok else 'FAIL'} {label:40} {arith:8} "
                          f"{'scaled' if scaled else 'plain '}  {detail}")
    print(f"\n{failures} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
