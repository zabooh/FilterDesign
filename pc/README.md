# FilterDesign PC tools

This folder contains PC versions of the FilterDesign library:

- **`fdesign/`** is a command line tool written in C. It wraps the design library and builds
  with a Makefile.
- **`gui/`** is a browser GUI written in Python with NiceGUI. It uses `fdesign` for the design,
  plots the usual filter characteristics and generates C code for the filter.

```
pc/
├── fdesign/
│   ├── Makefile
│   ├── lib/filter_design.c/.h   design library (copy of ../FilterDesign, see "Changes")
│   ├── src/fdesign.c/.h         specification front end (band edges -> F1..F4, SOS output)
│   ├── src/main.c               command line interface
│   └── test/verify.py           verification against SciPy
└── gui/
    ├── app.py                   NiceGUI application
    ├── fdcore.py                runs fdesign, analysis, fixed-point model
    ├── codegen.py               C code generator
    ├── test_codegen.py          compiles the generated code and compares it with Python
    └── requirements.txt
```

## Command line tool `fdesign`

### Build

You need a C99 compiler and GNU make. On Windows this means MinGW-w64 (for example WinLibs) or
MSYS2.

```sh
cd pc/fdesign
make                 # Windows/MinGW: mingw32-make
make REAL=float      # optional: float arithmetic like the PIC32 firmware
make test            # optional: verify ~370 designs against SciPy (needs Python + scipy, ~2 min)
```

The PC build computes in `double` by default. With `float`, narrow filters can miss the
specification. For example, a Chebyshev lowpass at 100 Hz with fs = 48 kHz reaches 0.51 dB
ripple instead of 0.50 dB.

### Usage

```
fdesign -t TYPE -c CHAR -r FS -p FPASS -s FSTOP [--ap DB] [--as DB] [-f text|json|sos] [--no-normalize]
```

| Option | Meaning |
|---|---|
| `-t` | `lowpass`, `highpass`, `bandpass`, `bandstop` (or `lp`, `hp`, `bp`, `bs`) |
| `-c` | `butterworth`, `chebyshev`, `elliptic` (or `butter`, `cheby`, `ellip`, `cauer`) |
| `-r` | Sampling rate in Hz |
| `-p` | Passband edge(s) in Hz. Use two values `f1,f2` for bandpass and bandstop. |
| `-s` | Stopband edge(s) in Hz. Use two values `f1,f2` for bandpass and bandstop. |
| `--ap` | Maximum passband attenuation in dB (default 1) |
| `--as` | Minimum stopband attenuation in dB (default 40) |
| `-f` | `text` (default), `json` (all results), `sos` (one line per section: `b0 b1 b2 a0 a1 a2`) |
| `--no-normalize` | Keep the gain from the library. By default the peak gain is set to 0 dB. |

Exit codes: 0 means OK, 1 a usage error, 2 an invalid specification and 3 a failed design (for
example when the order exceeds the limit).

```
$ fdesign -t bp -c elliptic -r 8000 -p 800,1200 -s 400,1800 --ap 3 --as 40
Filter        : elliptic bandpass
Sampling rate : 8000 Hz
Passband edge : 800 1200 Hz (max. 3 dB)
Stopband edge : 400 1800 Hz (min. 40 dB)
Order         : 6 (3 sections)
Achieved      : pass 1.47 1.47 dB, stop 57.14 52.81 dB -> specification met

  #  ord          b0               b1               b2               a1               a2
  1  2   1.287123799e-01  0.000000000e+00 -1.287123799e-01 -1.339694022e+00  8.712876201e-01
  ...
```

The result is a cascade of sections
`H(z) = Π (b0 + b1 z⁻¹ + b2 z⁻²) / (1 + a1 z⁻¹ + a2 z⁻²)`. The `sos` output uses the same row
layout as `scipy.signal.sosfilt`.

Limits: up to 16 sections. That is order 28 for lowpass and highpass, and order 2×16 for
bandpass and bandstop.

## GUI

```sh
cd pc/gui
pip install -r requirements.txt
python app.py                 # opens http://127.0.0.1:8090 in the browser
python app.py --port 9000 --no-browser
python app.py --host 0.0.0.0  # allow access from other machines
python app.py --light         # start in light mode (default: dark)
```

The sun/moon button in the header switches between dark and light mode. Plots and colors
follow the chosen mode.

The GUI looks for the `fdesign` executable in `../fdesign`. You can set the environment variable
`FDESIGN` to use another path.

You enter the specification (type, characteristic, sampling rate, band edges, ripple, stopband
attenuation) and the implementation (arithmetic, section scaling, Q format, C name). The design
is recomputed on every change.

| Tab | Content |
|---|---|
| Magnitude | Magnitude response in dB or linear, with linear or log frequency axis. The tolerance scheme is shaded red, and the quantized implementation is drawn as a dashed line. |
| Phase | Phase response, unwrapped or wrapped |
| Group delay | In samples or in ms |
| Pole / zero | z-plane with unit circle, multiplicities, and a list of radius and frequency |
| Impulse / step | Time responses of the design and of the implementation (fixed point is simulated bit-exactly) |
| Sections | Coefficients (fixed point as integers), pole radius, frequency and Q, peak gain and L1 norm at each section output |
| C code | Generated `.h` and `.c` for download |

The summary above the tabs shows the order, the achieved attenuation at the band edges, the
maximum pole radius, the settling time, the group delay and the attenuation reached by the
quantized implementation. It also shows the number of saturations in a fixed-point impulse test.

Zoom: use the mouse wheel on the frequency axis, box zoom from the toolbox, or the slider below
the plot. `?tab=phase` (or `group`, `pole`, `impulse`, `sections`, `c`) in the URL opens a tab
directly.

### Generated C code

| Arithmetic | Structure | Types |
|---|---|---|
| `float`, `double` | Transposed direct form II | Coefficients and state in the chosen type |
| `fixed32`, `fixed16` | Direct form I with rounding and saturation | `int32_t` / `int16_t` for samples and coefficients, `int64_t` for the accumulator |

Each filter provides the API `<name>_init()`, `<name>_process()` (one sample) and
`<name>_process_block()`. For fixed point, the code assumes an arithmetic right shift of negative
`int64_t` values. GCC, Clang, XC16, XC32, IAR and Keil all do this.

**Section scaling** (on by default) distributes the gain over the cascade so that the peak gain
at every section output is 0 dB. Without it, internal nodes of high-order elliptic filters
overflow in fixed point. The test showed several hundred saturations, and none with scaling.
**Auto Q** chooses the largest number of fractional bits for which all coefficients fit.

## Verification

- `fdesign/test/verify.py` designs 4 fixed and 120 random specifications with all three
  characteristics. Each design is checked independently with `scipy.signal.sosfreqz` for
  passband ripple, stopband attenuation, stability and 0 dB peak gain, and the order is
  compared with `buttord`, `cheb1ord` and `ellipord`. Result: 358 passed, 0 failed. 14 designs
  were rejected because their order exceeds the limit.
- `gui/test_codegen.py` generates code for 5 designs × 4 arithmetics × with and without scaling,
  compiles it with `gcc -Wall -Wextra -Werror` and compares the output. float and double are
  compared with `scipy.signal.sosfilt` (relative error below 1e-5 and exactly 0). Fixed point is
  compared bit-exactly with the Python model the GUI uses.

## Changes compared with the firmware library

`fdesign/lib/filter_design.c` is a copy of `../FilterDesign/filter_design.c` with these fixes:

1. Butterworth and Chebyshev bandstop computed the stopband edge with the bandpass formula, so
   the order was too low and the stopband attenuation was missed (36–38 dB instead of 40 dB).
2. The Cauer order was rounded up one step too far (`ceil(m + .5)` changed to `ceil(m)`).
3. The Cauer zero reordering read outside the array. The loop condition `i <= Sub_Filter || max >= 1`
   is now `max >= 1`.
4. Array overflows were reported as success (`return FALSE` = 0). There are now error codes, an
   order limit (`FLD_MAX_ORDER`), and errors from the BP/BS transformation are passed on.
5. All arithmetic uses the type `fld_real` (`double` on the PC), and duplicate prototypes were
   removed.

The front end (`src/fdesign.c`) adds:
- the mapping from passband and stopband edges to F1–F4
- for bandpass and bandstop, the stricter of the two stopband edges (the library itself only
  uses one edge and assumes geometric symmetry)
- removal of the pole/zero pair at z = −1 that the library uses to express first-order
  sections as biquads
- normalization of the peak gain to 0 dB

Known limitation: for bandstop filters, SciPy sometimes finds a lower order because it shifts
the passband edges by optimization. The designs here meet the specification, but are not always
minimal.
