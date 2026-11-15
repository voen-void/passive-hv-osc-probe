"""Regenerate the schematic and PCB images used in the README with kicad-cli.

    python render_kicad.py              # all outputs into ../images/kicad/
    python render_kicad.py --fast       # basic-quality 3D renders (seconds instead of minutes)
    python render_kicad.py --only sch   # just the schematic (sch, layout, assembly, 3d)

Outputs:
    schematic.pdf / schematic.svg       full schematic
    layout.svg                          copper, board outline, shield (User.1) and resin (User.2)
    assembly-top.pdf / assembly-bottom.pdf   fab + silkscreen layers (bottom mirrored)
    pcb-top.png / pcb-bottom.png / pcb-iso.png   3D renders

kicad-cli is looked up on PATH, in the usual install folders, or via
[kicad] cli = "..." in probe.toml.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROJECT = os.path.join(ROOT, "hv-probe-passive", "hv-probe-passive")


def find_cli():
    try:
        with open(os.path.join(HERE, "probe.toml"), "rb") as f:
            exe = tomllib.load(f).get("kicad", {}).get("cli", "")
        if exe and os.path.exists(exe):
            return exe
    except OSError:
        pass
    exe = shutil.which("kicad-cli")
    if exe:
        return exe
    def version(p):   # newest install first: "10.0" > "9.0" numerically, not alphabetically
        v = os.path.basename(os.path.dirname(os.path.dirname(p)))
        return tuple(int(x) if x.isdigit() else 0 for x in v.split("."))
    cands = sorted(glob.glob(r"C:\Program Files\KiCad\*\bin\kicad-cli.exe"), key=version, reverse=True)
    cands += ["/Applications/KiCad/KiCad.app/Contents/MacOS/kicad-cli", "/usr/bin/kicad-cli"]
    return next((c for c in cands if os.path.exists(c)), None)


def run(cli, args, what):
    print(f"  {what} ...", end=" ", flush=True)
    r = subprocess.run([cli, *args], capture_output=True, text=True)
    ok = r.returncode == 0
    print("ok" if ok else "FAILED")
    if not ok:
        print("   ", (r.stderr or r.stdout).strip().replace("\n", "\n    "))
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(ROOT, "images", "kicad"))
    ap.add_argument("--fast", action="store_true", help="basic-quality 3D renders")
    ap.add_argument("--only", choices=["sch", "layout", "assembly", "3d"], nargs="+")
    a = ap.parse_args()

    cli = find_cli()
    if not cli:
        sys.exit("kicad-cli not found; set [kicad] cli = '...' in probe.toml")
    os.makedirs(a.out, exist_ok=True)
    sch, pcb = PROJECT + ".kicad_sch", PROJECT + ".kicad_pcb"
    out = lambda name: os.path.join(a.out, name)
    want = lambda k: not a.only or k in a.only
    print(f"kicad-cli: {cli}\noutput:    {a.out}")

    if want("sch"):
        run(cli, ["sch", "export", "pdf", "-o", out("schematic.pdf"), sch], "schematic PDF")
        with tempfile.TemporaryDirectory() as tmp:
            if run(cli, ["sch", "export", "svg", "-e", "-o", tmp, sch], "schematic SVG"):
                svgs = glob.glob(os.path.join(tmp, "*.svg"))
                if svgs:
                    shutil.move(svgs[0], out("schematic.svg"))

    if want("layout"):
        run(cli, ["pcb", "export", "svg", "--mode-single", "--page-size-mode", "2", "--exclude-drawing-sheet",
                  "-l", "Edge.Cuts,F.Cu,B.Cu,User.1,User.2", "-o", out("layout.svg"), pcb],
            "layout SVG (copper, outline, shield, resin)")

    if want("assembly"):
        run(cli, ["pcb", "export", "pdf", "--mode-single", "--ibt", "--sp",
                  "-l", "F.Fab,F.SilkS,Edge.Cuts", "-o", out("assembly-top.pdf"), pcb], "assembly drawing, top")
        run(cli, ["pcb", "export", "pdf", "--mode-single", "--ibt", "--sp", "-m",
                  "-l", "B.Fab,B.SilkS,Edge.Cuts", "-o", out("assembly-bottom.pdf"), pcb], "assembly drawing, bottom")

    if want("3d"):
        q = ["--quality", "basic" if a.fast else "high"]
        # the board is long and thin: zoom in so it fills a wide frame
        common = ["--width", "3200", "--height", "560", "--zoom", "3.4", "--background", "opaque", *q]
        run(cli, ["pcb", "render", "--side", "top", *common, "-o", out("pcb-top.png"), pcb], "3D render, top")
        run(cli, ["pcb", "render", "--side", "bottom", *common, "-o", out("pcb-bottom.png"), pcb], "3D render, bottom")
        # looking at the component (bottom) side; angles must be positive for the CLI
        run(cli, ["pcb", "render", "--perspective", "--rotate", "125,0,25", "--zoom", "1.8",
                  "--width", "2400", "--height", "1000", "--background", "opaque", *q,
                  "-o", out("pcb-iso.png"), pcb], "3D render, isometric")


if __name__ == "__main__":
    main()
