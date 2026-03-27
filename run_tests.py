import subprocess
from pathlib import Path

root_dir = Path(__file__).resolve(strict=True).parent
ic_dir = root_dir.parent.parent
dataset_dir = ic_dir / "graphs" / "mtx" / "real"
sols_dir = ic_dir / "sols" / "prym" / "planar"
sols_dir.mkdir(parents=True, exist_ok=True)

for file in dataset_dir.iterdir():
    subprocess.run([root_dir / "build" / "bin" / "prym", "g6", file, "-o", str(sols_dir / file.stem) + ".sol"])