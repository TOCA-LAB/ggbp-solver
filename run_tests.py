import subprocess
from pathlib import Path

num_repetitions = 30

root_dir = Path(__file__).resolve(strict=True).parent
ic_dir = root_dir.parent.parent
dataset_dir = ic_dir / "instances" / "mtx" / "real"
sols_dir = ic_dir / "sols" / "prym" / "planar"
sols_dir.mkdir(parents=True, exist_ok=True)

for file in dataset_dir.iterdir():
    for i in range(num_repetitions):
        print(str(sols_dir / file.stem) + (f"_{i}" if num_repetitions > 1 else "") + ".sol")
        # subprocess.run([root_dir / "build" / "bin" / "prym", "g6", file, "-o", str(sols_dir / file.stem) + f"{i if i > 1 else ""}.sol"])