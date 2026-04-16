import subprocess
from pathlib import Path
import itertools
import pandas as pd
import pyarrow

num_repetitions = 2

root_dir = Path(__file__).resolve(strict=True).parent
ic_dir = root_dir.parent.parent
dataset_dir = ic_dir / "ggbp-instances" / "mtx" / "reduced"
sols_dir = ic_dir / "sols" / "brkga" / "real" / "r3"
sols_dir.mkdir(parents=True, exist_ok=True)

params = {
    "pop": [10, 25, 50, 75, 100, 125, 150, 175, 200],
    "gens": [10, 25, 50, 75, 100, 125, 150, 175, 200],
    "elite": [0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
    "mutant": [0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
    "rhoe": [0.50, 0.60, 0.70, 0.80, 0.90],
    "n_pops": [1, 2, 3],
}

df_dict = {"instance": [], "param": [], "value": [], "result": [], "time": []}

for file in dataset_dir.iterdir():

    for param in params.items():
        print(param[0])
        for value in param[1]:
            curr_sol_dir = sols_dir / file.stem / param[0]
            curr_sol_dir.mkdir(parents=True, exist_ok=True)

            for i in range(num_repetitions):
                result = subprocess.run(
                    [root_dir / "build" / "bin" / "RelWithDebInfo" / "brkga.exe", "mtx",
                    file, "3",
                    f"--{param[0]}", f"{value}",
                    "-o", str(curr_sol_dir) + f"/param_{value}_run{i}.sol"
                    ],
                    capture_output=True,
                    text=True
                )

                out_str = result.stdout.split("\n")
            
                df_dict["instance"].append(file.stem)
                df_dict["param"].append(param[0])
                df_dict["value"].append(value)
                df_dict["result"].append(int(out_str[0]))
                df_dict["time"].append(float(out_str[1]))

df = pd.DataFrame.from_dict(df_dict)
df.to_parquet("r3.parquet.gzip", compression="gzip")
print(pd.read_parquet("r3.parquet.gzip"))