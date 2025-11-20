# Filename:     data_get.v2.py

# Description:  Builds SAMOS URL for .csv files necessary for flux variables.
#               Data is already quality controlled, and can be passed an array of ships
#               as well as a year range. Downloaded in

# ----------------------------------------------------------------------
# 1. Necessary python libraries/modules
# ----------------------------------------------------------------------

import urllib.parse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests  # pip install requests
from tqdm import tqdm   # pip install tqdm

# ----------------------------------------------------------------------
# 2. INPUT
# ----------------------------------------------------------------------


ships = ["WTEP"]   # include as many ships as you want in this array
start_year = 2025  # set a start year
end_year = 2025  # and end year
output_dir = 'output/'

pairs = [(s, y) for s in ships for y in range(start_year, end_year + 1)]

# ----------------------------------------------------------------------
# 3. URL builder (this builds the URL for SAMOS with given variables)
# ----------------------------------------------------------------------


def build_working_url(ship: str, year: int) -> str:
    base = "https://erddap-samos.coaps.fsu.edu/erddap/tabledap/SAMOS_Fluxes_B23_v301.csv?"
    variables = "time%2Cplatform_call_sign%2Clatitude%2Clongitude%2Cin_T%2Cin_RH%2Cin_TS%2Cin_SPD%2Cin_P%2Cshf%2Clhf"
    parts = [
        variables,
        f"time%3E={year}-01-01",
        f"time%3C={year}-12-31T23%3A59%3A00Z",
        f"platform_call_sign=%22{urllib.parse.quote(ship)}%22",
        "latitude%3E=-78.64944",
        "latitude%3C=89.99979",
        "longitude%3E=0",
        "longitude%3C=359.9999",
        "in_T_qc=1", "in_RH_qc=1", "in_TS_qc=1", "in_SPD_qc=1", "in_P_qc=1",
        "shf_qc=1", "lhf_qc=1"
    ]
    return base + "&".join(parts)

# ----------------------------------------------------------------------
# 4. Download one file
# ----------------------------------------------------------------------


def download_one(pair: tuple[str, int], out_dir: Path):
    ship, year = pair
    url = build_working_url(ship, year)
    out_file = out_dir / f"{ship}_{year}.csv"

    resp = requests.get(url, timeout=90)
    resp.raise_for_status()
    out_file.write_bytes(resp.content)
    return out_file, resp.content.count(b"\n")

# ----------------------------------------------------------------------
# 5. Parallel driver
# ----------------------------------------------------------------------


def download_all(pairs, out_dir: Path | str = output_dir, workers: int = 12):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    with ThreadPoolExecutor(max_workers=workers) as exe:
        future_to_pair = {exe.submit(
            download_one, p, out_dir): p for p in pairs}
        for fut in tqdm(as_completed(future_to_pair), total=len(pairs), desc="Downloading"):
            path, lines = fut.result()
            results.append((path.name, lines))

    print("\nSummary:")
    for name, lines in sorted(results):
        status = "empty" if lines <= 1 else f"{lines:,} rows"
        print(f"  {name} → {status}")
    print(f"\nAll files saved to: {out_dir.resolve()}")


# ----------------------------------------------------------------------
# 6. Run
# ----------------------------------------------------------------------
if __name__ == "__main__":
    download_all(pairs)
