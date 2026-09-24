from pathlib import Path
from urllib.request import urlretrieve

URL = (
    "https://raw.githubusercontent.com/sitonglab/CombinGym/"
    "c5dccd514a9acfe1089a5f1555faea98afd78019/"
    "Data/DMS/Clean/eqFP611_clean.xlsx"
)
OUT = Path("eqFP611_clean.xlsx")

urlretrieve(URL, OUT)
print(f"downloaded={OUT}")
print(f"bytes={OUT.stat().st_size}")
