from pathlib import Path
from urllib.request import urlretrieve

URL = (
    "https://raw.githubusercontent.com/sitonglab/CombinGym/"
    "c5dccd514a9acfe1089a5f1555faea98afd78019/"
    "Data/DMS/Clean/bnAbs_CR9114_H1_clean.xlsx"
)
OUT = Path("bnAbs_CR9114_H1_clean.xlsx")

urlretrieve(URL, OUT)
print(f"downloaded={OUT}")
print(f"bytes={OUT.stat().st_size}")
