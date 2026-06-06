# ------------------------------------------------------------
# Get Python package versions for reporting summary
# ------------------------------------------------------------

import sys
import importlib.metadata as metadata
import pandas as pd

packages = [
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "brainspace"
]

versions = []

versions.append({
    "software_or_package": "Python",
    "version": sys.version.split()[0]
})

for pkg in packages:
    try:
        version = metadata.version(pkg)
    except metadata.PackageNotFoundError:
        version = None

    versions.append({
        "software_or_package": pkg,
        "version": version
    })

df_versions = pd.DataFrame(versions)

print(df_versions)

df_versions.to_csv("/Users/zhangyuan/Google Drive/2023_math_reading_neurotransmitter/GitHub/Python_package_versions.csv", index=False)