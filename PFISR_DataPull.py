# /usr/bin/env python3
"""
if you haven't used Python Madrigal before, one time do:

python3 -m pip install madrigalWeb
"""
from datetime import datetime
from pathlib import Path

import madrigalWeb.madrigalWeb

url = "http://cedar.openmadrigal.org"

username = "Steven+Cheng"
email = "stcheng@bu.edu"
org = "'Boston+University'"

parms = "ut1_unix,gdalt,elm,f10.7,nel,dnel"
data_dir = Path("~/data_gnss_py").expanduser()
t0 = datetime(2021, 1, 1, 0, 0, 0)
t1 = datetime(2022, 12, 31, 0, 0, 0)
inst_code = 61
format = "Hdf5"
filter = "filter=gdalt,80,500"
timeout = 10.0

data_dir.mkdir(exist_ok=True, parents=True)

madObj = madrigalWeb.madrigalWeb.MadrigalData(url)

print("Getting experiments")
exps = madObj.getExperiments(
    inst_code,
    t0.year,
    t0.month,
    t0.day,
    t0.hour,
    t0.minute,
    t0.second,
    t1.year,
    t1.month,
    t1.day,
    t1.hour,
    t1.minute,
    t1.second,
)
exps.sort()
print(exps[0])

print("Getting experiment files")
expfiles = madObj.getExperimentFiles(exps[0].id)

for exp in exps:
    print(f"Downloading files for {exp.id} {exp.name}")
    for expf in expfiles:
        remote_filename = expf.name
        local_filename = data_dir / expf.name.split("/")[-1]
        print(f"{remote_filename} => {local_filename}")
        madObj.downloadFile(
            remote_filename, local_filename, username, email, org, format="hdf5"
        )
