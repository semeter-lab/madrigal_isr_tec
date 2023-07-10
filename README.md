# Madrigal TEC analysis

The Matlab script TECvsTime_Uni.m is the top-level program.
It focuses on indexing the tabular data using vectors--database-like operations--which can be 100,000x faster than iterating over the rows.

## Madrigal Data Downloads

I tried fixing the seemingly unmaintained Matlab Madrigal driver, but found the maintained Python Madrigal driver was about 100x faster and much more reliable.

```sh
python -m pip install madrigalWeb
# one time setup

python PFISR_DataPull.py
```
