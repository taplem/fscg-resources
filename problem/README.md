## Prerequisities

Around ~200GB of disk space will be needed.

Install the following software:
* Clone [m3](https://github.com/menu-hanke/m3) and follow the build instructions.
* Clone [lim3](https://github.com/menu-hanke/lim3). No build step is needed.
* Install [proj](https://github.com/georust/proj)

You also need [LuaJIT](https://github.com/luajit/luajit), [SQLite3](https://sqlite.org/), and a Python interpreter ([CPython](https://github.com/python/cpython), or preferably, [PyPy](https://github.com/pypy/pypy)).

## Simulating treatment schedules

Start by decompressing the North Karelia forest inventory dataset:
```
$ xz -d north-karelia.sqlite3.xz
```

This dataset already contains the initial tree lists, so you do not need to run MELA2.0.

Pre-create the schedule database and set it to [WAL](https://sqlite.org/wal.html) mode.
This is not required, but makes the next step faster:
```
$ sqlite3 schedules.sqlite3 'PRAGMA journal_mode=WAL'
```

Copy the provided simulator script to the `lim3` directory and run `m3`. **This takes a while and creates a ~100GB SQLite database**:
```
$ cp sim.lua lim3
$ cd lim3
$ m3 sim.lua ../north-karelia.sqlite3 ../schedules.sqlite3
```

Run the `sqlite2xda.py` converter script. **This creates another large file.**
Using `pypy` instead of `python` can speed this up.
```
$ python sqlite2xda.py schedules.sqlite3 pk1.xda
$ python sqlite2xda.py schedules.sqlite3 pk1.cda
```
You may also create samples by setting the `XDA_SCALE` environment variable:
```
$ XDA_SCALE=10 python sqlite2xda.py schedules.sqlite3 pk10.xda
$ XDA_SCALE=10 python sqlite2xda.py schedules.sqlite3 pk10.cda
```
After you have created the `xda` and `cda` files, `schedules.sqlite3` is no longer needed, and can be deleted.

## Generating distance matrices

Run:
```
luajit mkfactories.lua
```

This will generate:
* `pk1_tiem_ka_kp.txt`: Management unit - demand point distance matrix
* `tiem_kp_kp.txt`: Demand point - demand point distance matrix
* `factories_prices.txt`: Demand point data

Note that `tiem_kp_kp.txt` and `factories_prices.txt` are already included in `pk100.tar.xz`, so only `pk1_tiem_ka_kp.txt` is needed.

To run `pk10` instances (or other sampled instances), generate a sampled distance matrix using the `mkfiltertiem.lua` script:
```
$ luajit mkfiltertiem.lua pk10.cda pk10_tiem_ka_kp.txt
```
