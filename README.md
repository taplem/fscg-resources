# Supplementary material for *"A column generation approach for integrated forest management and wood supply chain planning."*

This repository contains:
* Detailed results for all runs (`results.csv`)
* The solver sources (`solver/`)
* A driver script for running the solver (`solve.lua`)
* A script for generating MPS files for LP solvers (`mps.lua`)
* Harvesting constraint definitions (`con-*.txt`)
* Pre-generated small problem instances (`problem/pk100.tar.xz`)
* Data and scripts for generating large problem instances (`problem/*`)
	- Management unit data and initial tree lists (`north-karelia.sqlite3.xz`)
	- A script for simulating the treatment schedules (`sim.lua`)
	- The demand point data (`factories.lua`)

In order to run any of the Lua scripts, you **need** [LuaJIT](https://github.com/luajit/luajit).
The scripts depend on LuaJIT features and do not work on standard Lua.

## Detailed results

The results for each individual run are in `results.csv`.

The result file contains the following columns:
* `problem`: the instance name
* `rmp_solver`: the RMP solver (`highs` or `gurobi`)
* `init_mode`: initialization strategy
	- `g`: greedy initialization
	- `p`: heuristic initialization
	- `P`: sampling initialization (primal)
	- `D`: sampling initialization (dual)
* `run`: the run number (0,1,2,3,4)
* `iter`: number of column generation iterations
* `cols`: number of generated columns
* `init_time`: initialization time (s)
* `pricing_time`: total pricing time (s)
* `rmp_time`: total RMP time (s)

## Problems

### Small problem instances

The small problem instances (`pk100*`) are included pre-generated in this repository.
You can simply decompress them:
```
$ tar xJf problem/pk100.tar.xz
```

### Large problem instances

If you want to solve large problem instances, they must be generated first.
See [problem/README.md](problem/README.md) for instructions.

**NOTE: This requires around ~200GB of disk space.**

## Building

You need a Linux environment, a C compiler, and either [HiGHS](https://github.com/ERGO-Code/HiGHS) or [Gurobi](https://www.gurobi.com/).
Ensure that your C compiler can find the HiGHS or Gurobi header and shared library.
Build `solver/fscg.c` with either `-DFSCG_HIGHS` or `-DFSCG_GUROBI`.

Example build with HiGHS. Replace `{HIGHS_PATH}` with the path containing the HiGHS repository.
```
$ gcc solver/fscg.c -O3 -DNDEBUG -march=native -fPIC -DFSCG_HIGHS -lluajit-5.1 -lhighs -lm -I{HIGHS_PATH}/build -I{HIGHS_PATH}/highs -I{HIGHS_PATH}/highs/interfaces -L{HIGHS_PATH}/build/lib -shared -o fscg.so
```

Example build with Gurobi. Replace `{GUROBI_PATH}` with the path containing the Gurobi install.
```
$ gcc solver/fscg.c -O3 -DNDEBUG -march=native -fPIC -DFSCG_GUROBI -lluajit-5.1 -lgurobi120 -lm -I{GUROBI_PATH}/linux64/include -L{GUROBI_PATH}/linux64/lib -shared -o fscg.so
```

You can check whether the build succeeded with:
```
$ luajit -e 'require "fscg"'
```

To obtain detailed output, compile with `-DFSCG_TRACE`.

## Running

Run:
```
$ luajit solve.lua pk{M}{sy|rc}f{F}{mode}
```
where
* `M` defines the problem size (corresponding `pkM.{xda,cda}` data)
* `sy` or `rc` defines the harvesting constraints (`con-*.txt`)
* `F` defines the demand points (`72` uses all demand points, `34` excludes heating plants and CHPs)
* `mode` is one of:
	- `g`: greedy initialization
	- `p`: heuristic initialization
	- `P`: sampling initialization (primal)
	- `D`: sampling initialization (dual)

For example:
```
$ luajit solve.lua pk1thf72P
```

## Generating MPS files

Run:
```
$ luajit mps.lua pk{M}{sy|rc}f{F}
```
where the problem name syntax is the same as for `solve.lua`.
This will output a file called `pk{M}{sy|rc}f{F}.mps`.
