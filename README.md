# Yield Derivative
Code for Yield Derivative Method, presented:

> Yu Harabuchi, Tomohiko Yokoyama, Wataru Matsuoka, Taihei Oki, Satoru Iwata, and Satoshi Maeda. Differentiating the yield of chemical reactions using parameters in first-order kinetic equations to identify elementary steps that control the reactivity from complicated reaction path networks. *The Journal of Physical Chemistry*, 128(14):2883&mdash;2890, 2024.

## Folder structure
- `cpp`: C++ code for experiments
- `data`: data files
- `result`: experiment results will be saved here

## Running the C++ code
### Requirements
- C++ compiler with C++17 support (tested with Apple Clang 14.0.3)
- [CMake](https://cmake.org/) (version >= 3.30)
- [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/) for Linux and macOS and [pkg-config-lite](https://sourceforge.net/projects/pkgconfiglite/) for Windows
- [Boost](https://www.boost.org/) (version >= 1.7)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) (version >= 3.4)
- [GMP](https://gmplib.org/) (tested with version 6.3.0)

### Compile
```sh
cmake -B build -S cpp -DCMAKE_BUILD_TYPE=Release
make -j -C build
```

### Run

**Population calculation** (`pop`)
```sh
build/bin/pop <DATA> -i <REACTANT_EQ> --time <TIME [sec]> --temperature <TEMPERATURE [K]>
```

**Derivative calculation** (`deriv`)
```sh
build/bin/deriv <DATA> -i <REACTANT_EQ> -p <PRODUCT_EQ> --time <TIME [sec]> --temperature <TEMPERATURE [K]>
```

In both `pop` and `deriv`, the default values are 0 for `-i`, 86400 (sec) for `--time`, and 300 (K) for `--temperature`. The code loads the file `./data/<DATA>` and saves the result to `./result/<DATA>_pop.txt` for `pop` and `./result/<DATA>_deriv.txt"` for `deriv`. Supported input file formats are `*_MinPATH.rrm` and `*_kINP.rrm`.

### Output file format

#### `pop`

```
time,0,1,..,n-1
t_0,x_1(t_0),x_1(t_0),..,x_{n-1}(t_0)
t_1,x_1(t_1),x_1(t_1),..,x_{n-1}(t_1)
...
```

where
- $n$ is the number of EQs
- $t_k$ is the $k$-th reference time $(k=0, 1, \dotsc)$
- $x_i(t_k)$ is the the yield of the $i$-th EQ at time $t_k$ $(0 \le i < n, \, k = 0, 1, \dotsc)$

#### `deriv`

```
n
E(EQ0) D(EQ0)
E(EQ1) D(EQ1)
..
E(EQn-1) D(EQn-1)
m
A0 B0 E(TS0) D(TS0)
A1 B1 E(TS1) D(TS1)
..
Am-1 Bm-1 E(TSm-1) D(TSm-1)
```

where
- $n$ is the number of EQs
- $m$ is the number of TSs
- $E({\mathrm{EQ}i})$ is the energy of the $i$-th EQ $(0 \le i < n)$
- $E({\mathrm{TS}j})$ is the energy of the $j$-th TS $(0 \le j < m)$
- $D({\mathrm{EQ}i})$ is the derivative w.r.t. the $i$-th EQ $(0 \le i < n)$
- $D({\mathrm{TS}j})$ is the derivative w.r.t. the $j$-th TS $(0 \le j < m)$
- $j$-th TS is between the $A_j$-th and $B_j$-th EQs



The default values are 86400 (sec) for `--time` and 300 (K) for `--temperature`. The code loads the file `./data/<DATA>` and saves the result to `./result/data/<DATA>_deriv.txt"`. Supported input file formats are `*_MinPATH.rrm` and `*_kINP.rrm`. The output file is formatted as:

```
n
E(EQ0) D(EQ0)
E(EQ1) D(EQ1)
...
E(EQn-1) D(EQn-1)
m
A0 B0 E(TS0) D(TS0)
A1 B1 E(TS1) D(TS1)
..
Am-1 Bm-1 E(TSm-1) D(TSm-1)
```