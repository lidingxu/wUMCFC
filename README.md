

# wUMCFC: a solver for wireless unsplittable multi-commodity flow routing with network coding

`wUMCFC` is a specialized **branch-and-price** (BP) based solver for **unsplittable multi-commodity flow routing with network coding**.  Unsplittable multi-commodity flow is a fundamentel problem in network optimization. `wUMCFC` uses network coding to decrease the *routing cost* and increase the *capacity usage* for wireless networks. `wUMCFC` is written in C++ based on [SCIP](https://www.scipopt.org/).



## Instances
We create a data format `ntrt` to describe the network and interference. `ntrt` format has 5 groups of descriptors, and these 5 groups are listed one by one in the data file.

The first line (group) lists the numbers of vertices, edges, demands and interference cliques.
```
#nv #ne #nf #nc
```
In the second group, each line describes the cost of the  i-th vertex, from the 1st to  the nv-th vertex.  Currently, the vertex cost is not used.

In the third group,  each line describes the tail vertex, head vertex, cost, capacity of the i-th edge, from the 1st to the ne-th vertex. 

In the fourth group, every two lines desribe the size of interference clique and the vertices in the i-th clique, from the 1st to the nc-th clique.

In the fifth group, each line describes the source vertex, the target vertex and the flow of the i-th demand, from the 1st to the nf-th demand. 

Instances are located in `data/benchmark1`,  `data/benchmark2` and  `data/benchmark1_rlp`.
Files in `data/benchmark1` and `data/benchmark1`are in `ntrt` format. `data/benchmark1_rlp` consists of `rlp` format files converted from `data/benchmark1`, `data/benchmark1` and `data/benchmark1_rlp` are used to test `wUMCFC` with other solvers, e.g., SCIP and CPLEX; `data/benchmark2` is used to test the effect of network coding.


## Installation
`wUMCFC` has been developped and tested in *Linux* System. 

*Building* from sources requires:
- a SCIP installation.
- CMake.
- C++ compiler.

### Build
You may use the command line to configure and build `wUMCFC` by creating a `build` directory at the project directory and then building the configuration:
```
mkdir build
```
```
cd  build
```
If your SCIP 's build type is *Release*,  the following command builds and sets the compiler's optimization level to **O3** : 
```
camke .. -DCMAKE_BUILD_TYPE=Release
```
Otherwise your SCIP's build type is *Debug*,  the following command builds but does not enable compilor's optimization : 
```
camke .. -DCMAKE_BUILD_TYPE=Debug
```
The complier generates a binary exectuable `wUMCFC`.

## Usage

You can run `wUMCFC`  in the `build` directory via the following command:
```
./wUMCFC
```
Set the #timelimit:
```
set  limits/time #timelimit
```
The network coding functionality is enabled by default, and you can disable network coding:
```
set  ntrt/netcoding FALSE
```
Then, `wUMCFC` solves a conventional wireless UMCF problem without network coding.

Read a test #instance in `ntrt` format:
```
read #instance.ntrt
```
Optimize:
```
optimize
```

## Example
There is an instance in `build` directory, you can read this example
```
read random.30-1.6-24-2019-09-14-11_38_20.ntrt
```
and 
```
optimize
```
Energy cost of the optimal routing is 3.115.

Now, disable network routing
```
set  ntrt/netcoding FALSE
```
Read and run the solver again
```
read random.30-1.6-24-2019-09-14-11_38_20.ntrt
```
and 
```
optimize
```
The problem is infeasible without network routing!

## Format conversion
Directory `ntrt_bb` contains source code for a standalone branch-and-bound solver `ntrt_bb` based on SCIP.

Installation of `ntrt_bb` is the same as `wUMCFC`.

`ntrt_bb` can convert models (with network coding) from `ntrt` format instances into files of other formats.

Read an #instance in `ntrt` format:
```
read #instance.ntrt
```
Write the problem
```
write problem #instance.ext
```
where the output problem's format is given by file extension, e.g., ext in {lp,rlp,cip,mps}.

## References

If you find `wUMCFC` useful in your work, we kindly request that you cite the following paper draft, which is recommended reading for advanced users:


        @article{xubranch,
                author = {Xu, Liding and Haddad Vanier, Sonia},
                title = {Branch-and-price for energy optimization in multi-hop wireless sensor networks},
                journal = {Networks},
                volume = {80},
                number = {1},
                pages = {123-148},
                doi = {https://doi.org/10.1002/net.22083},
                year = {2022}
        }

