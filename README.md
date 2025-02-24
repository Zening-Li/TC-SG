# TC-SG
This is the source code of paper "Triangle Counting over Signed Graphs with Differential Privacy".

## Getting Started
### Required Libraries
* [GCE-Math](https://www.kthohr.com/gcem.html)
* [StatsLib](https://www.kthohr.com/statslib.html)
* [Parallel Hashmap](https://greg7mdp.github.io/parallel-hashmap/)

The directory structure of "./include" is as follows:
* gcem_incl/
* stats_incl/
* parallel_hashmap/
* gcem.hpp
* stats.hpp
* mt19937ar.h

### Execution Environment
* Ubuntu 22.04.2 LTS
* g++ 11.3.0
* python 3.10.8

### Datasets

Please download the datasets and put them in "./data".

The datasets wiki-vote(wikielections), epinions, and wikisigned(wikipolitics) can be downloaded from [link](https://github.com/egalimberti/polarized_communities).

The datasets youtube and pokec can be downloaded from [link](https://snap.stanford.edu/data/).

The dataset dbpedia can be downloaded from [link](http://konect.cc/networks/).

### Usage
```shell
# processes data
python data_prep.py --dataset <dataset_name>
# compile
make
```

#### Centralized Model
Evaluate the privacy and utility of triangle counting algorithms under centralized DP.
```shell
# # methods: CentralGS, CentralSS, CentralSU
# for dbpedia, this method may run out of memory
./DP.sh
```

```shell
# methods: CentralGS, CentralSU^{\star}
./DPSaveMem.sh
```

#### Local Model
Evaluate the privacy and utility of triangle counting algorithms under local DP.
```shell
./LDP.sh
```
