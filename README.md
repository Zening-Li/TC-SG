# TC-SG
This is the source code of paper "Privacy-Preserving Triangle Counting in Signed Graphs".

## Getting Started
### Required Libraries
* [StatsLib](https://www.kthohr.com/statslib.html)

The directory structure of "./include" is as follows:
* gcem_incl/
* stats_incl/
* gcem.hpp
* stats.hpp
* mt19937ar.h

### Execution Environment
* Ubuntu 22.04.2 LTS
* gcc 11.3.0
* python 3.10.8

### Datasets

Please download the datasets and put them in "./data".

The datasets Wikielections, Epinions, and Wikipolitics can be downloaded from [link](https://github.com/egalimberti/polarized_communities).

The datasets Youtube and Pokec can be downloaded from [link](https://snap.stanford.edu/data/).

### Usage
```shell
# processes data
python data_prep.py
# compile
make
```

#### Centralized Model
Evaluate the privacy and utility of triangle counting algorithms under centralized DP.
```shell
./DP.sh
```

#### Local Model
Evaluate the privacy and utility of triangle counting algorithms under local DP.
```shell
./LDP.sh
```