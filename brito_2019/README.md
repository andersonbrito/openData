# Host-virus cophylogenetic analysis

This repository contains data and code used during the analyses resulting in the article entitled *"Time-calibrated tree reconciliations reveal frequent losses, intrahost speciations, and host switches in the evolution of herpesviruses"*.


## Getting Started

The script `beast2jane.py` converts two beast trees (host and virus trees) into a Jane cophylogenetic file. A tab-delimited auxiliary file containing virus-host pairs (one pair pair line) named as in the original trees is required. If any pair is missing from such list, the code will prune such pair from the final converted tree file, as well as output individual versions of the trees without such pairs. The users must define the size of the time zones (in Million years). As in Jane 4 all time zones must be filled with at least one internal node, zones lacking assigned internal nodes will be automatically filled with nodes from an artificial outgroup clade added for this purpose.


### Prerequisites

For running `beast2jane.py`, the following Python 3 package must be installed:


* Phylo (from Biopython)


## Running the script

Basic command line:

```
python beast2Jane.py /Path/to/input/files/ host_beastMCC.tree virus_beastMCC.tree vhPairs.txt 5
```
where the `.tree` files are maximum clade credibility beast trees; `vhPairs.txt` is a tab-delimited auxiliary file containing virus-host pairs (see example directory for more information about file format), and `5` is the size of the time zones (in Million years).


As a result, `beast2jane.py` will output a nexus file in Jane format (`example_GraftTimed.nex`), and pruned versions of the original trees (in case any viral or host taxa are not included in `vhPairs.txt`).

## Authors

* **Anderson Brito** - [GitHub Page](https://github.com/andersonbrito)

## License

This project is licensed under the MIT License.

<!---
--->