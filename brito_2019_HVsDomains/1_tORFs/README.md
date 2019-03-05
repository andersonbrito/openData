# Fetching translated ORFs

`tORFs.py` retrieves peptides encoded by genomes deposited on the NBCI repository. No matter if annotated in the genome or not, all open reading frames (ORFs) longer than 40 amino acids are retrieved, regardless of their initial codon.


## Running `tORFs.py`

### Prerequisites

For running `tORFs.py`, the following Python3 package must be installed:

* Biopython


`tORFs.py` is the first script of the pipeline. It can be executed using a command similar to:

```
python tORFs.py '/path/to/input/files/' 'auxFile'
```

where the `auxFile` is a tab-delimited auxiliary file with lines showing genome accession numbers (e.g. NC_018874) and an species accronyms (e.g. AbHV). Please see the "input" directory for getting an example input file.

Generic file format:
```
NC_018874 \t AbHV
NC_002531 \t AlHV1
NC_024382 \t AlHV2
...
NC_002794 \t TuHV1
```

After running `tORFs.py`, for each genome in `auxFile` a fasta file containing peptide sequences is generated in the working directory. See the `examples` directory for more details about input and output files.

## Author

* **Anderson Brito** - [WebPage](https://andersonbrito.github.io/) - andersonfbrito@gmail.com

## License

This project is licensed under the MIT License.


<!---
### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
--->
