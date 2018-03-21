# ClusterDuck
Disease clustering from phenotypic literature data through Document Understanding, Comprehension and Knowledge

(doi: )


### What's the problem?
Typically, SNPs are studied in terms of a "one disease - one SNP" relationship. This results in researchers and clinicians with deep knowledge of a disease but often incomplete knowledge of all potentially relevant SNPs.


### Why should we solve it?
Knowledge of a larger set of potentially relevant SNPs to a collection of phenotypes would allow finding a novel set of relevant publications. 

### What is ClusterDuck?
ClusterDuck is a tool to automatically identify genetically-relevant publications and returns relevant

## How to use ClusterDuck?

### Prerequisite
- Python 3

### Installation
* Install python packages required: `pip3 install -r requirements.txt`

* Download the pubmed database: `python3 setup.py`

* Install data of [nltk](https://www.nltk.org/index.html), enter following codes:
```
python
import nltk
nltk.download('punkt')
nltk.download('stopwords')
```



## ClusterDuck Workflow
![Workflow Pipeline](https://github.com/neromike/DiseaseClusters/blob/master/pipeline.png "Workflow Pipeline")

## Test Suite
```
python3 ./dc/test_utils.py
```

### Input
Set of phenotypic terms from HPO ontology.

### Workflow
* A 'phenotypic' corpus of literature is extracted from PubMed using the user-input HPO phenotypic terms.
* All SNPs mentioned in the 'phenotypic clusters are idenfified.
* PubMed is queried using the phenotypically-relevant SNPs to extract a second 'phenotypic + genetic' corpus.
* Topic modeling is run on each corpus separately.
* Topic distributions are compared to discover new genetically-inspired and relevant topics.

### Output
A list of novel genetically-related topics to the initial phenotypic input.

### Planned Features
* Synonyms search from user-input
   HPO provides a synonym list for each of their controlled vocabulary terms. This can be incorporated as a preprocessor with the user input to allow
* Make use of hierarchy
   HPO is an ontology of terms and user-input terms are likely to have sub- and super-class terms.
* Filtering different types of research articles
   Optionally add a [PT] query filter to the PubMed query to limit the [types of publications](https://www.ncbi.nlm.nih.gov/books/NBK3827/table/pubmedhelp.T.publication_types/?report=objectonly) returned.
* Use of EMR-type data to build corpus as oppose to PubMed
   An EMR-based corpus is more likely to be associated with diseases (especially to ICD terms) than a PubMed-based corpus.



## People/Team
* Jennifer Dong
* Larry Gray
* Joseph Halstead
* Yi Hsiao
* Wayne Pereanu
* Neelay Trivedi
* Nathan Wan
* Donghui Wu



## Presentations
* (Day 1)[https://docs.google.com/presentation/d/1OeYWhXnbjgy0pLFU8URxye0xEFQdAGDFmN038SomnbU/edit?usp=sharing]
* (Day 2)[https://docs.google.com/presentation/d/1Dgd9E-IKHj1mSZOfUHR4GfukmkeYuXqyeWCGdBRG6Lg/edit?usp=sharing]
* (Day 3)[]
