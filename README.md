# DiseaseClusters

Disease Clustering from Literature Based on Minimal Training Data.

## Presentations

- [Day 1](https://docs.google.com/presentation/d/1OeYWhXnbjgy0pLFU8URxye0xEFQdAGDFmN038SomnbU/edit?usp=sharing)
- [Day 2](https://docs.google.com/presentation/d/1Dgd9E-IKHj1mSZOfUHR4GfukmkeYuXqyeWCGdBRG6Lg/edit?usp=sharing)

## Prerequisite

- Python 3

## Installation

- Install python packages required: `pip3 install -r requirements.txt`

- Download the pubmed database: `python3 setup.py`

- Install data of [nltk](https://www.nltk.org/index.html), enter following codes:
``` python
import nltk
nltk.download('punkt')
nltk.download('stopwords')
```

## Test Suite
`python3 ./dc/test_utils.py`
