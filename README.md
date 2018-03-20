# DiseaseClusters

Disease Clustering from Literature Based on Minimal Training Data

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