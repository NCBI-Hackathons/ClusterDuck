from nltk.corpus import stopwords
from nltk.stem import WordNetLemmatizer
from nltk.tokenize import RegexpTokenizer


STOPWORDS = stopwords.words('en')
wnl = WordNetLemmatizer()

abstract = "We found gene variants which may become early predictors of the \
        therapy outcome and allow development of new early prognostic tests \
        for estimation of therapy efficacy in CML patients. Normal genetic \
        variation may influence therapy efficacy during targeted treatment \
        of cancers."


def preprocess(abstract):
    # Lowercase all words
    abstract = abstract.lower()

    # tokenize words, remove punctuation and stopwords, lemmatize
    tokenizer = RegexpTokenizer(r'\w+')
    tokens = tokenizer.tokenize(abstract)
    words = [wnl(word) for word in tokens if word in set(STOPWORDS)]
    return words
