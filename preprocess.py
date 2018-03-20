from nltk.corpus import stopwords
from nltk.tokenize import RegexpTokenizer
from nltk.stem import WordNetLemmatizer
from nltk.tokenize import word_tokenize

STOPWORDS = set(stopwords.words('english'))


def preprocess(abstract):
    # Lowercase all words
    wnl = WordNetLemmatizer()
    abstract = abstract.lower()

    # tokenize words, remove punctuation and stopwords, lemmatize
    tokenizer = RegexpTokenizer(r'[a-zA-Z]\w+\'?\w*|\-\w+')
    # tokenizer = RegexpTokenizer(r'[\w-]')
    tokens = tokenizer.tokenize(abstract)
    # tokens = word_tokenize(abstract)
    words = [wnl.lemmatize(word) for word in tokens if word not in STOPWORDS]
    return words
