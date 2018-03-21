from nltk.corpus import stopwords
from nltk.tokenize import RegexpTokenizer
from nltk.stem import WordNetLemmatizer

STOPWORDS = set(stopwords.words('english'))

# Instantiate Lemmanizer
WNL = WordNetLemmatizer()


def preprocess(abstract, keywords=None):
    """
    Convert an abstract to word tokens.  This is done by lowering the case
    of the text, tokenizing the text, removing english stopwords and
    punctuation,and finally lemmatizing the words. 
    
    Args:
        abstract: (str)
    
    Return:
        str
    """
    # Lowercase all words
    abstract = abstract.lower()
    
    # tokenize words, remove punctuation
    tokenizer = RegexpTokenizer(r'\w[\w-]+')
    tokens = tokenizer.tokenize(abstract)
    
    # Remove stopwords and lemmatize tokens
    words = [WNL.lemmatize(word) for word in tokens if word not in STOPWORDS]
    return words
