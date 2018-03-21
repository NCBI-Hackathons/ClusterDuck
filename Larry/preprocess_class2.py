"""
Converts Abstract text into a list of terms
"""

import pandas as pd
from nltk.corpus import stopwords
from nltk.tokenize import RegexpTokenizer
from nltk.stem import WordNetLemmatizer


class Preprocessor(object):
    """
    Use nlp techniques to process abstract text.
    The case of the text is lowered and the punctuation is removed
    The words are lemmatize
    Keywords replaced with *keywords*
    """
    def __init__(self):
        pass
    
    def load(self):
        """
        Read in corpus
        """
        abstracts = pd.read_csv(self.corpus, index_col=0)
        abstracts['AB'] = abstracts['AB'].str.lower()
        abstracts['AB'] = abstracts['AB'].str.replace(self.keywords, '')
        return abstracts
    
    def preprocess(self, abstract):
        """
        These are the steps for normalizing the abstract text
        Convert an abstract to word tokens.  This is done by tokenizing the text,
        removing english stopwords and punctuation,and finally lemmatizing the words. 
        
        Args:
            abstract(str): Indivdual abstracts
        
        Return:
            words(list): list of normalized words
        
        """
        STOPWORDS = set(stopwords.words('english'))
        
        # Instantiate Lemmanizer
        WNL = WordNetLemmatizer()

        # tokenize words, remove punctuation
        tokenizer = RegexpTokenizer(r'\w[\w-]+')
        tokens = tokenizer.tokenize(abstract)
        # print(tokens)

        # Remove stopwords and lemmatize tokens
        # Append HPO term to lists
        words = [WNL.lemmatize(word) for word in tokens if word not in STOPWORDS]
        words.append('*{}*'.format(self.keywords))
        return words

    def process(self, keywords, corpus):
        """
        Entry point into Class
        Generates list of terms in each document
        
        Args:
            keywords(str): The keyword term from HPO phenotype
            corpus(str): The path to corpus
        
        Returns:
            documents(Pandas.Series.Series): Series of words in each document
            
        """
        self.corpus = corpus
        self.keywords = keywords
        self.abstracts = self.load()
        documents = self.abstracts.AB.apply(self.preprocess)
        return documents