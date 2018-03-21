
# coding: utf-8

# In[1]:


import pandas as pd
from spacy.lang.en.stop_words import STOP_WORDS
# Read and dedupe corpora
initial =  pd.read_csv('original.csv')
initial['tmp'] = initial['AB']
initial['AB'] = initial['PMID']
initial['PMID'] = initial['tmp']
initial = initial[['PMID', 'AB']].groupby('PMID').first()['AB']


# In[2]:



combined =  pd.read_csv('expand.csv')
combined['tmp'] = combined['AB']
combined['AB'] = combined['PMID']
combined['PMID'] = combined['tmp']
combined = combined[['PMID', 'AB']].groupby('PMID').first()['AB']


# In[3]:


# Preprocess - https://github.com/TropComplique/lda2vec-pytorch/blob/master/utils/preprocess.py

from collections import Counter
from tqdm import tqdm


def preprocess(docs, nlp, min_length, min_counts, max_counts):
    """Tokenize, clean, and encode documents.
    Arguments:
        docs: A list of tuples (index, string), each string is a document.
        nlp: A spaCy object, like nlp = spacy.load('en').
        min_length: An integer, minimum document length.
        min_counts: An integer, minimum count of a word.
        max_counts: An integer, maximum count of a word.
    Returns:
        encoded_docs: A list of tuples (index, list), each list is a document
            with words encoded by integer values.
        decoder: A dict, integer -> word.
        word_counts: A list of integers, counts of words that are in decoder.
            word_counts[i] is the number of occurrences of word decoder[i]
            in all documents in docs.
    """

    def clean_and_tokenize(doc):
        text = ' '.join(doc.split())  # remove excessive spaces
#         text = nlp(text, tag=True, parse=False, entity=False)
        text = nlp(text)
        return [t.lemma_ for t in text
                if t.is_alpha and len(t) > 2 and not t.is_stop and not t in STOP_WORDS]

    tokenized_docs = [(i, clean_and_tokenize(doc)) for i, doc in tqdm(docs)]

    # remove short documents
    n_short_docs = sum(1 for i, doc in tokenized_docs if len(doc) < min_length)
    tokenized_docs = [(i, doc) for i, doc in tokenized_docs if len(doc) >= min_length]
    print('number of removed short documents:', n_short_docs)

    # remove some tokens
    counts = _count_unique_tokens(tokenized_docs)
    tokenized_docs = _remove_tokens(tokenized_docs, counts, min_counts, max_counts)
    n_short_docs = sum(1 for i, doc in tokenized_docs if len(doc) < min_length)
    tokenized_docs = [(i, doc) for i, doc in tokenized_docs if len(doc) >= min_length]
    print('number of additionally removed short documents:', n_short_docs)

    counts = _count_unique_tokens(tokenized_docs)
    encoder, decoder, word_counts = _create_token_encoder(counts)

    print('\nminimum word count number:', word_counts[-1])
    print('this number can be less than MIN_COUNTS because of document removal')

    encoded_docs = _encode(tokenized_docs, encoder)
    return encoded_docs, decoder, word_counts


def _count_unique_tokens(tokenized_docs):
    tokens = []
    for i, doc in tokenized_docs:
        tokens += doc
    return Counter(tokens)


def _encode(tokenized_docs, encoder):
    return [(i, [encoder[t] for t in doc]) for i, doc in tokenized_docs]


def _remove_tokens(tokenized_docs, counts, min_counts, max_counts):
    """
    Words with count < min_counts or count > max_counts
    will be removed.
    """
    total_tokens_count = sum(
        count for token, count in counts.most_common()
    )
    print('total number of tokens:', total_tokens_count)

    unknown_tokens_count = sum(
        count for token, count in counts.most_common()
        if count < min_counts or count > max_counts
    )
    print('number of tokens to be removed:', unknown_tokens_count)

    keep = {}
    for token, count in counts.most_common():
        keep[token] = count >= min_counts and count <= max_counts

    return [(i, [t for t in doc if keep[t]]) for i, doc in tokenized_docs]


def _create_token_encoder(counts):

    total_tokens_count = sum(
        count for token, count in counts.most_common()
    )
    print('total number of tokens:', total_tokens_count)

    encoder = {}
    decoder = {}
    word_counts = []
    i = 0

    for token, count in counts.most_common():
        # counts.most_common() is in decreasing count order
        encoder[token] = i
        decoder[i] = token
        word_counts.append(count)
        i += 1

    return encoder, decoder, word_counts


# In[4]:



MIN_COUNTS = 20
MAX_COUNTS = 1800
# words with count < MIN_COUNTS
# and count > MAX_COUNTS
# will be removed

MIN_LENGTH = 15
# minimum document length 
# (number of words)
# after preprocessing

import spacy

nlp = spacy.load('en')





# In[5]:


# encode initial corpus
initial_abstracts = zip(initial.index, initial)

# (encoded_docs, decoder, word_counts)
initial_preprocess = preprocess(
    initial_abstracts, nlp, MIN_LENGTH, MIN_COUNTS, MAX_COUNTS
)


# In[6]:


# encode combined corpus
combined_abstracts = zip(combined.index, combined)

# (encoded_docs, decoder, word_counts)
combined_preprocess = preprocess(
    combined_abstracts, nlp, MIN_LENGTH, MIN_COUNTS, MAX_COUNTS
)


# In[7]:


from gensim import corpora, models

n_topics = 20


# In[8]:


# make single dictionary
(encoded_docs, decoder, word_counts) = initial_preprocess
initial_texts = [[decoder[j] for j in doc] for i, doc in encoded_docs]
(encoded_docs, decoder, word_counts) = combined_preprocess
combined_texts = [[decoder[j] for j in doc] for i, doc in encoded_docs]

dictionary = corpora.Dictionary(initial_texts + combined_texts)

alpha = 0.5
eta = 0.2


# In[9]:


get_ipython().run_cell_magic('time', '', '\n# train initial lda\n(encoded_docs, decoder, word_counts) = initial_preprocess\ntexts = [[decoder[j] for j in doc] for i, doc in encoded_docs]\ncorpus = [dictionary.doc2bow(text) for text in texts]\n\ninitial_lda = models.LdaModel(corpus, alpha=alpha, eta=eta, id2word=dictionary, num_topics=n_topics)')


# In[10]:


get_ipython().run_cell_magic('time', '', '\n# train combined lda\n(encoded_docs, decoder, word_counts) = combined_preprocess\ntexts = [[decoder[j] for j in doc] for i, doc in encoded_docs]\ncorpus = [dictionary.doc2bow(text) for text in texts]\n\ncombined_lda = models.LdaModel(corpus, alpha=alpha, eta=eta, id2word=dictionary, num_topics=n_topics)')


# In[11]:



for i, topics in initial_lda.show_topics(n_topics, formatted=False):
    print('topic', i, ':', ' '.join([t for t, _ in topics]))


# In[12]:



for i, topics in combined_lda.show_topics(n_topics, formatted=False):
    print('topic', i, ':', ' '.join([t for t, _ in topics]))


# In[13]:


# accumulate weights of top 5 topic words
from collections import defaultdict
cum1 = defaultdict(float)
for i in range(5):
    for w, s in initial_lda.show_topic(i):
        cum1[w] += s
cum1


# In[14]:


cum2 = defaultdict(float)
for i in range(5):
    for w, s in combined_lda.show_topic(i):
        cum2[w] += s
cum2


# In[15]:


import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from pprint import pprint

get_ipython().run_line_magic('matplotlib', 'inline')

plot_dictionary = corpora.Dictionary([cum1.keys(), cum2.keys()])
num_words = len(plot_dictionary)


# In[16]:


xs = sorted([k for k in plot_dictionary.token2id])
initial_scores = np.array([cum1[k] for k in xs])
combined_scores = np.array([cum2[k] for k in xs])
difference = combined_scores - initial_scores
difference[difference < 0] = 0


# In[17]:


inds = np.arange(num_words)
plt.bar(inds, initial_scores, align='center')
plt.xticks(inds, xs, rotation=70)


# In[18]:


plt.bar(inds, combined_scores, align='center')
plt.xticks(inds, xs, rotation=70)


# In[19]:


plt.figure(figsize=(16,12))
plt.bar(inds, difference, align='center')
plt.xticks(inds, xs, rotation=70)

