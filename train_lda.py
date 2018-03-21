
# words with count < MIN_COUNTS and count > MAX_COUNTS will be removed
MIN_COUNTS = 20
MAX_COUNTS = 1800

# minimum document length (number of words) after preprocessing
MIN_LENGTH = 15

N_TOPICS = 20
ALPHA = 0.5
ETA = 0.5


# Preprocessing - https://github.com/TropComplique/lda2vec-pytorch/blob/master/utils/preprocess.py

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
                if t.is_alpha and len(t) > 2 and not t.is_stop]

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


import spacy
from gensim import corpora, models


def train_ldas(initial_abstracts, combined_abstracts, n_topics=20, alpha='auto', eta='auto'):
    """Return two LDA models after training on each on the two input abstract lists.
    
    `initial_abstracts` and `combined_abstracts` are lists of tuples containing (PMID, abstract_text).
    """
    nlp = spacy.load('en')

    initial_preprocess = preprocess(
        initial_abstracts, nlp, MIN_LENGTH, MIN_COUNTS, MAX_COUNTS
    )
    combined_preprocess = preprocess(
        combined_abstracts, nlp, MIN_LENGTH, MIN_COUNTS, MAX_COUNTS
    )

    # make single dictionary
    (encoded_docs, decoder, _) = initial_preprocess
    initial_texts = [[decoder[j] for j in doc] for i, doc in encoded_docs]
    (encoded_docs, decoder, _) = combined_preprocess
    combined_texts = [[decoder[j] for j in doc] for i, doc in encoded_docs]

    dictionary = corpora.Dictionary(initial_texts + combined_texts)

    # train initial lda
    corpus = [dictionary.doc2bow(text) for text in initial_texts]
    initial_lda = models.LdaModel(corpus, alpha=alpha, eta=eta, id2word=dictionary, num_topics=n_topics)

    # train combined lda
    corpus = [dictionary.doc2bow(text) for text in combined_texts]
    combined_lda = models.LdaModel(corpus, alpha=alpha, eta=eta, id2word=dictionary, num_topics=n_topics)

    return initial_lda, combined_lda



if __name__ == '__main__':
    print("Training LDA's")
    # Sample usage:
    import pandas as pd
    # Read and dedupe corpora
    initial =  pd.read_csv('initial_corpus.csv').groupby('PMID').first()['AB']
    combined =  pd.read_csv('diabetes_corpus_combined.csv').groupby('PMID').first()['AB']


    # Create list of tuples
    initial_abstracts = zip(initial.index, initial)
    combined_abstracts = zip(combined.index, combined)

    initial_lda, combined_lda = train_ldas(initial_abstracts, combined_abstracts, n_topics=N_TOPICS, alpha=ALPHA, eta=ETA)

    # Display top topics for two different models
    print("Initial corpus topics:")
    for i, topics in initial_lda.show_topics(N_TOPICS, formatted=False):
        print('topic', i, ':', ' '.join([t for t, _ in topics]))

    print("Combined corpus topics:")
    for i in range(N_TOPICS):
       print(initial_lda.show_topic(i))