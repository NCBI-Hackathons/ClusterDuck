from dc.utils import get_pubmed_ids_from_phenotypes, get_abstracts_from_pubmed_ids, get_pubmed_ids_from_rsids, get_rsids_from_pubmed_id
test_queries = [
                'Autistic behavior',
                'Restrictive behavior',
                'Impaired social interactions',
                'Poor eye contact',
                'Impaired ability to form peer relationships',
                'No social interaction',
                'Impaired use of nonverbal behaviors',
                'Lack of peer relationships',
                'Stereotypy'
]

pubmed_ids = get_pubmed_ids_from_phenotypes(test_queries)

rsids = []
for id in pubmed_ids:
    rsid = get_rsids_from_pubmed_id(id)
    rsids.extend(rsid)


expand_pubmed_ids = get_pubmed_ids_from_rsids(rsids)

abstracts = get_abstracts_from_pubmed_ids(pubmed_ids)

expand_abstracts = get_abstracts_from_pubmed_ids(expand_pubmed_ids)

with open('original.csv', 'w') as f:
    for p, a in zip(pubmed_ids, abstracts):
        abstract = a.get('AB')
        if abstract is None:
            continue
        else:
            f.write(p + ',')
            if ',' in a['AB']:
                f.write('"' + a['AB'] + '"\n')
            else:
                f.write(a['AB'] + '\n')

with open('expand.csv', 'w') as f:
    for p, a in zip(expand_pubmed_ids, expand_abstracts):
        abstract = a.get('AB')
        if abstract is None:
            continue
        else:
            f.write(p + ',')
            if ',' in a['AB']:
                f.write('"' + a['AB'] + '"\n')
            else:
                f.write(a['AB'] + '\n')
