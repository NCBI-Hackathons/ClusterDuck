from dc.utils import get_pubmed_ids_from_phenotypes, get_pubmed_ids_from_rsids, get_rsids_from_pubmed_id
import argparse

# handle command line arguments
parser = argparse.ArgumentParser(description='Find new association based on phenotypes (HPO terms)')
parser.add_argument('phenotpyes', type=str, nargs='+', help='a list of phenotypes')


if __name__ == '__main__':
    # handle processing logics
    args = parser.parse_args()
    # get initial pubmed ids
    pubmed_ids = get_pubmed_ids_from_phenotypes(args.phenotpyes)
    # expand pubmed ids
    rsids = []
    for id in pubmed_ids:
        rsids.extend(get_rsids_from_pubmed_id(id))
    expand_pubmed_ids = get_pubmed_ids_from_rsids(rsids)
    # TODO: fetch these pubmed documents and do the clustering