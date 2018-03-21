from dc.utils import get_pubmed_ids_from_phenotypes, get_pubmed_ids_from_rsids, get_rsids_from_pubmed_id
import xml.etree.ElementTree as ET
import argparse
from glob import glob
from os.path import dirname, abspath, join

PROJECT_ROOT = dirname(abspath(__file__))

# parser for handling command line arguments
parser = argparse.ArgumentParser(description='Find new association based on phenotypes (HPO terms)')
parser.add_argument('phenotypes', type=str, nargs='+', help='a list of phenotypes')


if __name__ == '__main__':
    # processing arguments
    args = parser.parse_args()
    # get initial pubmed ids
    pubmed_ids = get_pubmed_ids_from_phenotypes(args.phenotypes)
    # expand pubmed ids
    rsids = []
    for id in pubmed_ids:
        rsids.extend(get_rsids_from_pubmed_id(id))
    expand_pubmed_ids = get_pubmed_ids_from_rsids(rsids)
    # load pubmed documents
    abstract = {}
    xml_files = glob(join(PROJECT_ROOT, 'data', '*.xml'))
    for f in xml_files:
        tree = ET.parse('./data/pubmed18n0001.xml')
        root = tree.getroot()
        for article in root.findall('PubmedArticle'):
            try:
                id = article.find('MedlineCitation').find('PMID').text
                abstr_text = article.find('MedlineCitation').find('Article').find('Abstract').find('AbstractText').text
                abstract[id] = abstr_text
            except AttributeError:
                continue
    extracted_abstr = {}
    for id in expand_pubmed_ids:
        abstr_text = abstract.get(id)
        if abstr_text is not None:
            extracted_abstr[id] = abstr_text
    # TODO: do the clustering
    print(extracted_abstr.keys())
