import six.moves.urllib
from Bio import Entrez, Medline
import re

Entrez.email = "hsiaoyi0504@gmail.com"


def get_pubmed_ids_from_phenotypes(phenotypes):
    pubmed_ids = []
    for phen in phenotypes:
        query = phen
        handle = Entrez.esearch(db='pubmed', term=query, rettype='medline', retmode='text')
        record = Entrez.read(handle)
        handle.close()
        idList = record['IdList']
        pubmed_ids.extend(idList)
    return pubmed_ids


def get_pubmed_ids_from_rsids(rsids):
    pubmed_ids = []
    for id in rsids:
        query = 'rs' + id + 'AND pubmed_snp_cited[sb]'
        handle = Entrez.esearch(db='pubmed', term=query, rettype="medline", retmode="text")
        record = Entrez.read(handle)
        handle.close()
        idList = record['IdList']
        pubmed_ids.extend(idList)
    return pubmed_ids


def get_rsids_from_pubmed_id(pubmed_id):
    result = Entrez.read(Entrez.elink(dbfrom='pubmed', db='snp', linkname='pubmed_snp_cited', id=pubmed_id))
    try:
        rsids = [r['Id'] for r in result[0]['LinkSetDb'][0]['Link']]
    except IndexError:
        rsids = []
    return rsids

