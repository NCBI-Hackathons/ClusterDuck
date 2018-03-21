import six.moves.urllib
from Bio import Entrez, Medline
import re

Entrez.email = "hsiaoyi0504@gmail.com"


def get_pubmed_ids_from_phenotypes(phenotypes):
    pubmed_ids = []
    for phen in phenotypes:
        query = phen
        handle = Entrez.esearch(db='pubmed', term=query, rettype='medline', retmode='text', retmax=200)
        record = Entrez.read(handle)
        handle.close()
        idList = record['IdList']
        pubmed_ids.extend(idList)
    return pubmed_ids


def get_pubmed_ids_from_rsids(rsids):
    pubmed_ids = []
    for id in rsids:
        query = 'rs' + id + 'AND pubmed_snp_cited[sb]'
        handle = Entrez.esearch(db='pubmed', term=query, rettype="medline", retmode="text", retmax=200)
        record = Entrez.read(handle)
        handle.close()
        idList = record['IdList']
        pubmed_ids.extend(idList)
    return pubmed_ids


def get_rsids_from_pubmed_id(pubmed_id):
    result = Entrez.read(Entrez.elink(dbfrom='pubmed', db='snp', linkname='pubmed_snp_cited', id=pubmed_id))
    rsids = []
    if len(result[0]['LinkSetDb']) > 0:
        for r in result[0]['LinkSetDb'][0]['Link']:
            rsids.append(r['Id'])
    return rsids


def get_abstracts_from_pubmed_ids(pubmed_ids):
    abstracts = {}
    handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="medline", retmode=type)
    records = list(Medline.parse(handle))
    handle.close()
    return records
