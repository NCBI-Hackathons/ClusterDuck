import urllib2, json
from Bio import Entrez
from Bio import Medline
import re

Entrez.email = "trivneel211@gmail.com"
max_res = 1000
format = 'PubTator'
bioconcept = "Disease" # can be "Gene,Mutation,Disease"
accep_pub_types = ["Journal Article", "Clinical Trial"] #add more if needed
hallmark_queries = ['proliferation receptor',#add fuzzy matching
                    'growth factor',
                    'cell cycle',
                    'contact inhibition',
                    'apoptosis',
                    'necrosis',
                    'autophagy',
                    'senescence',
                    'immortalization',
                    'angiogenesis',
                    'angiogenic factor',
                    'metastasis',
                    'mutation',
                    'DNA repair',
                    'adducts',
                    'DNA damage',
                    'inflammation',
                    'oxidative stress',
                    'warburg effect',
                    'growth',
                    'activation',
                    'immune system']
test_queries = ['Autistic behavior',
                'Restrictive behavior',
                'Impaired social interactions',
                'Poor eye contact',
                'Impaired ability to form peer relationships',
                'No social interaction',
                'Impaired use of nonverbal behaviors',
                'Lack of peer relationships',
                'Stereotypy']


def fetch_medline_records(idlist, type):
    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode=type)
    records = Medline.parse(handle)
    
    records = list(records)
    return records


def pubmed_query(max_res, query=""):
    handle = Entrez.esearch(db="pubmed", term=query, rettype="medline", retmode="text", retmax=max_res)
    record = Entrez.read(handle)
    handle.close()
    idlist = record["IdList"]
    return idlist


def read_url(url, as_str=True):
    if as_str:
        urllib_result = urllib2.urlopen(url)
        res = urllib_result.read()
        return res
    else:
        return urllib2.urlopen(url)

def get_records(query):
    handle = Entrez.esearch(db="pubmed", term=query, rettype="medline", retmode="text", retmax=max_res)
    record = Entrez.read(handle)
    handle.close()
    idList = record['IdList']
    records = fetch_medline_records(idList, "text")#this is a Python list
    return records
#for record in records:
#print(record.get("AB", "?"))
#print("\n")

def get_pmids_from_query(query):
    handle = Entrez.esearch(db="pubmed", term=query, rettype="medline", retmode="text", retmax=max_res)
    record = Entrez.read(handle)
    handle.close()
    idList = record['IdList']
    return list(idList)


def disease_extract(records):# uses PubTator (MEDIC disease dictionary)
    disease_pattern = re.compile(r"Disease\tD\w\w\w\w\w\w")
    for record in records:
        #if record.get("PT", "?") in accep_pub_types:
        pmid = record.get("PMID", "?")
        url_Submit = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + bioconcept + "/" + pmid + "/" + format + "/"
        url_result = urllib2.urlopen(url_Submit)
        res = url_result.read()
        raw_mesh = re.findall(disease_pattern, res)
        cooked_mesh = [mention.replace("Disease\t", "") for mention in raw_mesh]  # this is called a list comprehension
        cooked_mesh = list(set(cooked_mesh))  # only keep unique disease ids
        print(cooked_mesh)

def rsid_extract(records):
    rsid_pattern = re.compile(r"rs\d+")
    rsids = list()
    for record in records:
        pmid = '21844098'#record.get("PMID", "?")
        url_Submit = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + "Mutation" + "/" + pmid + "/" + format + "/"
        url_result = urllib2.urlopen(url_Submit)
        res = url_result.read()
        raw_mesh = re.findall(rsid_pattern, res)
        cooked_mesh = [mention.replace("rs", "") for mention in raw_mesh]
        cooked_mesh = list(set(cooked_mesh))
        rsids.append(cooked_mesh)
    rsids = [item for sublist in rsids for item in sublist]
    return rsids

def get_abstracts(pmids):
    records = fetch_medline_records(pmids, "Text")
    abstracts = list()
    for record in records:
        abstracts.append(record.get("PMID", "?"), record.get("AB", "?"))
    return abstracts

def get_rsids_from_phenotypes(pmids):
    rsid_pattern = re.compile(r"rs\d+")
    rsids = list()
    for pmid in pmids:
        url_Submit = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/" + "Mutation" + "/" + pmid + "/" + format + "/"
        url_result = urllib2.urlopen(url_Submit)
        res = url_result.read()
        raw_mesh = re.findall(rsid_pattern, res)
        cooked_mesh = [mention.replace("rs", "") for mention in raw_mesh]
        cooked_mesh = list(set(cooked_mesh))
        rsids.append(cooked_mesh)
    rsids = [item for sublist in rsids for item in sublist]
    return rsids


def get_pmids(rsids):
    pmids = list()
    for id in rsids:
        query = 'rs' + id + 'AND pubmed_snp_cited[sb]'
        handle = Entrez.esearch(db="pubmed", term=query, rettype="medline", retmode="text", retmax=max_res)
        record = Entrez.read(handle)
        handle.close()
        idList = record['IdList']
        pmids.append(idList)
    pmids = [item for sublist in pmids for item in sublist]
    return pmids

def get_pubmed_ids_from_phenotypes(phenotypes):
    idList = list()
    for phenotype in phenotypes:
        temp_ids = get_pmids_from_query(phenotype)
        idList.append(temp_ids)
    idList = [item for sublist in idList for item in sublist]
    return idList


#Test Code
#rsids = ['7041', '4588']
#get_pmids(rsids)

#for query in test_queries:
#  print(query)
# records = get_records(query)
# rsid_extract(records)
