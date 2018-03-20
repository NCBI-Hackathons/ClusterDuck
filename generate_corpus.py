from get_abstracts import *

autism = ['Autistic behavior',
          'Restrictive behavior',
          'Impaired social interactions',
          'Poor eye contact',
          'Impaired ability to form peer relationships',
          'No social interaction',
          'Impaired use of nonverbal behaviors',
          'Lack of peer relationships',
          'Stereotypy']

alzheimers = ['Middle age onset',
              'Memory impairment',
              'Progressive forgetfulness',
              'Transient global amnesia',
              'Deficit in phonologic short-term memory',
              'Progressive forgetfulness',
              'Abnormality of the nervous system',
              'Abnormality of metabolism/homeostasis',
              'Alzheimer disease',
              'Sleep disturbance',
              'Sleep-wake cycle disturbance',
              'Autosomal dominant inheritance',
              'Cerebral amyloid angiopathy',
              'Dementia']

diabetes = ['Type I diabetes mellitus',
            'Insulin-resistant diabetes mellitus at puberty',
            'Neonatal insulin-dependent diabetes mellitus',
            'Type I diabetes mellitus',
            'Diabetes mellitus',
            'Maternal diabetes',
            'Maturity-onset diabetes of the young',
            'Type II diabetes mellitus',
            'Insulin-resistant diabetes mellitus',
            'Neonatal insulin-dependent diabetes mellitus',
            'Transient neonatal diabetes mellitus',
            'Diabetic ketoacidosis']

'''autism_pmids = get_pubmed_ids_from_phenotypes(autism)
autism_rsids = get_rsids_from_phenotypes(autism_pmids)
print(autism_rsids)
autism_enlarged = autism_pmids + get_pmids(autism_rsids)
print(autism_enlarged)

alzheimers_pmids = get_pubmed_ids_from_phenotypes(alzheimers)
alzheimers_rsids = get_rsids_from_phenotypes(alzheimers_pmids)
alzheimers_enlarged = alzheimers_pmids + get_pmids(alzheimers_rsids)'''

diabetes_pmids = get_pubmed_ids_from_phenotypes(diabetes)
diabetes_abstracts = get_abstracts(diabetes_pmids)
print(diabetes_abstracts[0])
#diabetes_rsids = get_rsids_from_phenotypes(diabetes_pmids)
#diabetes_enlarged = diabetes_pmids + get_pmids(diabetes_rsids)