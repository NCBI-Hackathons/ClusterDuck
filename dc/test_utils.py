from utils import get_pubmed_ids_from_phenotypes, get_pubmed_ids_from_rsids, get_rsids_from_pubmed_id
import unittest


class test_get_pubmed_ids_from_rsids(unittest.TestCase):
    def test(self):
        rsids = ['7041', '4588']
        self.assertEqual(len(get_pubmed_ids_from_rsids(rsids)), 40)


class test_get_rsids_from_pubmed_id(unittest.TestCase):
    def test(self):
        pubmed_id = '25626708'
        self.assertEqual(len(get_rsids_from_pubmed_id(pubmed_id)), 29)


class test_get_pubmed_ids_from_phenotypes(unittest.TestCase):
    def tests(self):
        phenotypes = ['Neoplasm', 'Growth abnormality']  # HP:0002664 and HP:0001507
        self.assertEqual(len(get_pubmed_ids_from_phenotypes(phenotypes)), 40)


if __name__ == '__main__':
    unittest.main()
