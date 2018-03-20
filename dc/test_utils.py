from utils import get_pubmed_ids_from_phenotypes, get_pubmed_ids_from_rsids, get_rsids_from_pubmed_id
import unittest


class test_get_pubmed_ids_from_rsids(unittest.TestCase):
    def test(self):
        rsids = ['7041', '4588']
        pubmed_ids = get_pubmed_ids_from_rsids(rsids)
        self.assertEqual(len(pubmed_ids), 40)
        self.assertIsInstance(pubmed_ids, list)
        for id in pubmed_ids:
            self.assertIsInstance(id, str)


class test_get_rsids_from_pubmed_id(unittest.TestCase):
    def test_valid(self):
        pubmed_id = '25626708'
        rsids = get_rsids_from_pubmed_id(pubmed_id)
        self.assertEqual(len(rsids), 29)
        self.assertIsInstance(rsids, list)
        for id in rsids:
            self.assertIsInstance(id, str)

    def test_no_match(self):
        pubmed_id = '29554659'
        rsids = get_rsids_from_pubmed_id(pubmed_id)
        self.assertEqual(rsids, [])


class test_get_pubmed_ids_from_phenotypes(unittest.TestCase):
    def test(self):
        phenotypes = ['Neoplasm', 'Growth abnormality']  # HP:0002664 and HP:0001507
        pubmed_ids = get_pubmed_ids_from_phenotypes(phenotypes)
        self.assertEqual(len(pubmed_ids), 40)
        self.assertIsInstance(pubmed_ids, list)
        for id in pubmed_ids:
            self.assertIsInstance(id, str)


if __name__ == '__main__':
    unittest.main()
