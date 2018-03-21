from six.moves.urllib import request
from os.path import dirname, abspath, exists, join
from os import mkdir
import gzip
import shutil


PROJECT_ROOT = dirname(abspath(__file__))

data_path = join(PROJECT_ROOT, 'data')
if not exists(data_path):
    mkdir(data_path)

for i in range(1, 928):
    num = str(i).zfill(4)
    print('Downloading pubmed18n{}.xml.gz ...'.format(num))
    request.urlretrieve(
        'ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed18n{}.xml.gz'.format(num),
        join(data_path, 'pubmed18n{}.xml.gz'.format(num)))
    print('Extracting pubmed18n{}.xml.gz ...'.format(num))
    with gzip.open(join(data_path, 'pubmed18n{}.xml.gz'.format(num)), 'rb') as f_in:
        with open(join(data_path, 'pubmed18n{}.xml'.format(num)), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
