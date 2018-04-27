"""
CLI entrypoint to rnaseqer etl commands.

Accesses the RNAseqer API and issues celldb client commands.

"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import sys
import requests
import Queue
import threading
import urllib
from itertools import izip, imap
import csv
from pmap import pmap

from celldb import client as celldb

FILE_KEY = 'GENES_TPM_COUNTS_FTP_LOCATION'

# transform queue
tq = Queue.Queue()

# load queue
lq = Queue.Queue()

def load(connection, samples, features, values, study_id):
    """

    :param host:
    :return:
    """
    print('loading')
    print(connection)
    print('sample list', samples)
    print('sample features', features[0:10])
    print('sample values:', values[0][0:5])
    print('added cohort {}'.format(study_id))
    celldb.upsert_feature_set(connection, '', features)
    celldb.upsert_cohort(connection, study_id, samples)
    return celldb.upsert_samples(connection, samples, features, values)


def transform(filehandle_study_id):
    """
    Takes a filehandle for an FTP location and transforms it an iterator
    that emits well formed samples to be added to celldb.

    :param filehandle:
    :return:
    """
    filehandle = filehandle_study_id[0]
    study_id = filehandle_study_id[1]
    print('transforming')
    transposed = zip(*csv.reader(filehandle, delimiter=str('\t')))
    print(len(transposed[0]), len(transposed))
    feature_ids = transposed[0][1:]
    sample_ids = [x[0] for x in transposed[1:]]
    values = [x[1:] for x in transposed[1:]]
    return sample_ids, feature_ids, values, study_id


def download(filepath_study_id):
    """
    Transforms an FTP file location into a file handle

    :param filepath:
    :return:
    """
    filepath = filepath_study_id[0]
    print('opening')
    print(filepath)
    study_id = filepath_study_id[1]
    return urllib.urlopen(filepath), study_id

def extract(apipath, organism, offset=0, limit=-1):
    """
    Takes an apipath and organism and downloads the FILE_KEY for each.

    :param apipath:
    :param organism:
    :param offset:
    :param limit:
    :return:
    """
    print('extracting')
    url = "{}/json/getStudiesByOrganism/{}".format(apipath, organism)
    studies_json = requests.get(url).json()
    if limit == -1:
        end = len(studies_json)
    else:
        end = offset + limit
    paths = map(lambda x: x[FILE_KEY], studies_json)[offset:end]
    studies = map(lambda x: x['STUDY_ID'], studies_json)[offset:end]
    print(studies_json)
    print(len(paths))
    return imap(download, zip(paths, studies))


def main(args=None):
    parser = argparse.ArgumentParser(
        description='Add RNASeqer data to celldb')
    parser.add_argument("host", default="localhost", type=str,
        help="The host for a celldb instance.")
    parser.add_argument(
        "apipath", type=str, default="http://www.ebi.ac.uk/fg/rnaseq/api/",
        help="The URL to an RNA-seqer API instance.")
    parser.add_argument(
        "organism", type=str,
        help="homo_sapiens")
    parser.add_argument(
        "--limit", default=-1, type=int,
        help="The number of samples to attempt to upsert starting from the"
             "offset. -1 for no limit")
    parser.add_argument(
        "--offset", default=0, type=int,
        help="The number of items to skip before attempting upsert.")
    parsed = parser.parse_args(args)
    connection = celldb.connect(parsed.host)
    print(parsed)
    handles = extract(
        parsed.apipath, parsed.organism, parsed.offset, parsed.limit)
    loaded = pmap(
        lambda x: load(connection, x[0], x[1], x[2], x[3]), pmap(transform, handles))
    for k, l in enumerate(loaded):
        print("Loaded {}".format(k))
        print(l)


if __name__ == "__main__":
    print('hello')
    main(sys.argv)
