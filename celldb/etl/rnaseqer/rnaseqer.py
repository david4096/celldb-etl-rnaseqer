"""
CLI entrypoint to rnaseqer etl commands.

Accesses the RNAseqer API and issues celldb client commands.

"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
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

def load(connection, samples, features, values):
    """

    :param host:
    :return:
    """
    print('loading')
    print(connection)
    print(samples)
    print(features[0:10])
    print(values[0][0:5])
    return celldb.upsert_samples(connection, samples, features, values)


def transform(filehandle):
    """
    Takes a filehandle for an FTP location and transforms it an iterator
    that emits well formed samples to be added to celldb.

    :param filehandle:
    :return:
    """
    print('transforming')
    transposed = zip(*csv.reader(filehandle, delimiter=str('\t')))
    print(len(transposed[0]), len(transposed))
    feature_ids = transposed[0][1:]
    sample_ids = [x[0] for x in transposed[1:]]
    values = [x[1:] for x in transposed[1:]]
    return sample_ids, feature_ids, values


def download(filepath):
    """
    Transforms an FTP file location into a file handle

    :param filepath:
    :return:
    """
    print('opening')
    print(filepath)
    return urllib.urlopen(filepath)

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
    print(len(paths))
    return imap(download, paths)


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
        lambda x: load(connection, x[0], x[1], x[2]), pmap(transform, handles))
    for k, l in enumerate(loaded):
        print("Loaded {}".format(k))
        print(l)