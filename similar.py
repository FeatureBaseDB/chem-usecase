from __future__ import division

import argparse
import sys

import time

import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from pilosa import Client, Bitmap, Intersect

cluster = Client()


def get_similarity(smile, threshold, host, db, inverse_db, frame):
    time1 = time.time()
    mol = Chem.MolFromSmiles(smile)
    fp = list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096).GetOnBits())
    bit_maps = [Bitmap(f, frame) for f in fp]
    mole_ids = cluster.query("inverse-mol", Intersect(*bit_maps)).values()[0]["bits"]
    existed_mol = False
    found = None
    for m in mole_ids:
        mol = cluster.query(db, Bitmap(m, frame)).values()[0]["bits"]
        print mol
        if len(mol) == len(fp):
            found = m
            existed_mol = True
            break
    if not existed_mol:
        print "Could not find molecular %s" % smile
        return
    time2 = time.time()
    print 'finding mol took %0.3f ms' % ((time2 - time1) * 1000.0)

    if found:
        query_string = 'TopN(Bitmap(id=%s, frame="%s"), frame="%s", n=500000, tanimotoThreshold=%s)' % (found, frame, frame, threshold)
        topn = requests.post("http://%s/query?db=%s" % (host, db), data=query_string)
        if topn.status_code != 200:
            print "Error query similar molecules: %s" % topn.json()
            return []
        return topn.json()["results"][0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest="smile")
    parser.add_argument("-t", dest="threshold", default="90")
    parser.add_argument("-host", dest="host", default="127.0.0.1:15000")
    parser.add_argument("-db", dest="db", default="mol")
    parser.add_argument("-inverse-db", dest="inverse_db", default="inverse-mol")
    parser.add_argument("-f", dest="frame", default="mole.n")
    args = parser.parse_args()
    smile = args.smile
    threshold = args.threshold
    host = args.host
    db = args.db
    inverse_db = args.inverse_db
    frame = args.frame
    similar_list = get_similarity(smile, threshold, host, db, inverse_db, frame)
    print similar_list