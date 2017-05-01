from __future__ import division

import argparse
import time

import requests
from rdkit import Chem
from rdkit.Chem import AllChem


def get_similarity(smile, threshold, host, db, frame):
    time1 = time.time()
    mol = Chem.MolFromSmiles(smile)
    fp = list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096).GetOnBits())
    bit_maps = ["Bitmap(position_id=%s, frame=%s, inversed=%s)" % (f, frame, True) for f in fp]
    bitmap_string = ', '.join(bit_maps)
    intersection = "Intersect(%s)" % bitmap_string
    mole_ids = requests.post("http://%s/index/%s/query" % (host, db), data=intersection).json()["results"][0]["bits"]
    existed_mol = False
    found = None
    for m in mole_ids:
        mol = requests.post("http://%s/index/%s/query" % (host, db), data="Bitmap(chembl_id=%s, frame=%s)" % (m, frame)).json()["results"][0]["bits"]
        existed_mol = False
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
        query_string = 'TopN(Bitmap(chembl_id=%s, frame="%s"), frame="%s", n=2000000, tanimotoThreshold=%s)' % (found, frame, frame, threshold)
        topn = requests.post("http://%s/index/%s/query" % (host, db), data=query_string)
        if topn.status_code != 200:
            print "Error query similar molecules: %s" % topn.json()
            return []
        return topn.json()["results"][0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest="smile")
    parser.add_argument("-t", dest="threshold", default="90")
    parser.add_argument("-host", dest="host", default="127.0.0.1:10101")
    parser.add_argument("-i", dest="index", default="mole")
    parser.add_argument("-f", dest="frame", default="fingerprint")
    args = parser.parse_args()
    smile = args.smile
    threshold = args.threshold
    host = args.host
    db = args.index
    frame = args.frame
    similar_list = get_similarity(smile, threshold, host, db, frame)
    print similar_list