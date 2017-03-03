from __future__ import division
import argparse
import time

import requests
from pilosa import Client
cluster = Client()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-id", dest="mol_id")
    parser.add_argument("-host", dest="host", default="127.0.0.1:15000")
    parser.add_argument("-d", dest="db")
    parser.add_argument("-f", dest="frame")
    args = parser.parse_args()
    id = args.mol_id
    host = args.host
    db = args.db
    frame = args.frame

    for threshold in [50, 70, 75, 80, 85, 90]:
        query_string = 'TopN(Bitmap(id=%s, frame="%s"), frame="%s", n=500000, tanimoto=%s)' % (id, frame, frame, threshold)
        time1 = time.time()
        topn = requests.post("http://%s/query?db=%s" % (host, db), data=query_string)
        time2 = time.time()
        print 'topN took %0.3f ms' % ((time2 - time1) * 1000.0)
        print 'Total topN: %s' % len(topn.json()["results"][0])



