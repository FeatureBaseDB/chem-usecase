from __future__ import division
import argparse
import time

import requests

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-id", dest="mol_id")
    parser.add_argument("-host", dest="host", default="127.0.0.1:10101")
    parser.add_argument("-i", dest="index", default="mole")
    parser.add_argument("-f", dest="frame", default="fingerprint")
    args = parser.parse_args()
    id = args.mol_id
    host = args.host
    db = args.index
    frame = args.frame

    for threshold in [50, 70, 75, 80, 85, 90]:
        query_string = 'TopN(Bitmap(chembl_id=%s, frame="%s"), frame="%s", n=2000000, tanimotoThreshold=%s)' % (id, frame, frame, threshold)
        time1 = time.time()
        topn = requests.post("http://%s/index/%s/query" % (host, db), data=query_string)
        time2 = time.time()
        print 'topN took %0.3f ms' % ((time2 - time1) * 1000.0)
        print 'Total topN: %s' % len(topn.json()["results"][0])



