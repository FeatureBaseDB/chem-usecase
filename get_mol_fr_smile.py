import argparse

import time
from rdkit import Chem
from rdkit.Chem import AllChem
from pilosa import Client, Bitmap, Intersect


def get_mole_id(smile, hosts, db, inverse_db, frame):
    cluster = Client(hosts=[hosts])
    time1 = time.time()
    mol = Chem.MolFromSmiles(smile)
    fp = list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096).GetOnBits())
    bit_maps = [Bitmap(f, frame) for f in fp]
    mole_ids = cluster.query(inverse_db, Intersect(*bit_maps)).values()[0]["bits"]
    existed_mol = False
    found = None
    for m in mole_ids:
        mol = cluster.query(db, Bitmap(m, frame)).values()[0]["bits"]
        if len(mol) == len(fp):
            found = m
            existed_mol = True
            break
    if not existed_mol:
        print "Could not find molecular %s" % smile
        return None
    time2 = time.time()
    print 'finding mol took %0.3f ms' % ((time2 - time1) * 1000.0)
    return found


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest="smile")
    parser.add_argument("-host", dest="host", default="127.0.0.1:10101")
    parser.add_argument("-db", dest="db", default="mol")
    parser.add_argument("-inverse_db", dest="inverse_db", default="inverse-mol")
    parser.add_argument("-f", dest="frame", default="mole.n")
    args = parser.parse_args()
    smile = args.smile
    host = args.host
    db = args.db
    frame = args.frame
    inverse_db = args.inverse_db
    mole_id = get_mole_id(smile, host, db, inverse_db, frame)
    if mole_id:
        print "Molecule ID from %s: %s" % (smile, mole_id)