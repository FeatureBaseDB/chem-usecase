import argparse
import csv
import sys
import mysql.connector
from rdkit import Chem
from rdkit.Chem import AllChem
from pilosa import Client

cluster = Client()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-file", dest="file")
    parser.add_argument("-host", dest="host", default="localhost")
    parser.add_argument("-d", dest="db")
    parser.add_argument("-u", dest="user")
    parser.add_argument("-p", dest="password")
    parser.add_argument("-i", dest="inverse", default=False)
    args = parser.parse_args()
    host = args.host
    db = args.db
    usr = args.user
    password = args.password
    file_name = args.file
    inverse = args.inverse

    if not db or not file_name:
        print "Expect database name and csv file name"
        sys.exit()

    if password is None:
        conn = mysql.connector.connect(host="%s" % host, user="%s" % usr, database="%s" % db)
    else:
        conn = mysql.connector.connect(host="%s" % host, user="%s" % usr, password="%s" % password, database="%s" % db)
    cursor = conn.cursor()
    cursor.execute("SELECT molregno, canonical_smiles FROM compound_structures")
    rows = cursor.fetchall()

    i = 0
    error_bits = []
    with open(file_name, 'wb') as csvfile:
        fieldnames = ['row', 'col']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        for mol in rows:
            bitmap_id = mol[0]
            try:
                m = Chem.MolFromSmiles(mol[1])
                finger_print = list(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=4096).GetOnBits())
                for bit in finger_print:
                    if inverse:
                        writer.writerow({"row": bit, "col": bitmap_id})
                    else:
                        writer.writerow({"row": bitmap_id, "col": bit})
                i += 1
            except Exception as ex:
                error_bits.append(bitmap_id)
                continue
    print "Total error bit: %s " % len(error_bits)