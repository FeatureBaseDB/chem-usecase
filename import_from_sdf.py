import argparse
import csv
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from pilosa import Client

cluster = Client()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-file", dest="file")
    parser.add_argument("-p", dest="path")
    parser.add_argument("-i", dest="inverse", default=False)
    args = parser.parse_args()
    path = args.path
    file_name = args.file
    inverse = args.inverse

    if not path or not file_name:
        print "Expect path to sdf and csv file name"
        sys.exit()

    suppl = Chem.ForwardSDMolSupplier(path)
    i = 0

    error_bits = []
    with open(file_name, 'wb') as csvfile:
        fieldnames = ['row', 'col']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        for mol in suppl:
            if not mol: continue
            bitmap_id = mol.GetProp("chembl_id")[6:]
            try:
                finger_print = list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096).GetOnBits())
            except Exception as ex:
                error_bits.append(bitmap_id)
                continue
            for bit in finger_print:
                if inverse:
                    writer.writerow({"row": bit, "col": bitmap_id})
                else:
                    writer.writerow({"row": bitmap_id, "col": bit})
            i += 1
    print "Total error bits: %s" % len(error_bits)
