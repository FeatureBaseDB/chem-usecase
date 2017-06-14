# chem-usecase
Chemical Similarity Usecase



## Background

Chemical Similarity usecase uses Tanimoto algorithm to calculate % similarity between 1 molecule with others. 
We need to scan the whole data set to filter closest molecules, which inputs are molecule and threshold.
For example, given a molecule with SMILE as "IC=C1/CCC(C(=O)O1)c2cccc3ccccc23" and threshold = 90%, returns set of molecules that have
similarity >= 90%

## Requirements

* RDKit Python: http://www.rdkit.org/docs/Install.html
* ChemBL Dataset: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/
* Pilosa: https://github.com/pilosa/pilosa


## Import dataset to pilosa
Export data from sdf or mysql to csv file, then import csv to pilosa. 
Since chembl_id in SD file always come with CHEMBL, e.g CHEMBL6329, which pilosa hasn’t support string id yet, we will remove CHEMBL and import chembl_id with integer after that.
If you export data from mysql, it uses molregno in compound_structures as unique integer ID for molecule, so you just need to import molregno as chembl_id

* Export to csv from sdf file:
        ```
        python import_from_sdf -p ~/Downloads/chembl_22.sdf -file id_fingerprint.csv 
        ```
* Export to csv from mysql:
        ```
        python import_from_mysql -d chembl_22 -u root -file id_fingerprint.csv
        ```
    
* Import from csv file to pilosa (mol and inverse-mol are default database, mole.n is default frame)
   * Open Pilosa server:
       ```
        pilosa server
        ```
   * Create mole database
   
       ```bash
        curl localhost:10101/index/mole -X POST -d '{"options": {"columnLabel": "position_id"}}'
        ```
  
   * Create frame fingerprint
       ```bash
         curl localhost:10101/index/mole/frame/fingerprint -X POST -d '{"options": {"rowLabel": "chembl_id", "inverseEnabled": true, "cacheSize": 2000000}}'
        ```
   
   * Run import script
        ```bash
         pilosa import -d mole -f fingerprint id_fingerprint.csv
        ```
       
    
    
    
## Queries

* Retrieve molecule_ids that have similarity >= 90%. (Default db="mole", frame="fingerprint", hosts=127.0.0.1:10101)
    ```
    python similar.py -s "I\C=C/1\CCC(C(=O)O1)c2cccc3ccccc23" -t 90
    ```

* Retrive molecule_id from a smile

    ```
    python get_mol_fr_smile.py -s "I\C=C/1\CCC(C(=O)O1)c2cccc3ccccc23"
    ```
    
* Benchmark running for thresholds = [50, 70, 75, 80, 85, 90]
    
     ```
     python benchmarks -id 24
    ```
    
    
## License

chem-usecase by Linh Vo

To the extent possible under law, the person who associated CC0 with pilosa chem-usecase has waived all copyright and related or neighboring rights to chem-usecase.

You should have received a copy of the CC0 legalcode along with this work. If not, see http://creativecommons.org/publicdomain/zero/1.0/.    
