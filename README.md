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
* Python-Pilosa: https://github.com/pilosa/python-pilosa


## Import dataset to pilosa
Export data from sdf to csv file, then import csv to mysql
Since chembl_id in SD file always come with CHEMBL, e.g CHEMBL6329, which pilosa hasn’t support string id yet, we will remove CHEMBL and import chembl_id with integer after that.


* Export to csv from sdf file:

    * Export to csv file for database to get fingerprint from chembl_id:
        ```
        python import_from_sdf -p ~/Downloads/chembl_22.sdf -file id_fingerprint.csv 
        ```
    * Export to csv file for database to get chembl_id from fingerprint
        ```
        python import_from_sdf -p ~/Downloads/chembl_22.sdf -file fingerprint_id.csv  -i True
        ```
    
* Import from csv file to pilosa (mol and inverse-mol are default database, mole.n is default frame)
   * Open Pilosa server:
       ```
        pilosa server
        ```
   * Create mol and inverse-mol database
       ```bash
        $ curl -XPOST localhost:10101/db -d '{"db": "mol", "options": {"columnLabel": "position_id"}}'
    
        $ curl -XPOST localhost:10101/db -d '{"db": "inverse-mol", "options": {"columnLabel": "chembl_id"}}'
        ```
  
   * Create frame mole.n for each database
       ```bash
         $ curl -XPOST localhost:10101/frame -d '{"db": "mol", "frame": "mole.n", "options": {"rowLabel": "chembl_id"}}'
        
         $ curl -XPOST localhost:10101/frame -d '{"db": "mol", "frame": "mole.n", "options": {"rowLabel": "position_id"}}'

        ```
   
   * Run import script
        ```bash
         pilosactl import -d mol -f mole.n id_fingerprint.csv
         pilosactl import -d inverse-mol -f mole.n fingerprint_id.csv
        ```
       
    
    
    
## Queries

* Retrieve molecule_ids that have similarity >= 90%. (Default db="mol", inverse-db="inverse-mol", frame="mole.n", hosts=127.0.0.1:10101)
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
    
