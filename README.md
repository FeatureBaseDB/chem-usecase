# chem-usecase
Chemical Similarity Usecase



##Background

Chemical Similarity usecase uses Tanimoto algorithm to calculate % similarity between 1 molecule with others. 
We need to scan the whole data set to filter closest molecules, which inputs are molecule and threshold.
For example, given a molecule with SMILE as "IC=C1/CCC(C(=O)O1)c2cccc3ccccc23" and threshold = 90%, returns set of molecules that have
similarity >= 90%

##Requirement

* RDKit Python: http://www.rdkit.org/docs/Install.html
* ChemBL Dataset: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/
* Pilosa: https://github.com/pilosa/pilosa



##Import dataset to pilosa
Export data from sdf or mysql to csv file, then import csv to mysql

* Export to csv from sdf file:

    * Export to csv file for database to get fingerprint from chembl_id:
    ```
    python import_from_sdf -p <path_to_file> -file <file_name> 
    ```
    * Export to csv file for database to get chembl_id from fingerprint
    ```
    python import_from_sdf -p <path_to_file> - file <file_name> -i True
    ```
    
* Export to csv from mysql:
    * Export to csv file for database to get fingerprint from chembl_id:
    ```
    python import_from_mysql -h <mysql_host> -u <username> -p <password> -db <mysql_database_name> <csv_file_1>
    ```
    * Export to csv file for database to get chembl_id from fingerprint
    ```
    python import_from_sdf -h <mysql_host> -u <username> -p <password> -db <mysql_database_name> <csv_file_2> True
    ```   
* Import from csv file to pilosa (mol and inverse-mol are default database, mole.n is default frame)
   ```
    pilosactl import -d mol -f mole.n <csv_file_1>
    pilosactl import -d inverse-mol -f mole.n <csv_file_2>
    ```
###Queries from pilosa
*Retrieve molecule_ids that have similarity >= 90%

```
python similar.py -s "I\C=C/1\CCC(C(=O)O1)c2cccc3ccccc23" -t 90

```
    
