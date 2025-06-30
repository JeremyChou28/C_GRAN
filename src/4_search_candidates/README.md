# Search candidates

you could run the example as follows:

```sh
cd src/4_search_candidates
python process_pubchem_database.py
python search_candidates.py --molecular_network_file test_files/source_target_cor_edit.csv --pubchem_database_path ./pubchem_database.pk --candidates_folder ./candidates/
```