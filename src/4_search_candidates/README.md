# Search candidates

if you need to prepare pubchem database from scratch, you should run this script first:
```sh
cd src/4_search_candidates

python process_pubchem_database.py
```
or you could download our prepared pubchem database from the [google drive](https://drive.google.com/file/d/17Qmie31AWyOmBy-D1hfhvRoRn2VmdQuC/view?usp=drive_link) or [baiduyun](https://pan.baidu.com/s/1SlKP6dTYZhWI0A5Q3_7thw?pwd=e9dq).

then, run this script for searching candidates:
```sh
python search_candidates.py --edited_molecular_network_file test_files/source_target_cor_edit.csv --pubchem_database_path ./pubchem_database.pk --candidates_folder ./candidates/ --ppm_threshold 2 --is_filter_element --element_set 'C,H,O,N,P,S,F,Cl,Br,I'
```

- `edited_molecular_network_file`: Path to the edited molecular network CSV file with Source, Target, correlation, RT, etc from Step 3.

- `pubchem_database_path`: Path to the preprocessed PubChem database (pickle format).

- `candidates_folder`: Output folder to save retrieved candidate compounds for each node.

- `ppm_threshold`: Mass accuracy threshold in parts per million (ppm) for candidate matching.

- `is_filter_element`: Flag to enable filtering candidates by allowed elements.

- `element_set`: Comma-separated list of allowed chemical elements in candidate formulas.
