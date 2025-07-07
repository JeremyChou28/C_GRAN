# Filter compounds with high correlation values

you could run the example as follows:

```sh
cd src/2_filter_high_correlation_compounds

python filter_high_correlation_compounds.py --correlation_file ../1_calculate_correlation/correlation_results.csv --seednode_file test_files/seednode.csv --correlation_threshold 0.7
```

- `correlation_file`: Path to the correlation_results from Step 1.

- `seednode_file`: Path to the seed node CSV file. This file should contain a list of initial compounds (including columns such as ID and SMILES) to be used for annotation.

- `correlation_threshold`: Correlation threshold (between 0 and 1).
