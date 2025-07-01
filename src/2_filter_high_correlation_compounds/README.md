# Filter compounds with high correlation values

you could run the example as follows:

```sh
cd src/2_filter_high_correlation_compounds
python filter_high_correlation_compounds.py --correlation_file ../1_calculate_correlation/correlation_results.csv --seednode_file test_files/seednode.csv --threshold 0.7
```