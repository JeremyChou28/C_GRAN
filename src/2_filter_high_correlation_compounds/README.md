# Filter compounds with high correlation values

you could run the example as follows:

```sh
cd src/2_filter_high_correlation_compounds
python filter_high_correlation_compounds.py --correlation_file ../1_calculate_correlation/corr_pval_final_CD_sediment_pos_3SD_20240828_true_p0.05.csv --seednode_file test_files/seednode.csv
```