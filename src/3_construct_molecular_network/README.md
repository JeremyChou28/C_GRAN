# Construct molecular network

you could run the example as follows:

```sh
cd src/3_construct_molecular_network

python construct_molecular_network.py --molecular_network_file test_files/source_target.csv --correlation_file ../1_calculate_correlation/correlation_results.csv --correlation_threshold 0.7 --RT_threshold 0.01
```

- `molecular_network_file`: Path to the molecular network file (CSV), containing columns such as Source, Target, and retention time (RT).

- `correlation_file`: Path to the correlation_results from Step 1.

- `correlation_threshold`: Correlation threshold (between 0 and 1).

- `RT_threshold`: Maximum allowed retention time difference between two nodes to include an edge.


