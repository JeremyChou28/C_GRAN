# Calculate correlation

you could run the example as follows:

```sh
cd src/1_calculate_correlation

python calculate_correlation.py --intensity_file test_files/test.txt --compounds_num 14 --samples_num 98 --correlation_result_filename correlation_results.csv
```

- `intensity_file`: Path to the input data file. 

- `compounds_num`: Number of compounds in the dataset (i.e., number of rows in the input file).

- `samples_num`: Number of samples per compound (i.e., number of columns in the input file).

- `correlation_result_filename`: Name of the output CSV file that will store the computed correlation coefficients.
