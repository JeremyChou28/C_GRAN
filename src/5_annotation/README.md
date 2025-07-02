# Annotation

## Naive annotation

you could run the example step by step as follows:

```sh
cd src/5_annotation
python preprocess.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates

python naive_prediction.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --threshold_tanimoto_similarity 0.5
```

or you run the iterative annotation as follows:

```sh
cd src/5_annotation

python run_naive_iterative_annotation.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --threshold_tanimoto_similarity 0.5 --max_iterations 100
```


## Annotation with CFM-ID
First, you should prepare the [CFMID](https://hub.docker.com/r/wishartlab/cfmid) environment, and then you could run the example step by step as follows:

```sh
cd src/5_annotation
python preprocess.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates

python cfmid_prediction.py --num_containers 10 --tolerance 0.1 --energy_level 0 --spectrum_file ./test_files/compounds_spectrum.mgf

python postprocess.py --seednode_file ./test_files/seednode.csv --threshold_modified_cosine_similarity 0.7
```

or you run the iterative annotation as follows:

```sh
cd src/5_annotation

python run_iterative_annotation.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --num_containers 10 --tolerance 0.1 --energy_level 0 --threshold_modified_cosine_similarity 0.5 --spectrum_file ./test_files/compounds_spectrum.mgf --max_iterations 100
```

