# Annotation

## Naive annotation

you could run the example step by step as follows:

```sh
cd src/5_annotation

python preprocess.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --top_k 10

python naive_prediction.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --threshold_tanimoto_similarity 0.5
```

or you run the iterative annotation as follows:

```sh
cd src/5_annotation

python run_naive_iterative_annotation.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --threshold_tanimoto_similarity 0.5 --max_iterations 100 --top_k 10
```


## Annotation with CFM-ID
First, you should prepare the [CFMID](https://hub.docker.com/r/wishartlab/cfmid) environment, and then you could run the example step by step as follows:

```sh
cd src/5_annotation

python preprocess.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --top_k 10

python cfmid_prediction.py --num_containers 10 --tolerance 0.1 --energy_level 0 --ion_mode positive --spectrum_file ./test_files/compounds_spectrum.mgf 

python postprocess.py --seednode_file ./test_files/seednode.csv --threshold_modified_cosine_similarity 0.7
```

or you run the iterative annotation as follows:

```sh
cd src/5_annotation

python run_iterative_annotation.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --num_containers 10 --tolerance 0.1 --energy_level 0 --ion_mode positive --threshold_modified_cosine_similarity 0.5 --spectrum_file ./test_files/compounds_spectrum.mgf --max_iterations 100 --top_k 10
```

Finally, you could download the molecular structure images as follows:

```sh
python download_mol_imgs.py --annotation_result_file final_annotation_results.csv --structure_image_folder mol_imgs/
```