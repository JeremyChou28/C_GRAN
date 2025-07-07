# Annotation

## Naive annotation

you could run the example step by step as follows:

```sh
cd src/5_annotation

python preprocess.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --top_k 10

python naive_prediction.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --tanimoto_similarity_threshold 0.5
```

- `molecular_network_file`: Path to the edited molecular network CSV file with Source, Target, correlation, RT, etc from Step 3.

- `seednode_file`: Path to the seed node CSV file. This file should contain a list of initial compounds (including columns such as ID and SMILES) to be used for annotation.

- `candidates_folder`: Path to the folder with candidate compounds per node from Step 4.

- `top_k`: Number of top candidates to retain per node based on structural similarity.

- `tanimoto_similarity_threshold`: Minimum Tanimoto similarity to accept candidate annotations.

or you run the iterative annotation as follows:

```sh
cd src/5_annotation

python run_naive_iterative_annotation.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --tanimoto_similarity_threshold 0.5 --max_iterations 100 --top_k 10
```

- `max_iterations`: Maximum number of annotation rounds during the iterative annotation process.

## Annotation with CFM-ID

First, you should prepare the [CFM-ID](https://hub.docker.com/r/wishartlab/cfmid) environment, and then you could run the example step by step as follows:

```sh
cd src/5_annotation

python preprocess.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --top_k 10

python cfmid_prediction.py --num_containers 10 --tolerance 0.1 --energy_level 0 --ion_mode positive --spectrum_file ./test_files/compounds_spectrum.mgf --modified_cosine_similarity_threshold 0.7

python postprocess.py --seednode_file ./test_files/seednode.csv 
```

- `num_containers`: Number of Docker containers to run in parallel for CFM-ID predictions.

- `tolerance`: Mass tolerance for matching predicted and experimental peaks.

- `energy_level`: Collision energy level (e.g., 0, 10, 20).

- `ion_mode`: Ionization mode, either positive or negative.

- `spectrum_file`: Path to the experimental spectrum file in MGF format.

- `modified_cosine_similarity_threshold`: Minimum modified cosine similarity to accept spectrum match.

or you run the iterative annotation as follows:

```sh
cd src/5_annotation

python run_iterative_annotation.py --molecular_network_file ./test_files/source_target_cor_edit.csv --seednode_file ./test_files/seednode.csv --candidates_folder ../4_search_candidates/candidates --num_containers 10 --tolerance 0.1 --energy_level 0 --ion_mode positive --modified_cosine_similarity_threshold 0.5 --spectrum_file ./test_files/compounds_spectrum.mgf --max_iterations 100 --top_k 10
```

- `max_iterations`: Maximum number of annotation rounds during the iterative annotation process.

Finally, you could download the molecular structure images as follows:

```sh
python download_mol_imgs.py --annotation_result_file final_naive_annotation_results.csv --structure_image_folder naive_mol_imgs/
```

- `annotation_result_file`: CSV file containing final annotation results, with SMILES or molecular identifiers.

- `structure_image_folder`: Output folder to save downloaded molecular structure images.
