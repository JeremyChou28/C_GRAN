# C_GRAN (Community-Guided Recursive Annotation Network)

##  📗 Table of Contents

- [📖 About the Project](#about-project)
- [💻 Getting Started](#getting-started)
  - [Requirements](#requirements)
  - [Quick Start](#quick-start)
  <!-- - [Install](#install) -->
  <!-- - [Usage](#usage) -->
  <!-- - [Run tests](#run-tests) -->
  <!-- - [Deployment](#deployment) -->
<!-- - [🤝 Contributing](#contributing) -->
- [📝 License](#license)
- [👥 Contact](#contact)
- [🔗 Citation](#citation)
<!-- - [🙏 Acknowledgements](#acknowledgements) -->

## 📖 About the Project <a name="about-project"></a>
C-GRAN (Community-Guided Recursive Annotation Network) is an open, network-based framework for the systematic discovery and annotation of emerging structural analogs in complex environmental samples. Designed for non-target screening (NTS) using tandem mass spectrometry (MS/MS), C-GRAN integrates molecular networking with sample-wise co-occurrence analysis to uncover structurally or functionally related compounds beyond spectral similarity constraints. Starting from high-confidence seed annotations, candidate compounds are expanded through a recursive database search strategy, incorporating exact-mass-based matching. Each candidate is ranked based on structural similarity to known analogs, fragment match quality, and occurrence correlation. By iteratively propagating annotations across molecular networks, C-GRAN enables high-coverage identification of structurally diverse compounds—especially those missed by traditional spectral-based tools. 

## 💻 Getting Started <a name="getting-started"></a>


### Requirements <a name="requirements"></a>

You should prepare the environment as follows:
```sh
pip install -r requirements.txt
```

### Quick start <a name="quick-start"></a>

#### Step 1. Calculate correlation

```sh
cd src/1_calculate_correlation

python calculate_correlation.py --input_data test_files/test.txt --compounds_num 14 --samples_num 98 --correlation_result_filename correlation_results.csv
```

- `input_data`: Path to the input data file. 

- `compounds_num`: Number of compounds in the dataset (i.e., number of rows in the input file).

- `samples_num`: Number of samples per compound (i.e., number of columns in the input file).

- `correlation_result_filename`: Name of the output CSV file that will store the computed correlation coefficients.


#### Step 2. Filter compounds with high correlation values

```sh
cd src/2_filter_high_correlation_compounds

python filter_high_correlation_compounds.py --correlation_file ../1_calculate_correlation/correlation_results.csv --seednode_file test_files/seednode.csv --correlation_threshold 0.7
```

- `correlation_file`: Path to the correlation_results from Step 1.

- `seednode_file`: Path to the seed node CSV file. This file should contain a list of initial compounds (including columns such as ID and SMILES) to be used for annotation.

- `correlation_threshold`: Correlation threshold (between 0 and 1).

#### Step 3. Construct molecular network

```sh
cd src/3_construct_molecular_network

python construct_molecular_network.py --source_target_file test_files/source_target.csv --correlation_file ../1_calculate_correlation/correlation_results.csv --correlation_threshold 0.7 --RT_threshold 0.01
```

- `source_target_file`: Path to the molecular network file (CSV), containing columns such as Source, Target, and retention time (RT).

- `correlation_file`: Path to the correlation_results from Step 1.

- `correlation_threshold`: Correlation threshold (between 0 and 1).

- `RT_threshold`: Maximum allowed retention time difference between two nodes to include an edge.

#### Step 4. Search candidates

if you need to prepare pubchem database from scratch, you should run this script first:
```sh
cd src/4_search_candidates

python process_pubchem_database.py
```
or you could download our prepared pubchem database from the [google drive](https://drive.google.com/file/d/17Qmie31AWyOmBy-D1hfhvRoRn2VmdQuC/view?usp=drive_link) or [baiduyun](https://pan.baidu.com/s/1SlKP6dTYZhWI0A5Q3_7thw?pwd=e9dq).

then, run this script for searching candidates:
```sh
python search_candidates.py --molecular_network_file test_files/source_target_cor_edit.csv --pubchem_database_path ./pubchem_database.pk --candidates_folder ./candidates/ --ppm_threshold 2 --is_filter_element --element_set 'C,H,O,N,P,S,F,Cl,Br,I'
```

- `edited_molecular_network_file`: Path to the edited molecular network CSV file with Source, Target, correlation, RT, etc from Step 3.

- `pubchem_database_path`: Path to the preprocessed PubChem database (pickle format).

- `candidates_folder`: Output folder to save retrieved candidate compounds for each node.

- `ppm_threshold`: Mass accuracy threshold in parts per million (ppm) for candidate matching.

- `is_filter_element`: Flag to enable filtering candidates by allowed elements.

- `element_set`: Comma-separated list of allowed chemical elements in candidate formulas.

#### Step 5. Annotation

##### 5.1 Naive annotation

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

##### 5.2 Annotation with CFM-ID

First, you should prepare the [CFMID](https://hub.docker.com/r/wishartlab/cfmid) environment, and then you could run the example step by step as follows:

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



<!-- <p align="right">(<a href="#readme-top">back to top</a>)</p> -->


## 📝 License <a name="license"></a>

This project is licensed under the [MIT License](LICENSE).

<!-- <p align="right">(<a href="#readme-top">back to top</a>)</p> -->



## 👥 Contact <a name="contact"></a>

- Shuping Zheng: 
- Jianping Zhou: jianpingzhou0927@gmail.com


<!-- <p align="right">(<a href="#readme-top">back to top</a>)</p> -->


## 🔗 Citation <a name="citation"></a>

If you find this repo useful, please cite our paper. Thanks for your attention.

<!-- ```
@inproceedings{zhou2024mtsci,
  title={MTSCI: A Conditional Diffusion Model for Multivariate Time Series Consistent Imputation},
  author={Zhou, Jianping and Li, Junhao and Zheng, Guanjie and Wang, Xinbing and Zhou, Chenghu},
  booktitle={Proceedings of the 33rd ACM International Conference on Information and Knowledge Management},
  pages={3474--3483},
  year={2024}
}
``` -->