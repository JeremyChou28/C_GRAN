# C_GRAN

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
This is the project of C_GRAN implemented by Python.

## 💻 Getting Started <a name="getting-started"></a>


### Requirements <a name="requirements"></a>

You should prepare the environment as follows:
```sh
pip install -r requirements.txt
```

### Quick start <a name="quick-start"></a>

1. Calculate correlation

```sh
cd src/1_calculate_correlation

python calculate_correlation.py --input-data test_files/test.txt --compounds-num 14 --samples-num 98 --correlation-result-filename correlation_results.csv
```

2. Filter compounds with high correlation values

```sh
cd src/2_filter_high_correlation_compounds

python filter_high_correlation_compounds.py --correlation_file ../1_calculate_correlation/correlation_results.csv --seednode_file test_files/seednode.csv --threshold 0.7
```

3. Construct molecular network

```sh
cd src/3_construct_molecular_network

python construct_molecular_network.py --source_target_file test_files/source_target.csv --correlation_file ../1_calculate_correlation/correlation_results.csv --threshold 0.7 --RT_threshold 0.01
```


4. Search candidates

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

5. Annotation

- 5.1 Naive annotation

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


- 5.2 Annotation with CFM-ID

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