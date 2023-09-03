# G-Aligner: a graph-based feature alignment method for untargeted LC-MS-based metabolomics


## Highlights
- **Novelty:** G-Aligner enables comprehensive analysis of all potential correspondences among features all runs for the first time. G-Aligner treats features and potential correspondences as nodes and edges of a multipartite graph, converts the feature matching problem as a multidimensional assignment problem (MAP), and proposes three combinatorial optimization methods to solve the MAP.
- **Accuracy:** G-Aligner achieved the best performance in comparison with popular feature alignment methods in MZmine2, OpenMS and XCMS on two public metabolomics benchmark datasets.
- **Reliability:** G-Aligner achieved the best performance on manually annotated feature lists and untargeted extracted features of MZmine2, OpenMS and XCMS, and helped all compared software obtaining more accurate result by integrating G-Aligner into their workflow.
- **Open source:** We open-sourced G-Aligner under a permissive license to promote the accuracy of MS data analysis more broadly.
- **Dataset:** We manually annotated a feature dataset for both public benchmark datasets, which contains m/z, RT, area information of library analytes and can be used in evaluations of feature detection, quantification and alignment accuracy.

## Datasets
Raw MS files of the metabolomics datasets can be downloaded at [Google Drive](https://drive.google.com/drive/folders/1PRDIvihGFgkmErp2fWe41UR2Qs2VY_5G).

The mzML files of the metabolomics datasets can be downloaded at [Zenodo](https://doi.org/10.5281/zenodo.7995789).

Targeted annotation results, evaluation results and evaluation methods can be downloaded at [Zenodo](https://doi.org/10.5281/zenodo.7995789).


## Setup
1. Prepare the python environment based on your system and hardware.
   
2. Install the dependencies. Here we use ROOT_PATH to represent the root path of G-Aligner.
    
    ```cd ROOT_PATH```
   
    ```pip install -r requirements.txt```



## Run G-Aligner

### Supported formats
Feature extraction rsults in csv format, containing m/z, RT and area columns.

### Demos
Our demos can help you reproduce the evaluation results.

Place the data download from the Zenodo repository as follows.
```
G-Aligner-master
├── data
│   ├── MTBLS562
│   ├── MTBLS562_results_metapro
│   ├── MTBLS562_results_mzmine2
│   ├── MTBLS562_results_openms
│   ├── MTBLS562_results_xcms
│   ├── QE_HF
│   ├── QE_HF_results_metapro
│   ├── QE_HF_results_mzmine2
│   ├── QE_HF_results_openms
│   ├── QE_HF_results_xcms
│   ├── TripleTOF_6600
│   ├── TripleTOF_6600_results_metapro
│   ├── TripleTOF_6600_results_mzmine2
│   ├── TripleTOF_6600_results_openms
│   ├── TripleTOF_6600_results_xcms
│   ├── evaluate_metapro_galigner.py
│   ├── evaluate_metapro_mzmine2.java
│   ├── evaluate_metapro_openms.py
│   ├── evaluate_metapro_xcms.Rmd
│   ├── evaluate_mzmine2_galigner.py
│   ├── evaluate_openms_galigner.py
│   ├── evaluate_xcms_galigner.py
│   ├── metapro_result_comparison.py
│   ├── software_result_comparison.py
```

- To run the benchmark scripts:

```cd ROOT_PATH```

```python data/metapro_result_comparison.py```

```python data/software_result_comparison.py```

- To analyze with G-Aligner:

```cd ROOT_PATH```

Change the parameters in data/evaluate_metapro_galigner.py

```python data/evaluate_metapro_galigner.py```

Feature alignment results are saved in ```experiment``` folder.

## Citation

Cite our paper at:
```
@article{}
```

## License

G-Aligner is an open-source tool, using [***Mulan Permissive Software License，Version 2 (Mulan PSL v2)***](http://license.coscl.org.cn/MulanPSL2)

