---
description: Creating a conda environment for running the pipeline
---

# Installation

## 1. Installing conda

Conda installation tutorial can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

## 2. Creating conda environment and installing R/python packages

```bash
conda create --name access_data_analysis python=3
conda activate access_data_analysis
conda install r-essentials r-base r-argparse r-ggpubr r-ggthemes
pip install genotype-variants
```



