# ChildhoodCancerDataInitiative-CCDI_to_SRAy

This repo contains a python script which takes data from a validated CCDI submission manifest and creates an SRA submission file specifically for a CCDI project.

## Table of Contents
- [Python environment management](#python-environment-management)
- [Usage instruction](#usage-instruction)

### Python environment management
A controlled **virtual environment** of Python is always recommanded for running any python package/script due to dependency management purpose. There are many tools that you can use to create a virtual environment, such as `pyenv`, `virtualenv` or `conda`. An instruction is included here on how to create a `conda env` with all the dependencies installed.

- **Conda install**

    Conda is an open source package management system and environment management system that runs on Windows, macOs, and Lunix. [Here](https://docs.conda.io/projects/miniconda/en/latest/) is the site of installation instruction. Please pick the right package based on your operation system.

- **Create a conda env**

    An environment yaml `conda_environment.yml` can be be found under folder `envs/`. To create the environment, simply run

    ```
    conda env create -f <path_to_env_yml>
    ```
    You should be able to find an environment called `CCDI_to_SRA_env` when you run 

    ```
    conda env list
    ```
- **Activate the environment**
  
    All the dependecies that the script requires should be succesfully installed within this environment. To activate the environemnt, simply run

    ```
    conda activate CCDI_to_SRA_env
    ```

    You should be able to see `(CCDI_to_SRA_env)` at the begining of your terminal prompt line after activation.

### Usage instruction

```
>> python CCDI_to_SRA.py --help
usage: CCDI_to_SRA.py [-h] -f MANIFEST -t TEMPLATE [-s PREVIOUS_SUBMISSION]

This script is a python version to generate an SRA submission file using a validated CCDI submission manifest

required arguments:
  -f MANIFEST, --manifest MANIFEST
                        A validated dataset file based on the template CCDI_submission_metadata_template (.xlsx)
  -t TEMPLATE, --template TEMPLATE
                        A dbGaP SRA metadata template, 'phsXXXXXX.xlsx'

optional arguments:
  -s PREVIOUS_SUBMISSION, --previous_submission PREVIOUS_SUBMISSION
                        A previous SRA submission file (xlsx) from the same phs_id study.
```

- **Inputs**

    The script requires a validated `CCDI manifest` and an `SRA submission template`. The previous SRA submission file is optional.

- **Outputs**

    - ***A log file*** named by the basename of CCDI manifest provided
    - (If the script finishes successfully) ***An SRA submission file*** named by the phs accession numnber extracted from CCDI manifest provided.

- **Debug**

    The script performs couple rounds of validation and verification against SRA tempaltes. Once an `Error` is captured, the error and additional info of errored rows will be documented in both terminal and the log file. Below is an example of error log

```
15:37:17 - ERROR - Multiple sample_IDs were found associated with same library_ID:
('plexus_wawled_7',)
15:37:17 - WARNING - 2 rows were removed due to the sample_ID and library_ID issue
15:37:17 - WARNING - Additional info:
╒════════════════════╤═════════════════╕
│ sample_ID          │ library_ID      │
╞════════════════════╪═════════════════╡
│ vinier_sideburns_4 │ plexus_wawled_7 │
├────────────────────┼─────────────────┤
│ outran_horsies_3   │ plexus_wawled_7 │
╘════════════════════╧═════════════════╛
```
