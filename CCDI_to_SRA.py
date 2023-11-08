import argparse
from typing import List, TypeVar, Dict
import warnings
import logging
from datetime import date
import logging
import pandas as pd
import sys
import os
from shutil import copy


ExcelReader = TypeVar("ExcelReader")
DataFrame = TypeVar("DataFrame")
ExcelFile = TypeVar("ExcelFile")
Series = TypeVar("Series")


def get_date() -> str:
    """Returns the current date while the script is running"""
    date_obj = date.today()
    return date_obj.isoformat()


def get_logger(loggername: str, log_level: str):
    """Returns a basic logger with a logger name using a std format

    log level can be set using one of the values in log_levels.
    """
    log_levels = {  # sorted level
        "notset": logging.NOTSET,  # 00
        "debug": logging.DEBUG,  # 10
        "info": logging.INFO,  # 20
        "warning": logging.WARNING,  # 30
        "error": logging.ERROR,  # 40
    }
    # logger_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logger_format = "%(asctime)s - %(levelname)s - %(message)s"
    logger_filename = loggername + "_ccdi_to_sra_" + get_date() + ".log"
    logger = logging.getLogger(loggername)
    logger.setLevel(log_levels[log_level])

    # set the stream handler
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(logging.Formatter(logger_format, "%Y-%m-%d %H:%M:%S"))
    stream_handler.setLevel(log_levels["info"])

    # set the file handler
    file_handler = logging.FileHandler(logger_filename, mode="w")
    file_handler.setFormatter(logging.Formatter(logger_format, "%Y-%m-%d %H:%M:%S"))

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    return logger


def sra_template_to_dict(excel_file: ExcelFile) -> Dict:
    """Reads SRA tempalte in Excel format and returns
    a dictionary with sheetnames as keys and pandas
    dataframes as values
    """
    warnings.simplefilter(action="ignore", category=UserWarning)
    sra_dict = {}
    sra_dict["Sequence_Data"] = pd.read_excel(
        excel_file, sheet_name="Sequence_Data", header=0
    )
    sra_dict["Terms"] = pd.read_excel(excel_file, sheet_name="Terms", header=None)
    excel_file.close()
    return sra_dict


def ccdi_manifest_to_dict(excel_file: ExcelFile) -> Dict:
    """Reads a validated CDDI manifest excel and retruns
    a dictionary with sheetnames as keys and pandas
    dataframes as values
    """
    sheets_to_avoid = ["README and INSTRUCTIONS", "Dictionary", "Terms and Value Sets"]
    ccdi_dict_raw = excel_sheets_to_dict(excel_file, no_names=sheets_to_avoid)
    ccdi_dict = {}
    for key, item_df in ccdi_dict_raw.items():
        # drop the column "type" from data frame
        item_df = item_df.drop(["type"], axis=1)
        # remove any line or column that has all na values
        item_df.dropna(axis=0, how="all", inplace=True)
        item_df.dropna(axis=1, how="all", inplace=True)

        # some more filtering criteria
        # test if the df is empty
        # test if all column names contain a '.', if yes, do not add it to dict
        item_df_names = item_df.columns
        if not item_df.empty:
            if len([j for j in item_df_names if "." in j]) != len(item_df_names):
                ccdi_dict[key] = item_df
            else:
                pass
        else:
            pass
    del ccdi_dict_raw
    return ccdi_dict


def excel_sheets_to_dict(excel_file: ExcelFile, no_names: List) -> Dict:
    """Returns a list of sheet names in the excel file input"""
    warnings.simplefilter(action="ignore", category=UserWarning)
    sheetnames = excel_file.sheet_names
    sheetnames_subset = [i for i in sheetnames if i not in no_names]
    excel_dict = {}
    for i in sheetnames_subset:
        i_df = pd.read_excel(excel_file, sheet_name=i, dtype=str)
        excel_dict[i] = i_df
    excel_file.close()
    return excel_dict


def get_acl(workbook_dict: Dict) -> str:
    """Takes a workbook dict and returns acl value"""
    if "acl" in workbook_dict["study"].columns:
        acl_list = workbook_dict["study"]["acl"].tolist()
        # we only expect to find one acl at a time
        acl = acl_list[0].strip("[]'")
    elif "acl" in workbook_dict["study_admin"].columns:
        acl_list = workbook_dict["study_admin"]["acl"].tolist()
        # only expects one acl
        acl = acl_list[0].strip("[]'")
    else:
        acl = ""
    return acl


def get_study_name(workbook_dict: Dict) -> str:
    """Returns study name in CCDI manifest
    Only expects one study name
    """
    df_study_name = workbook_dict["study"]["study_name"].tolist()[0]
    return df_study_name


def remove_redundant_cols(df: DataFrame) -> DataFrame:
    """Removes extended columns of "filetype", "filename", and "MD5_checksum"

    The final data only contains filetype and filetype.1,
    filename and filename.1, MD5_checksum and MD5_checksum.1
    """
    match_list = ["filetype.", "filename.", "MD5_checksum."]
    col_to_drop = []
    df_cols = df.columns.tolist()
    for h in df_cols:
        if any(x in h for x in match_list) and h[len(h) - 1] != "1":
            col_to_drop.append(h)
        else:
            pass
    df = df.drop(col_to_drop, axis=1)
    return df


def sra_match_manifest_seq(
    sra_seq_df: DataFrame, manifest_seq_df: DataFrame
) -> DataFrame:
    """Returns a SRA dataframe and fills the dataframe
    with information of "sequencing data" sheet of CCDI manifest
    """
    sra_seq_df["phs_accession"] = manifest_seq_df["acl"]
    sra_seq_df["sample_ID"] = manifest_seq_df["sample.sample_id"]
    sra_seq_df["library_ID"] = manifest_seq_df["library_id"]
    sra_seq_df["title/short description"] = manifest_seq_df["study_name"]
    sra_seq_df["library_strategy (click for details)"] = manifest_seq_df[
        "library_strategy"
    ]
    sra_seq_df["library_source (click for details)"] = manifest_seq_df["library_source"]
    sra_seq_df["library_selection (click for details)"] = manifest_seq_df[
        "library_selection"
    ]
    sra_seq_df["library_layout"] = manifest_seq_df["library_layout"]
    sra_seq_df["platform (click for details)"] = manifest_seq_df["platform"]
    sra_seq_df["instrument_model"] = manifest_seq_df["instrument_model"]
    sra_seq_df["design_description"] = manifest_seq_df["design_description"]
    sra_seq_df["reference_genome_assembly (or accession)"] = manifest_seq_df[
        "reference_genome_assembly"
    ]
    sra_seq_df["alignment_software"] = manifest_seq_df["sequence_alignment_software"]
    sra_seq_df["filetype"] = manifest_seq_df["file_type"]
    sra_seq_df["filename"] = manifest_seq_df["file_name"]
    sra_seq_df["MD5_checksum"] = manifest_seq_df["md5sum"]
    sra_seq_df["Bases"] = manifest_seq_df["number_of_bp"]
    sra_seq_df["Reads"] = manifest_seq_df["number_of_reads"]
    sra_seq_df["coverage"] = manifest_seq_df["coverage"]
    sra_seq_df["AvgReadLength"] = manifest_seq_df["avg_read_length"]
    sra_seq_df["active_location_URL"] = [
        os.path.dirname(i) for i in manifest_seq_df["file_url_in_cds"].tolist()
    ]
    return sra_seq_df


def fix_design_description(description: List) -> List:
    """Extend the length of design desciption to at least
    250 characters long.
    """
    fixed_list = []
    for k in description:
        if pd.isna(k):
            # if item is nan, add empty str with space
            fixed_list.append("".ljust(250) + ".")
        elif len(k) < 250:
            # if item less than 250, adjust the len to 250
            fixed_list.append(k.ljust(250) + ".")
        else:
            fixed_list.append(k)
    return fixed_list


def get_sra_terms_dict(terms_df: DataFrame) -> Dict:
    """Returns an SRA terms dictionary
    which can be used for verification purpose.

    Keys of sra terms dict: [
    "strategy", "source", "selection", "layout",
    "platform", "model", "type"
    ]
    """
    terms_dict = {}
    # extract strategy df
    strategy_df = terms_df.iloc[1:36, 0:2].reset_index(drop=True)
    strategy_df.rename(columns=strategy_df.iloc[0], inplace=True)
    strategy_df.drop(strategy_df.index[0], inplace=True)
    strategy_df = strategy_df.reset_index(drop=True)
    terms_dict["strategy"] = strategy_df
    # extract source df
    source_df = terms_df.iloc[37:45, 0:2].reset_index(drop=True)
    source_df.rename(columns=source_df.iloc[0], inplace=True)
    source_df.drop(source_df.index[0], inplace=True)
    source_df = source_df.reset_index(drop=True)
    terms_dict["source"] = source_df
    # extract selection df
    selection_df = terms_df.iloc[46:80, 0:2].reset_index(drop=True)
    selection_df.rename(columns=selection_df.iloc[0], inplace=True)
    selection_df.drop(selection_df.index[0], inplace=True)
    selection_df = selection_df.reset_index(drop=True)
    terms_dict["selection"] = selection_df
    # define layout df
    layout_df = pd.DataFrame({"Layout": ["paired", "single"]})
    terms_dict["layout"] = layout_df
    # extract platform df
    platform_df = terms_df.iloc[81:88, 0:1].reset_index(drop=True)
    platform_df.rename(columns=platform_df.iloc[0], inplace=True)
    platform_df.drop(platform_df.index[0], inplace=True)
    platform_df = platform_df.reset_index(drop=True)
    terms_dict["platform"] = platform_df
    # extract platform model df
    model_df = terms_df.iloc[81:103, 1:7].reset_index(drop=True)
    model_df.rename(columns=model_df.iloc[0], inplace=True)
    model_df.drop(model_df.index[0], inplace=True)
    model_df = model_df.reset_index(drop=True)
    terms_dict["model"] = model_df
    # define type df
    type_df = pd.DataFrame(
        {
            "filetype": [
                "bam",
                "fastq",
                "cram",
                "sff",
                "reference_fasta",
                "OxfordNanopore_native",
                "PacBio_HDF5",
                "csv",
                "tab",
                "bai",
                "crai",
                "vcf",
                "bcf",
                "vcf_index",
            ]
        }
    )
    terms_dict["type"] = type_df
    return terms_dict


def reformat_sra_values(sra_df: DataFrame) -> DataFrame:
    """Reformat values of SRA dataframe

    The value modification based on CCDI model v1.7.0

    Fields reworked: [
    library_strategy, platform, library layout,
    library source, library selection, filetype,
    design_description
    ]
    """
    # fix library strategy value
    sra_df["library_strategy (click for details)"][
        sra_df["library_strategy (click for details)"].str.contains(
            "Archer_Fusion", na=False
        )
    ] = "OTHER"
    # fix platform value
    sra_df["platform (click for details)"][
        sra_df["platform (click for details)"].str.contains("Illumina", na=False)
    ] = "ILLUMINA"
    sra_df["platform (click for details)"][
        sra_df["platform (click for details)"] == "Ion Torrent"
    ] = "ION_TORRENT"
    sra_df["platform (click for details)"][
        sra_df["platform (click for details)"] == "LS 454"
    ] = "_LS454"
    sra_df["platform (click for details)"][
        sra_df["platform (click for details)"] == "PacBio SMRT"
    ] = "PACBIO_SMRT"
    sra_df["platform (click for details)"][
        sra_df["platform (click for details)"] == "Oxford Nanopore"
    ] = "OXFORD_NANOPORE"

    # fix library layout value
    sra_df["library_layout"][
        sra_df["library_layout"].str.contains("Single end", na=False)
    ] = "single"
    sra_df["library_layout"][
        sra_df["library_layout"].str.contains("Paired end", na=False)
    ] = "paired"

    # fix library source value
    sra_df["library_source (click for details)"] = sra_df[
        "library_source (click for details)"
    ].str.upper()
    sra_df["library_source (click for details)"][
        sra_df["library_source (click for details)"].str.contains("DNA", na=False)
    ] = "GENOMIC"
    sra_df["library_source (click for details)"][
        sra_df["library_source (click for details)"].str.contains("GENOMIC", na=False)
    ] = "GENOMIC"
    sra_df["library_source (click for details)"][
        sra_df["library_source (click for details)"].str.contains("RNA", na=False)
    ] = "TRANSCRIPTOMIC"
    sra_df["library_source (click for details)"][
        sra_df["library_source (click for details)"].str.contains(
            "TRANSCRIPTOMIC", na=False
        )
    ] = "TRANSCRIPTOMIC"

    # fix library selection value
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "Random"
    ] = "RANDOM"
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "Random PCR"
    ] = "RANDOM PCR"
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "Other"
    ] = "other"
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "Unspecified"
    ] = "unspecified"
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "Repeat Fractionation"
    ] = "repeat fractionation"
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "Size Fractionation"
    ] = "size fractionation"
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "cDNA Oligo dT"
    ] = "cDNA_oligo_dT"
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "cDNA Random Priming"
    ] = "cDNA_randomPriming"
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"]
        == "Padlock Probes Capture Method"
    ] = "padlock probes capture method"

    # fix filetype value and convert all values to lower case
    sra_df["filetype"][sra_df["filetype"].str.contains("tbi", na=False)] = "vcf_index"
    sra_df["filetype"][sra_df["filetype"] == "cram_index"] = "crai"

    # fix design_description, extend description to 250 char long
    sra_df["design_description"] = fix_design_description(
        sra_df["design_description"].tolist()
    )
    return sra_df


def find_new_value_in_col(target_col: Series, ref_col: Series) -> List:
    """Returns a list of index of unmatched values
    between two dataframe cols
    """
    bool_series = target_col.dropna().isin(ref_col.dropna())
    index_false = bool_series.index[~bool_series].to_list()
    return index_false


def sra_value_verification(
    sra_df: DataFrame, sra_terms_dict: Dict, logger
) -> DataFrame:
    """Returns a dataframe that passed verification
    using SRA template terms

    Fields for verification: [
    library strategy, library source,
    library selection, library layout,
    platform, platform model, filetype
    ]

    It reports any column names with missing data, and
    any unknown values found in these fields according
    to SRA template terms.
    """
    logger.info("Begin verification against SRA template terms")
    # check missing information
    required_fields = [
        "library_strategy (click for details)",
        "library_source (click for details)",
        "library_selection (click for details)",
        "library_layout",
        "platform (click for details)",
        "instrument_model",
        "filetype",
    ]
    required_df = sra_df[required_fields]
    cols_missing_info = required_df.columns[required_df.isna().any()].tolist()
    if len(cols_missing_info) > 0:
        logger.warning(
            f"These required fields contain missing info: {*cols_missing_info,}"
        )
    else:
        logger.info(f"All required fields are populated.")

    # check if values are accepted for certain fields
    # check library_strategy
    unknown_library_strategy_index = find_new_value_in_col(
        sra_df["library_strategy (click for details)"],
        sra_terms_dict["strategy"]["Strategy"],
    )
    unknown_library_strategy = sra_df["library_strategy (click for details)"][
        unknown_library_strategy_index
    ].tolist()
    if len(unknown_library_strategy) > 0:
        logger.warning(
            f"The following library strategy values are not accepted: {*unknown_library_strategy,}"
        )
    else:
        logger.info("Library strategy verification PASSED")

    # check library source
    unknown_library_source_index = find_new_value_in_col(
        sra_df["library_source (click for details)"],
        sra_terms_dict["source"]["Source"],
    )
    unknown_library_source = sra_df["library_source (click for details)"][
        unknown_library_source_index
    ].tolist()
    if len(unknown_library_source) > 0:
        logger.warning(
            f"The following library source values are not accepted: {*unknown_library_source,}"
        )
    else:
        logger.info("Library source verification PASSED")

    # check library selection
    unknown_library_selection_index = find_new_value_in_col(
        sra_df["library_selection (click for details)"],
        sra_terms_dict["selection"]["Selection"],
    )
    unknown_library_selection = sra_df["library_selection (click for details)"][
        unknown_library_selection_index
    ].tolist()
    if len(unknown_library_selection) > 0:
        logger.warning(
            f"The following library selection values are not accepted: {*unknown_library_selection,}"
        )
    else:
        logger.info("Library selection verification PASSED")

    # check library layout
    unkown_library_layout_index = find_new_value_in_col(
        sra_df["library_layout"], sra_terms_dict["layout"]["Layout"]
    )
    unknown_library_layout = sra_df["library_layout"][
        unkown_library_layout_index
    ].tolist()
    if len(unknown_library_layout) > 0:
        logger.warning(
            f"The following library layout values are not accepted: {*unknown_library_layout,}"
        )
    else:
        logger.info("Library layout verification PASSED")

    # check file type
    unknown_file_type_index = find_new_value_in_col(
        sra_df["filetype"], sra_terms_dict["type"]["filetype"]
    )
    unknown_file_type = sra_df["filetype"][unknown_file_type_index].tolist()
    if len(unknown_file_type) > 0:
        logger.warning(
            f"The following file type values are not accepted: {*unknown_file_type,}"
        )
    else:
        logger.info("File type verification PASSED")

    # check library platform and instrument model
    # create a subset df using only platform and model column
    unknown_platform = []
    unknown_model = []
    platform_model_df = sra_df[
        ["platform (click for details)", "instrument_model"]
    ].dropna(
        subset=["platform (click for details)"]
    )  # remove anyline with empty platform
    platform_model_dict = dict(
        zip(
            platform_model_df["platform (click for details)"],
            platform_model_df["instrument_model"],
        )
    )
    for i in platform_model_dict.keys():
        if i in sra_terms_dict["platform"]["platforms"].to_list():
            model_i = platform_model_dict[i]
            if model_i not in sra_terms_dict["model"][i].dropna().to_list():
                unknown_model.append(model_i)
            else:
                pass
        else:
            unknown_platform.append(i)
    if len(unknown_platform) > 0:
        logger.warning(
            f"The following platform values are not accepted: {*unknown_platform,}"
        )
    else:
        logger.info("Platform verification PASSED")
    if len(unknown_model) > 0:
        logger.warning(
            f"The following model values are not accepted given the platform value: {*unknown_model,}"
        )
    else:
        logger.info("Model verification PASSED")
    unknown_platform_index = find_new_value_in_col(
        sra_df["platform (click for details)"], sra_terms_dict["platform"]["platforms"]
    )

    # a list of index of any row that contains an unknow value
    unknown_index_list = (
        unknown_library_strategy_index
        + unknown_library_source_index
        + unknown_library_selection_index
        + unkown_library_layout_index
        + unknown_platform_index
        + unknown_file_type_index
    )
    unknown_index_list_uniq = list(set(unknown_index_list))
    if len(unknown_index_list_uniq) > 0:
        logger.warning(
            f"{len(unknown_index_list_uniq)} rows were removed due to unknown values found"
        )
    else:
        logger.info(f"All rows were kept because no unknown values were found")

    sra_df = sra_df.drop(index=unknown_index_list_uniq).reset_index(drop=True)

    return sra_df


def spread_sra_df(sra_df: DataFrame) -> DataFrame:
    """Returns a spreaded dataframe if multiple files
    were found with same library ID

    The final return df only contains unqiue library ID
    for each row
    """
    # extract unique values of library_ID columns
    uniq_library = sra_df["library_ID"].unique()
    # create an empty dataframe
    return_df = pd.DataFrame(columns=sra_df.columns.tolist())
    for i in uniq_library:
        # subset a df for that library_ID and reset index starting from 0
        i_df = sra_df[sra_df["library_ID"] == i].reset_index(drop=True)
        # get the line of i_df
        i_df_row = i_df.shape[0]
        if i_df_row == 1:
            return_df = pd.concat([return_df, i_df], axis=0, ignore_index=True)
        else:
            i_df_firstrow = i_df.loc[[0], :]
            for j in range(1, i_df.shape[0]):
                j_filename = "filename." + str(j)
                j_filetype = "filetype." + str(j)
                j_md5 = "MD5_checksum." + str(j)
                i_df_firstrow[j_filename] = i_df.at[i_df.index[j], "filename"]
                i_df_firstrow[j_filetype] = i_df.at[i_df.index[j], "filetype"]
                i_df_firstrow[j_md5] = i_df.at[i_df.index[j], "MD5_checksum"]
            return_df = pd.concat([return_df, i_df_firstrow], axis=0, ignore_index=True)
    return return_df


def duplicated_library_ID(sra_df: DataFrame, logger) -> None:
    """Report if any library IDs were found in multiple lines
    after combining sra df and previsous sra submission df
    """
    library_id_size = (
        sra_df.groupby(["library_ID"])
        .size()
        .reset_index(name="counts")
        .sort_values("counts", ascending=False)
    )
    duplicated_ids = library_id_size[library_id_size["counts"] > 1][
        "library_ID"
    ].tolist()
    logger.error(
        f"These lirbary IDs have been submitted in previous submission and please fix before submission: {*duplicated_ids,}"
    )


def rename_colnames_output(sra_df: DataFrame) -> DataFrame:
    """Rename some column names in SRA dataframe

    "fileanme.#" -> "filename"
    "filetype.#" -> "filetype"
    "MD5_checksum.#" -> "MD5.checksum"
    """
    col_names = sra_df.columns.tolist()
    col_to_fix = [i for i in col_names if "." in i]
    col_rename = {}
    for i in col_to_fix:
        if "filetype" in i:
            i_tofix = "filetype"
        elif "filename" in i:
            i_tofix = "filename"
        elif "MD" in i:
            i_tofix = "MD5_checksum"
        col_rename[i] = i_tofix
    sra_df = sra_df.rename(columns=col_rename)
    return sra_df


def main():
    # set up arguments for this script
    parser = argparse.ArgumentParser(
        description="This script is a python version to generate an SRA submission file using a validated CCDI submission manifest"
    )
    parser._action_groups.pop()
    required_arg = parser.add_argument_group("required arguments")
    optional_arg = parser.add_argument_group("optional arguments")
    required_arg.add_argument(
        "-m",
        "--manifest",
        type=str,
        required=True,
        help="A validated dataset file  based on the template CCDI_submission_metadata_template (.xlsx)",
    )

    required_arg.add_argument(
        "-t",
        "--template",
        type=str,
        help="A dbGaP SRA metadata template, 'phsXXXXXX.xlsx'",
        required=True,
    )

    optional_arg.add_argument(
        "-s",
        "--previous_submission",
        type=str,
        required=False,
        help="A previous SRA submission file (xlsx) from the same phs_id study.",
    )

    args = parser.parse_args()

    manifest = args.manifest
    template = args.template

    # Initiate a logger object for the script
    manifest_base = os.path.splitext(os.path.basename(manifest))[0]
    logger = get_logger(loggername=manifest_base, log_level="info")

    # Check if manifest and template file can be found
    try:
        manifest_f = pd.ExcelFile(manifest)
        logger.info(f"Checking file {manifest}")
        # create a dict using the CCDI manifest
        workbook_dict = ccdi_manifest_to_dict(manifest_f)
        logger.info(f"Reading the validated CCDI manifest {manifest}")
    except FileNotFoundError as err:
        logger.error(err)
        sys.exit()
    except ValueError as err:
        logger.error(err)
        sys.exit()
    except:
        logger.error(f"Issue occurred while openning file {manifest}")
        sys.exit()

    try:
        template_f = pd.ExcelFile(template)
        logger.info(f"Checking file {template}")
        # Create a dict using the SRA template
        sra_dict = sra_template_to_dict(template_f)
        logger.info(f"Reading the SRA template {template}")
    except FileNotFoundError as err:
        logger.error(err)
        sys.exit()
    except ValueError as err:
        logger.error(err)
        sys.exit()
    except:
        logger.error(f"Issue occurred while openning file {template}")
        sys.exit()

    # Read previous submission if the file was provided
    if args.previous_submission:
        pre_sub = args.previous_submission
        try:
            pre_sub_f = pd.ExcelFile(pre_sub)
            logger.info(f"Reading a previous submission file {pre_sub}")
        except FileNotFoundError as err:
            logger.error(err)
        except ValueError as err:
            logger.error(err)
        except:
            logger.error(f"Issue occurred while openning file {pre_sub}")
    else:
        logger.warning("No previsous submission file was provided.")

    # If there is no seuqencing file record in CCDI manifest, exit execution
    if not "sequencing_file" in workbook_dict.keys():
        logger.info(
            "No seuqneincg files found in CCDI submission file, and no SRA submission file will be generated"
        )
        sys.exit()
    else:
        # create a sequencing df if "sequencing_file" exists in workbook_dict
        sequencing_df = workbook_dict["sequencing_file"]
        logger.info(f"Sequecing file records found in validated CCDI manifest")

    # create a dictionary of contents sra_dict["Terms"]
    sra_terms_dict = get_sra_terms_dict(sra_dict["Terms"])

    # extract study acl and name
    sequencing_df["acl"] = get_acl(workbook_dict)
    sequencing_df["study_name"] = get_study_name(workbook_dict)
    logger.info("Extracted study name and study acl")

    # create data frame with the columns of the SRA template
    sra_df = pd.DataFrame(columns=sra_dict["Sequence_Data"].columns.tolist())
    logger.info("Creating Dataframe using SRA Sequence Data sheet as template")
    # drop redundant columns when the used template has more than two filetype, filename, or MD5_checksum
    sra_df = remove_redundant_cols(df=sra_df)
    # fill the cols in sra_df using sequencing df
    sra_df = sra_match_manifest_seq(sra_seq_df=sra_df, manifest_seq_df=sequencing_df)
    logger.info(
        "Fill the sequence data dataframe with info from CCDI manifest sequencing_file sheet"
    )

    # remove a row if the the filename of that row is empty
    sra_df = sra_df[pd.notna(sra_df["filename"])]

    # special fixes (before verifiation)
    sra_df = reformat_sra_values(sra_df)
    logger.info("Reformat the value in the Seuquence Data Dataframe Done.")

    # verification against template
    sra_df = sra_value_verification(
        sra_df=sra_df, sra_terms_dict=sra_terms_dict, logger=logger
    )
    logger.info(f"Verifying values in the Sequence Data DataFrame Done.")
    # Check if there is any row left after verification.
    # abort the script if no row left
    if sra_df.shape[0] == 0:
        logger.error(
            "No row passed verification step. Please fix the values mentioned above in the manifest and rerun the script"
        )
        sys.exit()
    else:
        pass

    # data frame manipulation, spread sra_df if multiple sequencing files are
    # sharing same library_ID
    sra_df = spread_sra_df(sra_df=sra_df)

    # if previous sra submission provided, concatenated
    if args.previous_submission:
        pre_sub_df = sra_template_to_dict(pre_sub_f)
        pre_sub_df_seq = pre_sub_df["Sequence_Data"]
        if not pre_sub_df_seq.empty:
            sra_df = pd.concat([pre_sub_df_seq, sra_df], axis=0, ignore_index=True)
            # check if multiple lines sharing the same library_ID
            if sra_df.groupby(["library_ID"]).size().max() > 1:
                duplicated_library_ID(sra_df=sra_df, logger=logger)
            else:
                pass
        else:
            pass
    else:
        pass

    # rename any colnames, e.g. filetype.1 -> filetype
    sra_df = rename_colnames_output(sra_df=sra_df)
    logger.info(
        "Rename column names in SRA df, and same names for certain columns are expected."
    )

    # create file for final output
    acl_name = sra_df["phs_accession"].tolist()[0]
    sra_output_path = acl_name + "_" + get_date() + "_SRA_submission.xlsx"
    copy(src=template, dst=sra_output_path)
    # write sra_df into Sequence_Data sheet of output
    logger.info(f"Writing output to path: {sra_output_path}")
    with pd.ExcelWriter(
        sra_output_path, mode="a", engine="openpyxl", if_sheet_exists="overlay"
    ) as writer:
        sra_df.to_excel(writer, sheet_name="Sequence_Data", index=False, header=True)
    logger.info(f"Script finished!")


if __name__ == "__main__":
    main()
