"""
This script is a tool in python that takes validated CCDI manifest
and  SRA templates, and creates a  submission file for SRA.

Authors: Qiong Liu <qiong.liu@nih.gov>
"""
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


class Color:
    """A class for terminal color codes."""

    BOLD = "\033[1m"
    BLUE = "\033[94m"
    WHITE = "\033[97m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD_WHITE = BOLD + WHITE
    BOLD_BLUE = BOLD + BLUE
    BOLD_GREEN = BOLD + GREEN
    BOLD_YELLOW = BOLD + YELLOW
    BOLD_RED = BOLD + RED
    END = "\033[0m"


class ColorLogFormatter(logging.Formatter):
    """A class for formatting colored logs."""

    FORMAT = "%(asctime)s - %(prefix)s%(levelname)s%(suffix)s - %(message)s"

    LOG_LEVEL_COLOR = {
        "DEBUG": {"prefix": "", "suffix": ""},
        "INFO": {"prefix": Color.GREEN, "suffix": Color.END},
        "WARNING": {"prefix": Color.YELLOW, "suffix": Color.END},
        "ERROR": {"prefix": Color.RED, "suffix": Color.END},
        "CRITICAL": {"prefix": Color.BOLD_RED, "suffix": Color.END},
    }

    def format(self, record):
        """Format log records with a default prefix and suffix to terminal color codes that corresponds to the log level name."""
        if not hasattr(record, "prefix"):
            record.prefix = self.LOG_LEVEL_COLOR.get(record.levelname.upper()).get(
                "prefix"
            )

        if not hasattr(record, "suffix"):
            record.suffix = self.LOG_LEVEL_COLOR.get(record.levelname.upper()).get(
                "suffix"
            )

        formatter = logging.Formatter(self.FORMAT, "%H:%M:%S")
        return formatter.format(record)


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

    logger_filename = loggername + "_ccdi_to_sra_" + get_date() + ".log"
    logger = logging.getLogger(loggername)
    logger.setLevel(log_levels[log_level])

    # set the stream handler
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(ColorLogFormatter())
    stream_handler.setLevel(log_levels["info"])

    # set the file handler
    file_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
    file_handler = logging.FileHandler(logger_filename, mode="w")
    file_handler.setFormatter(logging.Formatter(file_FORMAT, "%H:%M:%S"))

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

    The sheet will be dropped if found empty
    """
    sheets_to_avoid = ["README and INSTRUCTIONS", "Dictionary", "Terms and Value Sets"]
    ccdi_dict_raw = excel_sheets_to_dict(excel_file, no_names=sheets_to_avoid)
    ccdi_dict = {}
    for key, item_df in ccdi_dict_raw.items():
        # drop the column "type" from data frame
        item_df = item_df.drop(["type"], axis=1)
        # remove any line or column that has all na values
        item_df.dropna(axis=0, how="all", inplace=True)
        # keep empty columnsat this step
        # item_df.dropna(axis=1, how="all", inplace=True)

        # some more filtering criteria
        # test if the df is empty
        # test if all column names contain a '.', if yes, do not add it to dict
        item_df_names = item_df.columns
        if len([j for j in item_df_names if "." in j]) != len(item_df_names):
            ccdi_dict[key] = item_df
        else:
            pass
    del ccdi_dict_raw
    return ccdi_dict


def concat_seq_single_seq(seq_df: DataFrame, single_df: DataFrame) -> DataFrame:
    """Returns a dataframe that combines sequencing_file
    and single_cell_sequencing_file sheets
    """
    cols_to_keep = [
        "sample.sample_id",
        "pdx.pdx_id",
        "cell_line.cell_line_id",
        "library_id",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "platform",
        "instrument_model",
        "design_description",
        "reference_genome_assembly",
        "sequence_alignment_software",
        "file_type",
        "file_name",
        "md5sum",
        "number_of_bp",
        "number_of_reads",
        "coverage",
        "avg_read_length",
        "file_url_in_cds",
    ]
    seq_df_subset = seq_df[cols_to_keep]
    single_df_subset = single_df[cols_to_keep]
    combined_df = pd.concat(
        [seq_df_subset, single_df_subset], axis=0, ignore_index=True
    )
    # create a sample id list that takes sample.sample_id value, if empty,
    # takes pdx.pdx_id value, if empty, takes cell_line.cell_line_id value
    sample_ID = []
    sample_id_df = combined_df[
        ["sample.sample_id", "pdx.pdx_id", "cell_line.cell_line_id"]
    ]
    for i in range(sample_id_df.shape[0]):
        i_row = sample_id_df.loc[i, :]
        if not pd.isna(i_row["sample.sample_id"]):
            sample_ID.append(i_row["sample.sample_id"])
        else:
            if not pd.isna(i_row["pdx.pdx_id"]):
                sample_ID.append(i_row["pdx.pdx_id"])
            else:
                sample_ID.append(i_row["cell_line.cell_line_id"])
    combined_df["sample_ID"] = sample_ID
    return combined_df


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
    sra_seq_df["sample_ID"] = manifest_seq_df["sample_ID"]
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
        os.path.dirname(i) + "/" for i in manifest_seq_df["file_url_in_cds"].tolist()
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
                "bai",
                "cram",
                "crai",
                "fastq",
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

    # fix instrument value
    sra_df["instrument_model"][
        sra_df["instrument_model"] == "Illumina HiSeq X Ten"
    ] = "HiSeq X Ten"
    sra_df["instrument_model"][
        sra_df["instrument_model"] == "Illumina HiSeq X Five"
    ] = "HiSeq X Five"
    sra_df["instrument_model"][
        sra_df["instrument_model"] == "Illumina NovaSeq"
    ] = "Illumina NovaSeq 6000"

    # fix library layout value
    sra_df["library_layout"][
        sra_df["library_layout"].str.contains("Single end", na=False)
    ] = "single"
    sra_df["library_layout"][
        sra_df["library_layout"].str.contains("single-end", na=False)
    ] = "single"
    sra_df["library_layout"][
        sra_df["library_layout"].str.contains("Paired end", na=False)
    ] = "paired"
    sra_df["library_layout"][
        sra_df["library_layout"].str.contains("paired-end", na=False)
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
    sra_df["library_selection (click for details)"][
        sra_df["library_selection (click for details)"] == "rRNA Depletion"
    ] = "Inverse rRNA"

    # fix filetype value and convert all values to lower case
    sra_df["filetype"][sra_df["filetype"].str.contains("tbi", na=False)] = "vcf_index"
    sra_df["filetype"][sra_df["filetype"] == "cram_index"] = "crai"
    sra_df["filetype"][sra_df["filetype"] == "bam_index"] = "bai"

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


def get_unknown_model_index(sra_df: DataFrame, unknown_model: List) -> List:
    index_to_remove = []
    for i in unknown_model:
        i_platform = i["platform"]
        i_model = i["model"]
        for j in sra_df.index:
            j_platform = sra_df.iloc[j]["platform (click for details)"]
            j_model = sra_df.iloc[j]["instrument_model"]
            if j_platform == i_platform and j_model == i_model:
                index_to_remove.append(j)
            else:
                pass
    return index_to_remove


def sra_value_verification(
    sra_df: DataFrame, sra_terms_dict: Dict, logger
) -> DataFrame:
    """Returns a dataframe that passed verification
    using SRA template terms

    Fields for verification: [
    sample_ID, library_ID,
    library strategy, library source,
    library selection, library layout,
    platform, platform model, filetype
    filename
    ]

    It reports any column names with missing data, and
    any unknown values found in these fields according
    to SRA template terms. Any row containing unknown
    values will be removed before return result
    """
    logger.info("Begin verification against SRA template terms")
    # check missing information
    required_fields = [
        "sample_ID",
        "library_ID",
        "library_strategy (click for details)",
        "library_source (click for details)",
        "library_selection (click for details)",
        "library_layout",
        "platform (click for details)",
        "instrument_model",
        "filetype",
        "filename",
    ]
    required_df = sra_df[required_fields]
    cols_missing_info = required_df.columns[required_df.isna().any()].tolist()
    missing_info_index = []
    if len(cols_missing_info) > 0:
        df_w_missing = required_df[required_df.isna().any(axis=1)]
        logger.error(
            f"{df_w_missing.shape[0]} rows were found with missing value in required fields of\n{*cols_missing_info,}"
        )
        logger.warning(
            "Additional info: sample_ID and library_ID of affected rows\n"
            + df_w_missing[["sample_ID", "library_ID"]].to_markdown(
                tablefmt="fancy_grid", index=False
            )
        )
        missing_info_index = df_w_missing.index.tolist()
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
        logger.error(
            f"The following library strategy values are not accepted: {*list(set(unknown_library_strategy)),}"
        )
        logger.warning(
            "Additional info:\n"
            + sra_df.iloc[unknown_library_strategy_index][
                ["sample_ID", "library_ID", "library_strategy (click for details)"]
            ].to_markdown(
                tablefmt="fancy_grid",
                index=False,
                headers=["sample_ID", "library_ID", "library_strategy"],
            )
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
        logger.error(
            f"The following library source values are not accepted: {*list(set(unknown_library_source)),}"
        )
        logger.warning(
            "Additional info:\n"
            + sra_df.iloc[unknown_library_source_index][
                ["sample_ID", "library_ID", "library_source (click for details)"]
            ].to_markdown(
                tablefmt="fancy_grid",
                index=False,
                headers=[
                    "sample_ID",
                    "library_ID",
                    "library_source",
                ],
            )
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
        logger.error(
            f"The following library selection values are not accepted: {*list(set(unknown_library_selection)),}"
        )
        logger.warning(
            "Additional info:\n"
            + sra_df.iloc[unknown_library_selection_index][
                ["sample_ID", "library_ID", "library_selection (click for details)"]
            ].to_markdown(
                tablefmt="fancy_grid",
                index=False,
                headers=[
                    "sample_ID",
                    "library_ID",
                    "library_selection",
                ],
            )
        )
    else:
        logger.info("Library selection verification PASSED")

    # check library layout
    unknown_library_layout_index = find_new_value_in_col(
        sra_df["library_layout"], sra_terms_dict["layout"]["Layout"]
    )
    unknown_library_layout = sra_df["library_layout"][
        unknown_library_layout_index
    ].tolist()
    if len(unknown_library_layout) > 0:
        logger.error(
            f"The following library layout values are not accepted: {*list(set(unknown_library_layout)),}"
        )
        logger.warning(
            "Additional info:\n"
            + sra_df.iloc[unknown_library_layout_index][
                ["sample_ID", "library_ID", "library_layout"]
            ].to_markdown(tablefmt="fancy_grid", index=False)
        )
    else:
        logger.info("Library layout verification PASSED")

    # check file type
    unknown_file_type_index = find_new_value_in_col(
        sra_df["filetype"], sra_terms_dict["type"]["filetype"]
    )
    unknown_file_type = sra_df["filetype"][unknown_file_type_index].tolist()
    if len(unknown_file_type) > 0:
        logger.error(
            f"The following file type values are not accepted: {*list(set(unknown_file_type)),}"
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
    for i in platform_model_df.index:
        i_platform = platform_model_df.iloc[i]["platform (click for details)"]
        i_model = platform_model_df.iloc[i]["instrument_model"]
        i_platform_model_dict = {"platform": i_platform, "model": i_model}
        if i_platform in sra_terms_dict["platform"]["platforms"].to_list():
            if i_model not in sra_terms_dict["model"][i_platform].dropna().to_list():
                unknown_model.append(i_platform_model_dict)
            else:
                pass
        else:
            unknown_platform.append(i_platform)

    unknown_platform_index = find_new_value_in_col(
        sra_df["platform (click for details)"], sra_terms_dict["platform"]["platforms"]
    )
    unknown_model_index = get_unknown_model_index(sra_df, unknown_model=unknown_model)

    if len(unknown_platform) > 0:
        logger.error(
            f"The following platform values are not accepted: {*list(set(unknown_platform)),}"
        )
        logger.warning(
            "Additional info:\n"
            + sra_df.iloc[unknown_platform_index][
                ["sample_ID", "library_ID", "platform (click for details)"]
            ].to_markdown(
                tablefmt="fancy_grid",
                index=False,
                headers=["sample_ID", "library_ID", "platform"],
            )
        )
    else:
        logger.info("Platform verification PASSED")
    if len(unknown_model) > 0:
        logger.error(
            f"The following model values are not accepted given the platform value: {*list(set([str(k)for k in unknown_model])),}"
        )
        logger.warning(
            "Additional info:\n"
            + sra_df.iloc[unknown_model_index][
                [
                    "sample_ID",
                    "library_ID",
                    "platform (click for details)",
                    "instrument_model",
                ]
            ].to_markdown(
                tablefmt="fancy_grid",
                index=False,
                headers=["sample_ID", "library_ID", "platform", "model"],
            )
        )
    else:
        logger.info("Model verification PASSED")

    # Create a list of index of any row that contains an unknow value
    index_to_remove_verification = (
        missing_info_index
        + unknown_library_strategy_index
        + unknown_library_source_index
        + unknown_library_selection_index
        + unknown_library_layout_index
        + unknown_platform_index
        + unknown_file_type_index
        + unknown_model_index
    )

    index_to_remove_verification_uniq = list(set(index_to_remove_verification))
    logger.info(f"Verifying values in the Sequence Data DataFrame Done.")

    return index_to_remove_verification_uniq


def check_value_constancy(library_id: str, unit_df: DataFrame, logger) -> None:
    """This function checks if certain values are same between files sharing
    the same library_ID and report library_IDs if this issue is found
    """
    check_fields = [
        "library_strategy (click for details)",
        "library_source (click for details)",
        "library_selection (click for details)",
        "library_layout",
        "platform (click for details)",
        "instrument_model",
        "design_description",
        "reference_genome_assembly (or accession)",
        "alignment_software",
        "Bases",
        "Reads",
        "coverage",
        "AvgReadLength",
        "active_location_URL",
    ]
    cols_to_report = []
    for i in check_fields:
        i_uniq_list = unit_df[i].dropna().unique().tolist()
        if len(i_uniq_list) == 0:
            pass
        elif len(i_uniq_list) > 1:
            cols_to_report.append(i)
        else:
            pass
    if len(cols_to_report) > 0:
        logger.warning(
            f"Files from library {library_id} have different values for {*cols_to_report,}"
        )
        logger.warning(
            f"Additional info:\n"
            + unit_df[["library_ID"] + cols_to_report].to_markdown(
                tablefmt="fancy_grid",
                index=False,
            )
        )
    else:
        pass


def spread_sra_df(sra_df: DataFrame, logger) -> DataFrame:
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
        i_df = sra_df[sra_df["library_ID"] == i]
        # reoder i_df based on filetype, prioritize bam and cram files
        i_df = i_df.sort_values(by=["filetype"], ascending=False).reset_index(drop=True)
        # get the line of i_df
        i_df_row = i_df.shape[0]
        if i_df_row == 1:
            return_df = pd.concat([return_df, i_df], axis=0, ignore_index=True)
        else:
            # check value constancy of i_df and report cols with issue
            check_value_constancy(library_id=i, unit_df=i_df, logger=logger)
            i_df_firstrow = i_df.loc[[0], :]
            for j in range(1, i_df.shape[0]):
                j_filename = "filename." + str(j)
                j_filetype = "filetype." + str(j)
                j_md5 = "MD5_checksum." + str(j)
                i_df_firstrow[j_filename] = i_df.at[i_df.index[j], "filename"]
                i_df_firstrow[j_filetype] = i_df.at[i_df.index[j], "filetype"]
                i_df_firstrow[j_md5] = i_df.at[i_df.index[j], "MD5_checksum"]
            return_df = pd.concat([return_df, i_df_firstrow], axis=0, ignore_index=True)
    # check missing values of the entire df
    report_missing_library_id = return_df[return_df.isna().any(axis=1)][
        "library_ID"
    ].tolist()
    if len(report_missing_library_id) > 0:
        logger.warning(
            f"These {len(report_missing_library_id)} library_IDs contains missing information in the final output and please address this issue accordingly:\n"
            + f"{*report_missing_library_id,}"
        )
    else:
        pass
    return return_df


def check_and_remove_duplicates(sra_df: DataFrame, logger) -> None:
    """Report if any filenames were found in multiple lines
    after combining sra df and previsous sra submission df

    If one row shares same library_ID and filename,
    then it is considered same sequencing record
    """
    filename_size = (
        sra_df.groupby(["library_ID", "filename"])
        .size()
        .reset_index(name="counts")
        .sort_values("counts", ascending=False)
    )
    duplicated_filename = filename_size[filename_size["counts"] > 1][
        "filename"
    ].tolist()
    logger.warning(
        f"These filenames have been submitted in previous submission and will be removed: {*duplicated_filename,}"
    )
    # remove duplicates and keep the last record
    # due to the way of concatenation, we keep the last record of
    # duplicate lines, which should be the record in ccdi manifest
    sra_df = sra_df.drop_duplicates(
        subset=["library_ID", "filename"], keep="last"
    ).reset_index(drop=True)
    return sra_df


def reorder_col_names(col_list: List) -> List:
    extended_file = [i for i in col_list if "." in i]
    extended_file_to_add = [
        k
        for k in extended_file
        if k not in ["filetype.1", "filename.1", "MD5_checksum.1"]
    ]
    original_cols = [
        "phs_accession",
        "sample_ID",
        "library_ID",
        "title/short description",
        "library_strategy (click for details)",
        "library_source (click for details)",
        "library_selection (click for details)",
        "library_layout",
        "platform (click for details)",
        "instrument_model",
        "design_description",
        "reference_genome_assembly (or accession)",
        "alignment_software",
        "filetype",
        "filename",
        "MD5_checksum",
        "filetype.1",
        "filename.1",
        "MD5_checksum.1",
        "active_location_URL",
        "Bases",
        "Reads",
        "coverage",
        "AvgReadLength",
    ]
    reordered_cols = original_cols + extended_file_to_add
    return reordered_cols


def rename_colnames_output(sra_df: DataFrame) -> DataFrame:
    """Rename some column names in SRA dataframe

    "fileanme.#" -> "filename"
    "filetype.#" -> "filetype"
    "MD5_checksum.#" -> "MD5.checksum"

    Reorder the columns of the dataframe

    Change the datatype of of ["Bases","Reads","coverage","AvgReadLength"] to numeric
    """
    cols_to_change_type = ["Bases", "Reads", "coverage", "AvgReadLength"]
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
    reordered_colnames = reorder_col_names(col_names)
    sra_df = sra_df[reordered_colnames]
    sra_df = sra_df.rename(columns=col_rename)
    # Change the few cols into numeric datatype
    sra_df[cols_to_change_type] = sra_df[cols_to_change_type].apply(
        pd.to_numeric, errors="coerce"
    )
    return sra_df


def reformat_previous_sra(p_sra_df: DataFrame) -> DataFrame:
    """Reformats the dataframe of SRA Sequence Data dataframe

    For lines with multiple files derived from same library,
    move the extended file info (filetype.1, filename.1, and MD5_checksum.1)
    into a new single line
    """
    sra_cols = p_sra_df.columns.tolist()
    additional_cols = [i for i in sra_cols if "." in i]
    extra_max = max([int(i[len(i) - 1]) for i in additional_cols])
    p_sra_df_nochange = p_sra_df[[i for i in sra_cols if i not in additional_cols]]
    nofile_cols = [
        i
        for i in p_sra_df_nochange.columns
        if i not in ["filetype", "filename", "MD5_checksum"]
    ]

    for i in range(1, extra_max + 1):
        i_filetype_col = "filetype." + str(i)
        i_filename_col = "filename." + str(i)
        i_md5_col = "MD5_checksum." + str(i)
        i_cols = [i_filetype_col, i_filename_col, i_md5_col]
        i_file_df = p_sra_df[pd.notna(p_sra_df[i_filetype_col])][i_cols]
        i_nofile_df = p_sra_df[pd.notna(p_sra_df[i_filetype_col])][nofile_cols]
        if i_file_df.empty:
            pass
        else:
            i_concat = pd.concat([i_nofile_df, i_file_df], axis=1)
            i_concat = i_concat.rename(
                columns={
                    i_filetype_col: "filetype",
                    i_filename_col: "filename",
                    i_md5_col: "MD5_checksum",
                }
            )
            p_sra_df_nochange = pd.concat(
                [p_sra_df_nochange, i_concat], axis=0, ignore_index=True
            ).reset_index(drop=True)

    return p_sra_df_nochange


def validate_sample_library(sra_df: DataFrame, logger) -> List:
    index_rows_to_remove = []
    logger.info(
        "Checking if any library_ID was found associated with multiple sample_IDs"
    )
    sample_library_size = (
        sra_df.groupby(["library_ID", "sample_ID"]).size().groupby("library_ID").size()
    )
    if sample_library_size.max() > 1:
        library_to_fix = sample_library_size[sample_library_size > 1].index.tolist()
        rows_to_remove = sum(sra_df["library_ID"].isin(library_to_fix))
        logger.error(
            f"{rows_to_remove} rows were found with one to many issue between library_ID and sample_ID\nlibrary_IDs involved\n{*library_to_fix,}"
        )
        df_rows_to_remove = sra_df[["sample_ID", "library_ID"]][
            sra_df["library_ID"].isin(library_to_fix)
        ]
        logger.warning(
            "Additional info: sample_ID and library_ID of affected rows\n"
            + df_rows_to_remove.to_markdown(tablefmt="fancy_grid", index=False)
        )
        # get index of rows to be removed
        index_rows_to_remove = df_rows_to_remove.index.tolist()
    else:
        logger.info("sample_ID and library_ID checking PASSED ")
    return index_rows_to_remove


def main():
    # set up arguments for this script
    parser = argparse.ArgumentParser(
        description="This script is a python version to generate an SRA submission file using a validated CCDI submission manifest"
    )
    parser._action_groups.pop()
    required_arg = parser.add_argument_group("required arguments")
    optional_arg = parser.add_argument_group("optional arguments")
    required_arg.add_argument(
        "-f",
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
    except BaseException as err:
        logger.error(f"Issue occurred while openning file {manifest}: {err}")
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
    except BaseException as err:
        logger.error(f"Issue occurred while openning file {template}: {err}")
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
        except BaseException as err:
            logger.error(f"Issue occurred while openning file {pre_sub}: {err}")
    else:
        logger.warning("No previsous submission file was provided.")

    # create a dictionary of contents sra_dict["Terms"]
    sra_terms_dict = get_sra_terms_dict(sra_dict["Terms"])

    # If there is no seuqencing record in CCDI manifest, exit execution
    sequencing_file_df = workbook_dict["sequencing_file"]
    single_sequencing_file_df = workbook_dict["single_cell_sequencing_file"]
    if sequencing_file_df.empty and single_sequencing_file_df.empty:
        logger.info(
            "No seuqneincg file or single cell sequencing file found in CCDI submission file, and no SRA submission file will be generated"
        )
        sys.exit()
    else:
        pass

    # Combine records from sheets sequencing_file and single_cell_sequencing_file
    sequencing_df = concat_seq_single_seq(
        seq_df=sequencing_file_df, single_df=single_sequencing_file_df
    )
    logger.info(
        f"A total of {sequencing_df.shape[0]} sequencing files were found in the provided manifest"
    )

    # extract study acl and name
    sequencing_df["acl"] = get_acl(workbook_dict)
    sequencing_df["study_name"] = get_study_name(workbook_dict)
    logger.info("Extracted study name and study acl")

    # create data frame with the columns of the SRA template
    sra_df = pd.DataFrame(columns=sra_dict["Sequence_Data"].columns.tolist())
    logger.info(
        "Created a Dataframe using SRA Sequence Data sheet as template and CCDI manifest as the source of information"
    )
    # drop redundant columns when the used template has more than two filetype, filename, or MD5_checksum
    sra_df = remove_redundant_cols(df=sra_df)
    # fill the cols in sra_df using sequencing df
    sra_df = sra_match_manifest_seq(sra_seq_df=sra_df, manifest_seq_df=sequencing_df)

    # special fixes (before verifiation)
    sra_df = reformat_sra_values(sra_df)
    logger.info("Reformatting the value in the Seuquence Data Dataframe Done.")

    # check sample_ID and library_ID is one to many relationship
    # Many library_ID can be derived from same sample_ID
    sample_library_issue_index = validate_sample_library(sra_df=sra_df, logger=logger)

    # verification against template and check if any required fields are empty
    verification_issue_index = sra_value_verification(
        sra_df=sra_df, sra_terms_dict=sra_terms_dict, logger=logger
    )

    # combine the index identified from sample_ID library_ID validation
    # and value verification steps
    rows_to_remove_index = list(
        set(sample_library_issue_index + verification_issue_index)
    )
    # remove errored rows identified with smaple_ID issue, missing info,
    # or invalid values
    if len(rows_to_remove_index) > 0:
        sra_df = sra_df.drop(index=rows_to_remove_index).reset_index(drop=True)
        logger.warning(
            f"{len(rows_to_remove_index)} rows were removed due to the errors described above"
        )
    else:
        pass

    # Check if there is any row left after verification.
    # abort the script if no row left
    if sra_df.shape[0] == 0:
        logger.error(
            "No row passed verification step. Please fix the values mentioned above in the manifest and rerun the script"
        )
        sys.exit()
    else:
        pass

    # if previous submission file found and not empty
    # reformat it so the filetype.1, filename.1 and MD5_checksum.1
    # will be empty and reformatted into the newline
    if args.previous_submission:
        pre_sub_df = sra_template_to_dict(pre_sub_f)
        pre_sub_df_seq = pre_sub_df["Sequence_Data"]
        if not pre_sub_df_seq.empty:
            pre_sub_reformatted = reformat_previous_sra(p_sra_df=pre_sub_df_seq)
            logger.info(
                f"Found {pre_sub_reformatted.shape[0]} records of sequencing files in the previous submission"
            )
            # combine previous submission with current submisison
            sra_df = pd.concat([pre_sub_reformatted, sra_df], axis=0, ignore_index=True)
            # check and remove duplicates
            if sra_df.groupby(["library_ID", "filename"]).size().max() > 1:
                # if multiple lines are sharing same library_ID and same filename
                # they are considered as same record
                sra_df = check_and_remove_duplicates(sra_df=sra_df, logger=logger)
            else:
                pass
        else:
            pass
    else:
        pass

    # data frame manipulation, spread sra_df if multiple sequencing files are
    # sharing same library_ID. This function won't result multiple row sharing
    # same libary_ID
    sra_df = spread_sra_df(sra_df=sra_df, logger=logger)
    logger.info(
        "Sequencing files sharing the same library_ID will be reorganized into single row"
    )

    # rename any colnames, e.g. filetype.1 -> filetype
    sra_df = rename_colnames_output(sra_df=sra_df)
    logger.info(
        "Rename column names in SRA df, and same names for certain columns are expected."
    )

    # create file for final output
    acl_name = sra_df["phs_accession"].tolist()[0]
    sra_output_path = acl_name + "_SRA_submission.xlsx"
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
