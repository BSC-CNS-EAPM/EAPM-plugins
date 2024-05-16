"""
A module that performs regression analysis on a dataset.
"""

import os

from HorusAPI import PluginVariable, SlurmBlock, VariableGroup, VariableList, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
excelFile = PluginVariable(
    name="Excel file",
    id="excel_file",
    description="The file to where the selected features are saved in excel format.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["xlsx"],
)


# ==========================#
# Variable outputs
# ==========================#
outputOutliers = PluginVariable(
    name="Outliers output",
    id="out_zip",
    description="The path to the output for the outliers.",
    type=VariableTypes.FILE,
)

##############################
#       Other variables      #
##############################
numThreads = PluginVariable(
    name="Number of threads",
    id="num_threads",
    description="The number of threads to use.",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
scalerVar = PluginVariable(
    name="Scaler",
    id="scaler",
    description="The scaler to use.",
    type=VariableTypes.STRING,
    defaultValue="StandardScaler",
    allowedValues=["StandardScaler", "MinMaxScaler", "RobustScaler"],
)
contaminationVar = PluginVariable(
    name="Contamination",
    id="contamination",
    description="The contamination value.",
    type=VariableTypes.FLOAT,
    defaultValue=None,
)
numFeatures = PluginVariable(
    name="Number of features",
    id="num_features",
    description="The number of features to use.",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)


def runClassificationBioml(block: SlurmBlock):

    input_excel = block.inputs.get("input_excel", None)
    if input_excel is None:
        raise Exception("No input excel provided")
    if not os.path.exists(input_excel):
        raise Exception(f"The input excel file does not exist: {input_excel}")

    input_hyperparameters = block.inputs.get("input_hyperparameters", None)
    if input_hyperparameters is None:
        raise Exception("No input hyperparameters provided")
    if not os.path.exists(input_hyperparameters):
        raise Exception(f"The input hyperparameters file does not exist: {input_hyperparameters}")

    input_sheets = block.inputs.get("input_sheets", None)
    if input_sheets is None:
        raise Exception("No input sheets provided")

    input_label = block.inputs.get("input_label", None)
    if input_label is None:
        raise Exception("No input label provided")
    if not os.path.exists(input_label):
        raise Exception(f"The input label file does not exist: {input_label}")

    command = "python -m BioML.Ensemble "
    command += f"--excel {input_excel} "
    command += f"--hyperparameter_path {input_hyperparameters} "
    command += f"--sheets {input_sheets} "
    command += f"--label {input_label} "

    jobs = [command]

    folderName = block.variables.get("folder_name", "ensembleBioml")
    block.extraData["folder_name"] = folderName
    removeExisting = block.variables.get("remove_existing_results", False)

    # If folder already exists, raise exception
    if removeExisting and os.path.exists(folderName):
        os.system("rm -rf " + folderName)

    if not removeExisting and os.path.exists(folderName):
        raise Exception(
            "The folder {} already exists. Please, choose another name or remove it.".format(
                folderName
            )
        )

    # Create an copy the inputs
    os.makedirs(folderName, exist_ok=True)
    os.system(f"cp {input_excel} {folderName}")
    os.system(f"cp {input_hyperparameters} {folderName}")
    os.system(f"cp {input_sheets} {folderName}")
    os.system(f"cp {input_label} {folderName}")

    from utils import launchCalculationAction

    launchCalculationAction(
        block,
        jobs,
        program="bioml",
        uploadFolders=[
            folderName,
        ],
    )


def finalAction(block: SlurmBlock):
    pass


from utils import BSC_JOB_VARIABLES

outliersBioMLBlock = SlurmBlock(
    name="Regression BioMl",
    initialAction=runClassificationBioml,
    finalAction=finalAction,
    description="Train Regression models.",
    inputs=[excelFile],
    variables=BSC_JOB_VARIABLES
    + [
        numThreads,
        scalerVar,
        contaminationVar,
        numFeatures,
    ],
    outputs=[outputOutliers],
)
