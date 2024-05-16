"""
Bioml Classification
    | Wrapper class for the bioml Classification module.
    | Train classification models.
"""

# TODO Add to the documentation

import os

from HorusAPI import PluginVariable, SlurmBlock, VariableGroup, VariableList, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
inputLabelFile = PluginVariable(
    name="Input Label",
    id="input_label",
    description="The path to the labels of the training set in a csv format or string if it is inside training features.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["csv"],
)
inputLabelString = PluginVariable(
    name="Input Label",
    id="input_label",
    description="The labels of the training set in a string format.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
trainingFeatures = PluginVariable(
    name="Training Features",
    id="training_features",
    description="The file to where the training features are saved in excel or csv format.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["csv", "xlsx"],
)
fileGroup = VariableGroup(
    id="fileType_input",
    name="Input File",
    description="The input is a file",
    variables=[inputLabelFile, trainingFeatures],
)
stringGroup = VariableGroup(
    id="stringType_input",
    name="Input String",
    description="The input is a string",
    variables=[inputLabelString, trainingFeatures],
)

# ==========================#
# Variable outputs
# ==========================#
outputClassification = PluginVariable(
    name="Classification output",
    id="out_zip",
    description="The zip file to the output for the classification models",
    type=VariableTypes.FILE,
)

##############################
#       Other variables      #
##############################
trainingOutput = PluginVariable(
    name="Training Output",
    id="training_output",
    description="The path where to save the models training results.",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
scalerVar = PluginVariable(
    name="Scaler",
    id="scaler",
    description="Choose one of the scaler available in scikit-learn, defaults to zscore.",
    type=VariableTypes.STRING_LIST,
    allowedValues=["robust", "zscore", "minmax"],
    defaultValue=None,
)
kfoldParameters = PluginVariable(
    name="Kfold Parameters",
    id="kfold_parameters",
    description="The parameters for the kfold in num_split:test_size format ('5:0.2').",
    type=VariableTypes.STRING,
    defaultValue=None,
)
outliersVar = PluginVariable(
    name="Outliers",
    id="outliers",
    description="A list of outliers if any, the name should be the same as in the excel file with the filtered features, you can also specify the path to a file in plain text format, each record should be in a new line",
    type=VariableTypes.STRING,
    defaultValue=None,
)
budgetTime = PluginVariable(
    name="Budget Time",
    id="budget_time",
    description="The time budget for the training in minutes, should be > 0 or None.",
    type=VariableTypes.FLOAT,
    defaultValue=None,
)
precisionWeight = PluginVariable(
    name="Precision Weight",
    id="precision_weight",
    description="Weights to specify how relevant is the precision for the ranking of the different features.",
    type=VariableTypes.FLOAT,
    defaultValue=None,
)
recallWeight = PluginVariable(
    name="Recall Weight",
    id="recall_weight",
    description="Weights to specify how relevant is the recall for the ranking of the different features.",
    type=VariableTypes.FLOAT,
    defaultValue=None,
)
reportWeight = PluginVariable(
    name="Report Weight",
    id="report_weight",
    description="Weights to specify how relevant is the f1, precision and recall for the ranking of the different features with respect to MCC which is a more general measures of the performance of a model.",
    type=VariableTypes.FLOAT,
    defaultValue=None,
)
differenceWeight = PluginVariable(
    name="Difference Weight",
    id="difference_weight",
    description="How important is to have similar training and test metrics.",
    type=VariableTypes.FLOAT,
    defaultValue=None,
)
bestModels = PluginVariable(
    name="Best Models",
    id="best_models",
    description="The number of best models to select, it affects the analysis and the saved hyperparameters.",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
seedVar = PluginVariable(
    name="Seed",
    id="seed",
    description="The seed for the random state.",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
dropVar = PluginVariable(
    name="Drop",
    id="drop",
    description="The models to drop.",
    type=VariableTypes.STRING_LIST,
    defaultValue=None,
    allowedValues=[
        "lr",
        "knn",
        "nb",
        "dt",
        "svm",
        "rbfsvm",
        "gpc",
        "mlp",
        "ridge",
        "rf",
        "qda",
        "ada",
        "gbc",
        "lda",
        "et",
        "xgboost",
        "lightgbm",
        "catboost",
        "dummy",
    ],
)
selectedVar = PluginVariable(
    name="Selected",
    id="selected",
    description="The models to select.",
    type=VariableTypes.STRING_LIST,
    defaultValue=None,
    allowedValues=[
        "lr",
        "knn",
        "nb",
        "dt",
        "svm",
        "rbfsvm",
        "gpc",
        "mlp",
        "ridge",
        "rf",
        "qda",
        "ada",
        "gbc",
        "lda",
        "et",
        "xgboost",
        "lightgbm",
        "catboost",
        "dummy",
    ],
)
tuneVar = PluginVariable(
    name="Tune",
    id="tune",
    description="If to tune the best models.",
    type=VariableTypes.BOOLEAN,
    defaultValue=None,
)
plotVar = PluginVariable(
    name="Plot",
    id="plot",
    description="The plots to save.",
    type=VariableTypes.STRING_LIST,
    defaultValue=None,
    allowedValues=["learning", "confusion_matrix", "class_report", "pr", "auc"],
)
optimizeVar = PluginVariable(
    name="Optimize",
    id="optimize",
    description="The metric to optimize for retuning the best models.",
    type=VariableTypes.STRING_LIST,
    defaultValue=None,
    allowedValues=["MCC", "Prec.", "Recall", "F1", "AUC", "Accuracy", "Average Precision Score"],
)
sheetName = PluginVariable(
    name="Sheet Name",
    id="sheet_name",
    description="The sheet name for the excel file if the training features is in excel format.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
numIter = PluginVariable(
    name="Number of Iterations",
    id="num_iter",
    description="The number of iterations for the hyperparameter search.",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
splitStrategy = PluginVariable(
    name="Split Strategy",
    id="split_strategy",
    description="The strategy to split the data.",
    type=VariableTypes.STRING_LIST,
    defaultValue=None,
    allowedValues=["mutations", "cluster", "stratifiedkfold", "kfold"],
)
clusterVar = PluginVariable(
    name="Cluster",
    id="cluster",
    description="The path to the cluster file generated by mmseqs2 or a custom group index file just like data/resultsDB_clu.tsv.",
    type=VariableTypes.FILE,
    defaultValue=None,
)
mutationsVar = PluginVariable(
    name="Mutations",
    id="mutations",
    description="The column name of the mutations in the training data.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
testNumMutations = PluginVariable(
    name="Test Number of Mutations",
    id="test_num_mutations",
    description="The threshold of number of mutations to be included in the test set.",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
greaterVar = PluginVariable(
    name="Greater",
    id="greater",
    description="Include in the test set, mutations that are greater of less than the threshold, default greater.",
    type=VariableTypes.BOOLEAN,
    defaultValue=None,
)
shuffleVar = PluginVariable(
    name="Shuffle",
    id="shuffle",
    description="If to shuffle the data before splitting.",
    type=VariableTypes.BOOLEAN,
    defaultValue=None,
)
crossValidation = PluginVariable(
    name="Cross Validation",
    id="cross_validation",
    description="If to use cross validation, default is True.",
    type=VariableTypes.BOOLEAN,
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

classificationBioMLBlock = SlurmBlock(
    name="Classification BioML",
    initialAction=runClassificationBioml,
    finalAction=finalAction,
    description="Train classification models.",
    inputGroups=[fileGroup, stringGroup],
    variables=BSC_JOB_VARIABLES
    + [
        selectedVar,
        dropVar,
        trainingOutput,
        scalerVar,
        kfoldParameters,
        outliersVar,
        budgetTime,
        precisionWeight,
        recallWeight,
        reportWeight,
        differenceWeight,
        bestModels,
        seedVar,
        tuneVar,
        plotVar,
        optimizeVar,
        sheetName,
        numIter,
        splitStrategy,
        clusterVar,
        mutationsVar,
        testNumMutations,
        greaterVar,
        shuffleVar,
        crossValidation,
    ],
    outputs=[outputClassification],
)
