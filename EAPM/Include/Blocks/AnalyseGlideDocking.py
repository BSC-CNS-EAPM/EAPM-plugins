from HorusAPI import PluginBlock, VariableTypes, PluginVariable, VariableList

import os

# Input variables
dockingFolderVariable = PluginVariable(
    id="docking_folder",
    name="Docking Folder",
    description="Folder containing the docking results",
    type=VariableTypes.FOLDER,
    defaultValue="docking",
)

# Variables
atom1SelectionVariable = PluginVariable(
    id="atom1_selection",
    name="Atom 1",
    description="Select the first atom to calculate the distance",
    type=VariableTypes.ATOM,
)

atom2SelectionVariable = PluginVariable(
    id="atom2_selection",
    name="Atom 2",
    description="Select the second atom to calculate the distance",
    type=VariableTypes.ATOM,
)

atomSelectionsGroupVariable = PluginVariable(
    id="atom_selections",
    name="Group",
    description="Name of the group selection to clusterize",
    type=VariableTypes.STRING,
)

selectionsListVariable = VariableList(
    id="selections_list",
    name="Selections",
    description="List of selections to analyse",
    prototypes=[
        atom1SelectionVariable,
        atom2SelectionVariable,
        atomSelectionsGroupVariable,
    ],
)


def analyseDockingAction(block: PluginBlock):
    folder_to_analyse = block.inputs.get("docking_folder")

    if folder_to_analyse is None or not os.path.isdir(folder_to_analyse):
        raise Exception("No valid docking folder selected")

    # Get the docking results


analyseGlideDocking = PluginBlock(
    name="Analyse Glide Docking",
    description="Analyse the docking results from Glide",
    inputs=[dockingFolderVariable],
    variables=[selectionsListVariable],
    outputs=[],
    action=analyseDockingAction,
)
