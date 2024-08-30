"""
Module containing the rDock parameter file block for the EAPM plugin
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
receptorFile = PluginVariable(
    name="Receptor File",
    id="receptor_file",
    description="Prepared receptor file.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["mol2"],
)
referenceLigandFile = PluginVariable(
    name="Reference Ligand File",
    id="ref_ligand_file",
    description="File for the prepared ligand in the bound position without the protein",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["sd", "sdf"],
)

# ==========================#
# Variable outputs
# ==========================#

parameterFile = PluginVariable(
    name="Parameter File",
    id="parameter_file",
    description="The parameter file required for rDock.",
    type=VariableTypes.FILE,
    defaultValue="parameter_file.prm",
)

##############################
#       Other variables      #
##############################

receptorFlex = PluginVariable(
    name="Receptor flexibility",
    id="receptor_flex",
    description="Defines terminal OH and NH3+ groups within this distance of docking volume as flexible (angstroms).",
    type=VariableTypes.FLOAT,
    defaultValue=3.0,
)
radius = PluginVariable(
    name="Radius",
    id="radius",
    description="Radius of cavity mapping region (angstroms).",
    type=VariableTypes.FLOAT,
    defaultValue=6.0,
)
minVolume = PluginVariable(
    name="Minimum volume",
    id="min_volume",
    description="Minimum volume accepted to try the docking (angstroms^3).",
    type=VariableTypes.FLOAT,
    defaultValue=100,
)
gridStep = PluginVariable(
    name="Grid step",
    id="grid_step",
    description="Grid resolution for mapping (angstroms).",
    type=VariableTypes.FLOAT,
    defaultValue=0.5,
)


# Align action block
def paramFileRbdock(block: PluginBlock):

    import os

    if block.remote.name != "Local":
        raise Exception("This block is only available for local execution.")

    # Loading plugin variables
    # 1. receptor file
    receptor_file = block.inputs.get(receptorFile.id, None)
    if receptor_file is None:
        raise Exception("No receptor file provided.")
    if not os.path.exists(receptor_file):
        raise Exception("Parameter file does not exist.")

    # 2. reference ligand file
    reference_ligand_file = block.inputs.get(referenceLigandFile.id, None)

    if reference_ligand_file is None:
        raise Exception("No reference ligand file provided.")
    if not os.path.exists(reference_ligand_file):
        raise Exception("Reference ligand file does not exist.")

    # 3. directories
    path_receptor = os.path.dirname(receptor_file)
    path_reference_ligand = os.path.dirname(reference_ligand_file)

    if path_receptor != path_reference_ligand:
        print("WARNING! The receptor and reference ligand are not in the same directory.")
        print("The parameter file (.prm) will be saved in the receptor directory.")

    out = "parameter_file.prm"

    # 4. files paths
    receptor_file_path = os.path.join(path_receptor, receptor_file)
    reference_ligand_file_path = os.path.join(path_reference_ligand, reference_ligand_file)
    parameter_file = os.path.join(path_receptor, block.outputs.get(parameterFile.id, out))

    # 5. writing the parameter file
    with open(parameter_file, "w") as f:
        f.write("RBT_PARAMETER_FILE_V1.00\n")
        f.write("TILTE rdock\n")
        f.write("\n")
        f.write(f"RECEPTOR_FILE {receptor_file_path}\n")
        f.write(f"RECEPTOR_FLEX {block.variables.get('receptor_flex', 6.0)}\n")
        f.write("\n")
        f.write("##################################################################\n")
        f.write("### CAVITY DEFINITION: REFERENCE LIGAND METHOD\n")
        f.write("##################################################################\n")
        f.write("SECTION MAPPER\n")
        f.write("    SITE_MAPPER RbtLigandSiteMapper\n")
        f.write(f"    REF_MOL {reference_ligand_file_path}\n")
        f.write(f"    RADIUS {block.variables.get('radius', 6.0)}\n")
        f.write("    SMALL_SPHERE 1.0\n")
        f.write(f"    MIN_VOLUME {block.variables.get('min_volume', 100)}\n")
        f.write("    MAX_CAVITIES 1\n")
        f.write("    VOL_INCR 0.0\n")
        f.write(f"   GRIDSTEP {block.variables.get('grid_step', 0.5)}\n")
        f.write("END_SECTION\n")
        f.write("\n")
        f.write("#################################\n")
        f.write("#CAVITY RESTRAINT PENALTY\n")
        f.write("#################################\n")
        f.write("\n")
        f.write("SECTION CAVITY\n")
        f.write("    SCORING_FUNCTION RbtCavityGridSF\n")
        f.write("    WEIGHT 1.0\n")
        f.write("END_SECTION\n")
        f.write("\n")

    # Set the output
    block.setOutput(parameterFile.id, parameter_file)


rbParameterFileBlock = PluginBlock(
    name="rDockParameterFile",
    id="rbParameterFile",
    description="Generate the parameter file for rDock (for local).",
    action=paramFileRbdock,
    variables=[
        receptorFlex,
        radius,
        minVolume,
        gridStep,
    ],
    inputs=[receptorFile, referenceLigandFile],
    outputs=[parameterFile],
)
