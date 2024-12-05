"""
Module containing the rDock parameter file block for the rDock plugin
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes, VariableGroup

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

boxCenter = PluginVariable(
    name="Cavity Sphere",
    id="cavity_box_center",
    description="Coordinates of the center of the cavity box.",
    type=VariableTypes.SPHERE,
    defaultValue=None,
)

radius = PluginVariable(
    name="Radius",
    id="radius",
    description="Radius of cavity mapping region (angstroms).",
    type=VariableTypes.FLOAT,
    defaultValue=10.0,
)

groupCavitySphere = VariableGroup(
    name="RbtSphereSiteMapper",
    id="RbtSphereSiteMapper",
    description="Cavity definition using a sphere.",
    variables=[receptorFile, boxCenter],
)

groupCavityLigand = VariableGroup(
    name="RbtLigandSiteMapper",
    id="RbtLigandSiteMapper",
    description="Cavity definition using a reference ligand.",
    variables=[receptorFile, referenceLigandFile, radius],
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
    allowedValues=["prm"],
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

    # 2. paths
    path_receptor = os.path.dirname(receptor_file)
    out = "parameter_file.prm"

    # 3. files paths
    receptor_file_path = os.path.join(path_receptor, receptor_file)
    parameter_file = block.outputs.get(parameterFile.id, out)

    # 4. writing the parameter file
    with open(parameter_file, "w") as f:
        f.write("RBT_PARAMETER_FILE_V1.00\n")
        f.write("TILTE rdock\n")
        f.write("\n")
        f.write(f"RECEPTOR_FILE {receptor_file_path}\n")
        f.write(f"RECEPTOR_FLEX {block.variables.get('receptor_flex', 6.0)}\n")
        f.write("\n")
        f.write("##################################################################\n")
        f.write("### CAVITY DEFINITION\n")
        f.write("##################################################################\n")
        f.write("SECTION MAPPER\n")
        if block.selectedInputGroup == groupCavitySphere.id:

            # Checking inputs
            sphere_position = block.inputs.get(boxCenter.id, None)
            print(sphere_position)

            if sphere_position is None:
                raise Exception("No cavity sphere provided.")

            # Getting center and radius of sphere
            center = sphere_position["center"]
            formatted_center = ",".join(center.values())
            radius_sphere = sphere_position["radius"]

            f.write("    SITE_MAPPER RbtSphereSiteMapper\n")
            f.write(f"    CENTER ({formatted_center})\n")
            f.write(f"    RADIUS {radius_sphere}\n")

        else:
            # Checking inputs
            reference_ligand_file = block.inputs.get(referenceLigandFile.id, None)
            if reference_ligand_file is None:
                raise Exception("No reference ligand file provided.")
            if not os.path.exists(reference_ligand_file):
                raise Exception("Reference ligand file does not exist.")

            # Declaring paths
            path_reference_ligand = os.path.dirname(reference_ligand_file)
            if path_receptor != path_reference_ligand:
                print("WARNING! The receptor and reference ligand are not in the same directory.")
                print("The parameter file (.prm) will be saved in the receptor directory.")

            reference_ligand_file_path = os.path.join(
                path_reference_ligand, reference_ligand_file
            )

            f.write("    SITE_MAPPER RbtLigandSiteMapper\n")
            f.write(f"    REF_MOL {reference_ligand_file_path}\n")
            f.write(f"    RADIUS {block.inputs.get('radius', 6.0)}\n")

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
        minVolume,
        gridStep,
    ],
    inputGroups=[groupCavityLigand, groupCavitySphere],
    outputs=[parameterFile],
)
