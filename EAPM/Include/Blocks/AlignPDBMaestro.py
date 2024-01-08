import os
from os.path import abspath, basename, join, dirname
from glob import glob
import pandas as pd
from HorusAPI import PluginVariable, VariableTypes, PluginBlock


input_ = PluginVariable(
    name="input pdb folder",
    id="inp_path",
    description="The folder containing the alpha_fold models",
    type=VariableTypes.FOLDER
)


output = PluginVariable(
    name="Aligned pdb folder",
    id="out_path",
    description="Output folder to keep the models that passed the confidence threshold",
    type=VariableTypes.FOLDER,
    defaultValue="aligned_pdb"
)

maestro_path = PluginVariable(
    name="Schrodinger path",
    id="schro_path",
    description="Path to the Schrodinger installation path",
    type=VariableTypes.FOLDER,
)

output_name = PluginVariable(
    name="Output folder Name",
    id="out_name",
    description="Output folder name",
    type=VariableTypes.STRING,
    defaultValue="aligned_pdb"
)


ref_pdb = PluginVariable(
    name="Reference pdb",
    id="ref_pdb",
    description="Reference pdb file to align the models to",
    type=VariableTypes.FILE,
    allowedValues=["pdb"],
)


def rename(folder: str):
    files = glob(f"{folder}/*.pdb")
    for f in files:
        os.rename(f, f"{dirname(f)}/{basename(f).replace('rot-', '')}")


def align(ref_pdb: str, input_folder: str, path_to_maestro:str, 
          output_folder: str=None, rmsd_file: str="rmsd.txt"):
    command = f"{path_to_maestro}/utilities/structalign {abspath(ref_pdb)} {abspath(input_folder)}/*.pdb > {rmsd_file}"
    current = os.getcwd()
    if output_folder:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        os.chdir(output_folder)
        os.system(command)
        os.chdir(current)   
    else:
        os.system(command)


def extract_rmsd(rmsd_inpt: str, rmsd_out: str):
    with open(rmsd_inpt) as f:
        rmsd = f.readlines()
    data = {}
    for x in rmsd:
        if x.startswith('#  ct1:'):
            y = basename(x.split()[2].strip())
            data[y] = []
        if x.startswith('           RMSD:'):
            rmd = float(x.split()[1])
            data[y].append(rmd)
        if x.startswith("Alignment Score:"):
            score = float(x.split()[2])
            data[y].append(score)
    a = pd.DataFrame(data).T
    a.columns = ["Alignment Score", "RMSD"]
    a.to_csv(rmsd_out)
    


def align_rename(block: PluginBlock):
    """
    Aligns the pdb files in the input folder to the reference pdb file using structalign

    Parameters
    ----------
    block : PluginBlock
    """
    input_af = block.inputs.get("inp_path")
    schro_path = block.variables.get("schro_path")
    ref_pdb = block.inputs.get("ref_pdb")
    output_name = block.variables.get("out_name")
    output_dir = join(os.getcwd(), output_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print("aligning the models")
    align(ref_pdb, input_af, schro_path, output_dir)
    print("renaming the models")
    rename(output_dir)
    if os.path.exists(f"{output_dir}/{input_af}"):
        os.remove(f"{output_dir}/{input_af}")
    print("extracting the rmsd values for each model and saving in a csv file")
    extract_rmsd(f"{output_dir}/rmsd.txt", f"{output_dir}/rmsd.csv")



align_pdb_schro_block = PluginBlock(
    name="Structalign",
    description="Aligns the pdb files in the input folder to the reference pdb file using structalign",
    action=align_rename,
    variables=[maestro_path, output_name],
    inputs=[input_, ref_pdb],
    outputs=[output]
)