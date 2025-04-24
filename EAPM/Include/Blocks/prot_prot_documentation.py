from HorusAPI import PluginBlock, Extensions
import shutil
import os


def open_docs(block: PluginBlock):

    pdf = "/home/perry/data/Programs/HorusServer/public_flows/prot_demo.pdf"

    pdf = shutil.copyfile(pdf, os.path.join(os.getcwd(), "workshop_tutorial.pdf"))

    Extensions().loadPDF(
        pdf,
        "Workshop Tutorial",
        store=False,
    )

    Extensions().loadPDF(
        pdf,
        "Workshop Tutorial",
    )


prot_docs = PluginBlock(
    id="prot_docs",
    name="Protein Protein Tutorial PDF",
    description="This block will open the tutorial PDF right into Horus",
    action=open_docs,
)
