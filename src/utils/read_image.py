from pathlib import Path

from PIL import Image

Image.MAX_IMAGE_PIXELS = None
WORKDIR = Path.joinpath(Path.home(), "workspace/svgbit_local_test")


def read_mba_image(gsm, sample):
    image_path = Path.joinpath(
        WORKDIR, f"data/2020_SciAdv_MouseBrainAtlus/HE/{gsm}_HE_{sample}.jpg")
    return Image.open(image_path)


def read_dlpfc_image(sample):
    image_path = Path.joinpath(WORKDIR,
                               f"data/spatialLIBD/{sample}_full_image.tif")
    return Image.open(image_path)
