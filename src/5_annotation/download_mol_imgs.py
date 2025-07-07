import os
import argparse
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Draw


def parse_args():
    parser = argparse.ArgumentParser(description="Molecular Networking Annotation.")
    parser.add_argument(
        "--annotation_result_file",
        default="final_naive_annotation_results.csv",
        type=str,
        required=True,
        help="the annotation result file containing ID, SMILES, and other information",
    )
    parser.add_argument(
        "--structure_image_folder",
        default="mol_imgs/",
        type=str,
        required=True,
        help="the folder to save the structure images",
    )
    parser.add_argument(
        "--size",
        type=int,
        nargs=2,
        metavar=("W", "H"),
        default=(400, 400),
        help="Image size, e.g. --size 400 400",
    )
    parser.add_argument(
        "--format",
        choices=["png", "svg"],
        default="png",
        help="Image format (png or svg)",
    )
    return parser.parse_args()


def mol_from_smiles(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    Chem.rdDepictor.Compute2DCoords(mol)
    return mol


def save_image(mol, out_path, size=(400, 400), fmt="png"):
    if fmt == "png":
        img = Draw.MolToImage(mol, size=size)
        img.save(out_path, dpi=(300, 300))
    else:  # svg
        svg = Draw.MolToImage(mol, size=size, useSVG=True)
        out_path.write_text(svg, encoding="utf-8")


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.structure_image_folder, exist_ok=True)
    annotation_results = pd.read_csv(args.annotation_result_file)

    structure_paths = []
    for _, row in tqdm(annotation_results.iterrows(), total=len(annotation_results)):
        smiles = row.get("SMILES")
        if pd.isna(smiles):
            continue

        mol = mol_from_smiles(smiles)
        if mol is None:
            continue

        id_ = row.get("ID", "unknown_id")
        mol_img_filename = args.structure_image_folder + f"{id_}.{args.format}"
        save_image(mol, mol_img_filename, size=args.size, fmt=args.format)

        structure_paths.append(mol_img_filename)

    # 生成最终输出文件
    annotation_results["Structure"] = structure_paths
    # 重新排列
    columns_order = [
        "Round",
        "ID",
        "Seed Node",
        "Structure",
        "CID",
        "MonoIsotopic Weight",
        "SMILES",
        "Formula",
        "Score",
    ]
    annotation_results = annotation_results[columns_order]
    final_output_file = args.annotation_result_file.replace(
        ".csv", f"_with_structures.csv"
    )
    annotation_results.to_csv(final_output_file, index=False)
