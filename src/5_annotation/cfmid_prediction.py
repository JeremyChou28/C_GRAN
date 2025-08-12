# # 对top10的smiles调用cfmid预测，并且利用matchms算出来modified cosine similarity写入文件
import os
import argparse
import subprocess
from multiprocessing import Pool, cpu_count
import pandas as pd
import numpy as np
from tqdm import tqdm
from matchms import Spectrum
from matchms.similarity import ModifiedCosine
from functools import partial


def parse_args():
    parser = argparse.ArgumentParser(description="CFM prediction.")
    # offline parameters
    parser.add_argument(
        "--num_containers",
        default=10,
        type=int,
        required=True,
        help="the number of Docker containers to run in parallel",
    )
    parser.add_argument(
        "--tolerance",
        default=0.1,
        type=float,
        required=True,
        help="the tolerance for the modified cosine similarity calculation",
    )
    parser.add_argument(
        "--energy_level",
        default=0,
        type=int,
        required=True,
        choices=[0, 1, 2],
        help="the energy level of the spectrum data predicted by CFM-ID, 0 represents low, 1 represents medium, and 2 represents high energy level",
    )
    parser.add_argument(
        "--ion_mode",
        default="positive",
        type=str,
        required=True,
        choices=["positive", "negative"],
        help="the ion mode of the spectrum data predicted by CFM-ID, positive or negative",
    )
    parser.add_argument(
        "--spectrum_file",
        default="./test_files/compounds_spectrum.mgf",
        type=str,
        required=True,
        help="the mgf file containing target node spectra",
    )
    parser.add_argument(
        "--modified_cosine_similarity_threshold",
        default=0.5,
        type=float,
        required=True,
        help="the threshold for modified cosine similarity to filter the results",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        required=True,
        default=10,
        help="the top k candidates for annotation",
    )
    parser.add_argument(
        "--round_num",
        type=int,
        default=0,
        help="the round of annotation, used to name the output file",
    )
    parser.add_argument(
        "--tmp_name",
        default="test",
        type=str,
        help="directory to save temporary files",
    )
    return parser.parse_args()


class CFMDockerController:
    def __init__(self, ion_mode, container_name="cfmid_runner"):
        self.container_name = container_name
        self.ion_mode = ion_mode

    def start_container(self):
        cmd = (
            f"docker run -dit --name {self.container_name} "
            f"-v $(pwd):/cfmid/public/ wishartlab/cfmid:latest sh"
        )
        subprocess.run(cmd, shell=True, check=True)
        print(f"✅ 容器 {self.container_name} 已启动")

    def predict(
        self,
        smiles,
        output_file,
        ppm=0.001,
        model_dir="/trained_models_cfmid4.0/[M+H]+",
        log_file="param_output.log",
        config_file="param_config.txt",
    ):
        model_dir = (
            model_dir
            if self.ion_mode == "positive"
            else model_dir.replace("[M+H]+", "[M-H]-")
        )
        cmd = (
            f"docker exec {self.container_name} "
            f'sh -c "cd /cfmid/public && '
            f"cfm-predict '{smiles}' {ppm} {model_dir}/{log_file} {model_dir}/{config_file} 0 {output_file}\""
        )
        # print(f"⏳ 正在预测: {smiles} → {output_file}")
        subprocess.run(cmd, shell=True, check=True)

    def stop_container(self):
        subprocess.run(f"docker stop {self.container_name}", shell=True, check=True)
        subprocess.run(f"docker rm {self.container_name}", shell=True, check=True)
        print(f"🛑 容器 {self.container_name} 已停止并删除")


def calculate_modified_cosine_similarity(
    prediction_ms, prediction_peaks, target_node_id, all_target_spectra, tolerance
):
    target_node_ms = all_target_spectra.get(target_node_id, None)["params"]["PEPMASS"]
    target_node_peaks = all_target_spectra.get(target_node_id, None)["peaks"]

    target_node_mz_list, target_node_intensity_list = [], []
    for mz, intensity in target_node_peaks:
        target_node_mz_list.append(mz)
        target_node_intensity_list.append(intensity)

    prediction_mz_list, prediction_intensity_list = [], []
    for mz, intensity in prediction_peaks:
        prediction_mz_list.append(mz)
        prediction_intensity_list.append(intensity)

    prediction_spectrum = Spectrum(
        mz=np.array(prediction_mz_list),
        intensities=np.array(prediction_intensity_list),
        metadata={"precursor_mz": prediction_ms},
    )
    target_spectrum = Spectrum(
        mz=np.array(target_node_mz_list),
        intensities=np.array(target_node_intensity_list),
        metadata={"precursor_mz": target_node_ms},
    )

    modified_cosine = ModifiedCosine(tolerance=tolerance)

    score = modified_cosine.pair(prediction_spectrum, target_spectrum)["score"]
    return score


def parse_groundtruth_spectrum(groundtruth_spectrum):
    spectra = {}
    current_spectrum = None

    with open(groundtruth_spectrum, "r") as f:
        for line in f:
            line = line.strip()
            if not line:  # 跳过空行
                continue

            if line.startswith("BEGIN IONS"):
                # 开始新的谱图
                current_spectrum = {"params": {}, "peaks": []}

            elif line.startswith("END IONS"):
                # 结束当前谱图
                if current_spectrum is not None:
                    spectra[current_spectrum["params"]["MATCH_ID"]] = current_spectrum
                current_spectrum = None

            else:
                # 如果在谱图块内，就处理具体的信息
                if current_spectrum is not None:
                    # 判断是否是 key=value 格式
                    if "=" in line:
                        key, value = line.split("=", 1)
                        current_spectrum["params"][key.strip()] = value.strip()
                    else:
                        # 否则认为是峰表数据
                        parts = line.split()
                        if len(parts) >= 2:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            current_spectrum["peaks"].append((mz, intensity))
                        # 如果有第三列可继续处理
        f.close()
    return spectra


def parse_prediction_spectrum(prediction_spectrum_file, energy_level):
    with open(prediction_spectrum_file, "r") as f:
        content = f.read()
        f.close()

    if not content.strip():
        os.remove(prediction_spectrum_file)
        print(f"已删除内容为空的文件: {prediction_spectrum_file}")
        return None, []

    lines = content.strip().splitlines()
    pmass = None
    peaks = []
    in_target_energy = False

    for line in lines:
        line = line.strip()
        if line.startswith("#PMass="):
            pmass = float(line.split("=")[-1])
        elif line.startswith("energy"):
            in_target_energy = int(line.split("energy")[-1]) == energy_level
        elif in_target_energy and line:
            try:
                mz, intensity = map(float, line.split())
                peaks.append((mz, intensity))
            except ValueError:
                pass  # Skip invalid lines
    return pmass, peaks


def process_smiles(args):
    (
        idx,
        energy_level,
        tolerance,
        modified_cosine_similarity_threshold,
        smile,
        output_file,
        container_idx,
        target_node_id,
        all_target_spectra,
    ) = args
    controller = controllers[container_idx]

    try:
        controller.predict(smile, output_file)
        prediction_ms, prediction_peaks = parse_prediction_spectrum(
            output_file, energy_level
        )

        if prediction_ms is None or not prediction_peaks:
            return idx, np.nan

        score = calculate_modified_cosine_similarity(
            prediction_ms,
            prediction_peaks,
            target_node_id,
            all_target_spectra,
            tolerance,
        )
        if score < modified_cosine_similarity_threshold:
            return idx, np.nan
        else:
            return idx, score
    except Exception as e:
        print(f"❌ 错误处理 {smile}: {e}")
        return idx, np.nan


def process_file(
    file,
    file_type,
    all_target_spectra,
    energy_level,
    tolerance,
    modified_cosine_similarity_threshold,
    topk,
    round_num,
    result_path,
):
    prefix = "not_unique" if file_type == "not_unique" else "unique"
    folder = (
        result_path
        + f"Seednode_and_Targetnode_Morgan_Similarity_score_split_{prefix}_Top{topk}_Round{round_num}/"
    )

    df = pd.read_csv(os.path.join(folder, file))
    smiles = df["SMILES"].values.tolist()
    target_node_id = file.split(".csv")[0]
    os.makedirs(result_path + f"cfmid_spectrum_results/{target_node_id}", exist_ok=True)

    # 准备任务
    tasks = []
    for i, smile in enumerate(smiles):
        output_file = result_path + f"cfmid_spectrum_results/{target_node_id}/{i}.txt"
        container_idx = i % NUM_CONTAINERS  # 简单轮询分配容器
        tasks.append(
            (
                i,
                energy_level,
                tolerance,
                modified_cosine_similarity_threshold,
                smile,
                output_file,
                container_idx,
                target_node_id,
                all_target_spectra,
            )
        )

    # 多进程并行
    with Pool(processes=NUM_CONTAINERS) as pool:
        results = pool.map(process_smiles, tasks)

    for idx, score in results:
        df.at[idx, "CFM-ID_score"] = score

    output_csv = result_path + f"cfmid_score_results/{target_node_id}.csv"
    df.to_csv(output_csv, index=False)


if __name__ == "__main__":
    args = parse_args()
    tmp_result_path = f"tmp/{args.tmp_name}/"

    # 读取not_unique和unique文件夹下的文件
    not_unique_files = os.listdir(
        tmp_result_path
        + f"Seednode_and_Targetnode_Morgan_Similarity_score_split_not_unique_Top{args.top_k}_Round{args.round_num}/"
    )
    unique_files = os.listdir(
        tmp_result_path
        + f"Seednode_and_Targetnode_Morgan_Similarity_score_split_unique_Top{args.top_k}_Round{args.round_num}/"
    )
    os.makedirs(tmp_result_path + "cfmid_score_results", exist_ok=True)

    # 读取目标节点的mgf文件并解析
    all_target_spectra = parse_groundtruth_spectrum(args.spectrum_file)

    NUM_CONTAINERS = min(
        args.num_containers, cpu_count()
    )  # 根据CPU核心数创建4个或更少容器
    controllers = []

    # 启动多个容器
    for i in range(NUM_CONTAINERS):
        cname = f"cfmid_runner_{i}"
        c = CFMDockerController(ion_mode=args.ion_mode, container_name=cname)
        c.start_container()
        controllers.append(c)

    for file in tqdm(not_unique_files, desc="Processing not_unique files"):
        if file.endswith(".csv"):
            process_file(
                file,
                "not_unique",
                all_target_spectra,
                args.energy_level,
                args.tolerance,
                args.modified_cosine_similarity_threshold,
                args.top_k,
                args.round_num,
                tmp_result_path,
            )

    for file in tqdm(unique_files, desc="Processing unique files"):
        if file.endswith(".csv"):
            process_file(
                file,
                "unique",
                all_target_spectra,
                args.energy_level,
                args.tolerance,
                args.modified_cosine_similarity_threshold,
                args.top_k,
                args.round_num,
                tmp_result_path,
            )

    for c in controllers:
        c.stop_container()
