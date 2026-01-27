#!/usr/bin/env python3
# -- coding: utf-8 --
"""
Created on December 11 5:08 PM 2025
Created in CLion
Created as mm_strip_reconstruction/process_run

@author: Dylan Neff, dylan
"""

import os
import re
from typing import Tuple, Optional


DECODE_EXECUTABLE = '/home/dylan/CLionProjects/mm_strip_reconstruction/cmake-build-debug/decoder/decode'
WAVEFORM_ANALYSIS_EXECUTABLE = '/home/dylan/CLionProjects/mm_strip_reconstruction/cmake-build-debug/waveform_analysis/analyze_waveforms'

decode = True
analyze = True

def main():
    # runs_dir = '/media/dylan/data/x17/nov_25_beam_test/dream_run/'
    runs_dir = '/media/dylan/data/x17/cosmic_bench/det_1/mx17_1-27-26/'
    runs = ['mx17_det1_1-27-26']
    # pedestal_dir = 'ped_thresh_1_12_25_18_30'
    # pedestal_dir = 'ped_thresh_30_11_25_14_40'
    # runs_dir = '/media/dylan/data/sps_beam_test_25/run_54/rotation_30_test_0/'
    # runs_dir = '/media/dylan/data/sps_beam_test_25/run_67/rotation_-60_banco_scan_0/'
    # runs = ['raw_daq_data']
    # pedestal_dir = 'dummy_peds'
    # pedestal_dir = ''
    pedestal_dir = 'same'

    # raw_dream_dir_name = 'raw_dream'
    raw_dream_dir_name = 'raw_daq_data'
    decoded_root_dir_name = 'decoded_root'
    hits_dir_name = 'hits_root'

    for run in runs:
        run_dir = os.path.join(runs_dir, run)
        run_dir_dirs = os.listdir(run_dir)
        for sub_run_name in run_dir_dirs:
            if not os.path.isdir(os.path.join(run_dir, sub_run_name)):
                continue

            sub_run_dir = os.path.join(run_dir, sub_run_name)
            raw_dream_dir = os.path.join(sub_run_dir, raw_dream_dir_name)
            if not os.path.exists(raw_dream_dir):
                print(f'No raw dream directory found at {raw_dream_dir}, skipping...')
                continue
            if analyze:
                make_dir_if_not_exists(f'{sub_run_dir}/{hits_dir_name}')
            print(f'Processing run directory: {raw_dream_dir}')
            fdf_files = [file for file in os.listdir(raw_dream_dir) if file.endswith('.fdf')]
            for file in fdf_files:
                file_path = os.path.join(raw_dream_dir, file)
                decoded_root_out_path = f'{sub_run_dir}/{decoded_root_dir_name}/{file.replace(".fdf", ".root")}'
                if decode:
                    print(f'\nDecoding {file_path}...')
                    decode_fdf(file_path, decoded_root_out_path)
                    print(f'Decoded to {decoded_root_out_path}')
                hits_out_path = f'{sub_run_dir}/{hits_dir_name}/{file.replace(".fdf", "_hits.root")}'
                if analyze:
                    print(f'\nAnalyzing waveforms in {decoded_root_out_path}...')
                    if pedestal_dir == 'same':
                        sub_run_ped_path = os.path.join(sub_run_dir, decoded_root_dir_name)
                    else:
                        sub_run_ped_path = os.path.join(runs_dir, sub_run_dir, pedestal_dir)
                    analyze_waveforms(decoded_root_out_path, sub_run_ped_path, hits_out_path)
                    print(f'Analyzed waveforms in {decoded_root_out_path} to {hits_out_path}\n')
    print('donzo')


def decode_fdf(file_path, out_root_path=None):
    """
    Decode .fdf file and extract relevant data.
    :param file_path: Path to the .fdf file.
    :param out_root_path: Optional path for the output .root file. If None, replaces .fdf with .root in the input path.
    """
    if out_root_path is None:
        out_root_path = file_path.replace('.fdf', '.root')
    else:
        make_dir_if_not_exists(os.path.dirname(out_root_path))
    command = f"{DECODE_EXECUTABLE} {file_path} {out_root_path}"
    os.system(command)


def analyze_waveforms(root_path, pedestal_dir, hits_out_path=None):
    """
    Analyze waveforms from the decoded .root file.
    :param root_path: Path to the decoded .root file.
    :param pedestal_dir: Directory containing pedestal files.
    :param hits_out_path: Optional path for the output hits .root file.
    If None, replaces .root with _hits.root in the input path.
    """
    file_num, feu_num = extract_file_numbers_tuple(os.path.basename(root_path))

    if pedestal_dir == 'same':
        pedestal_dir = os.path.dirname(root_path)

    ped_files = find_file_feu_file(pedestal_dir, file_num=None, feu_num=feu_num, file_extension='.root', flag='_pedthr_')
    if not ped_files:
        print(f'No pedestal file found for file number {file_num} and FEU number {feu_num}')
        ped_file_path = ''
    elif len(ped_files) > 1:
        print(f'Multiple pedestal files found for file number {file_num} and FEU number {feu_num}: {ped_files}')
        return
    else:
        ped_file_path = os.path.join(pedestal_dir, ped_files[0])
        print(f'Using pedestal file: {ped_file_path}')
    if hits_out_path is None:
        hits_out_path = root_path.replace('.root', '_hits.root')
    else:
        make_dir_if_not_exists(os.path.dirname(hits_out_path))
    print(f'\nhits_out_path: {hits_out_path}\n')
    command = f"{WAVEFORM_ANALYSIS_EXECUTABLE} {root_path} {hits_out_path} {ped_file_path}"
    print(f'Analyzing waveforms with command: {command}')
    os.system(command)


def find_file_feu_file(search_dir_path, file_num=None, feu_num=None, file_extension=None, flag=None) -> list:
    """
    Finds a list of files in the specified directory matching the given file number and FEU number.

    Args:
        search_dir_path (str): The directory to search in.
        file_num (int): The file number to match (xxx).
        feu_num (int): The FEU number to match (yy).
        file_extension (str): The file extension to filter by.
        Flag (str): A specific flag that should be present in the filename.
        Returns:
        List[str]: A list of matching filenames.
    """
    if file_num is None and feu_num is None and file_extension is None and flag is None:  # If no criteria provided, return all files
        return os.listdir(search_dir_path)

    matching_files = []
    for filename in os.listdir(search_dir_path):
        if file_extension is not None and not filename.endswith(file_extension):
            continue

        if flag is not None and flag not in filename:
            continue

        extracted = extract_file_numbers_tuple(filename)
        if extracted is not None:
            extracted_file_num, extracted_feu_num = extracted
            if file_num is not None and extracted_file_num != file_num:
                continue
            if feu_num is not None and extracted_feu_num != feu_num:
                continue
            matching_files.append(filename)
    return matching_files



def extract_file_numbers_tuple(filename: str) -> Optional[Tuple[int, int]]:
    """
    Extracts the file number (xxx) and FEU number (yy) from a filename
    following the pattern ..._xxx_yy.ext

    Args:
        filename (str): The name of the file.

    Returns:
        Optional[Tuple[int, int]]: A tuple (file_number, feu_number) as integers,
                                    or None if the pattern is not found.
    """
    # Pattern explanation:
    # r'.*': Matches any characters at the start (non-greedy)
    # '_': Matches the literal underscore
    # '(\d{3})': Capturing group 1: exactly three digits (xxx)
    # '_': Matches the literal underscore
    # '(\d{2})': Capturing group 2: exactly two digits (yy)
    # '\..*': Matches a dot (start of extension) followed by any characters
    pattern = r'.*_(\d{3})_(\d{2})\..*'

    match = re.match(pattern, filename)

    if match:
        # Group 1 is the file number (xxx), Group 2 is the FEU number (yy)
        file_num_str, feu_num_str = match.groups()
        # Convert the strings to integers and return
        return (int(file_num_str), int(feu_num_str))
    else:
        return None


def make_dir_if_not_exists(dir_path: str):
    """
    Creates the directory if it does not already exist.
    :param dir_path: Path to the directory.
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


if __name__ == '__main__':
    main()
