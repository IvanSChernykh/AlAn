import sys as ss
import os
import pandas as pd
import copy as cp
import numpy as np

if "-h" in ss.argv:
    print("AlAn, version 1.3.0")
    print("Description of parameters (all are optional):")
    print("-h - type this message and exit")
    print("-i - input directory (default - current)")
    print("-o - output directory (default - current)")
    print("-c - core precent threshold for ideal groups (default - 0.9)")
    print("-e - output files type (default - excel, otherwise - csv)")
else:
    try:
        indir_index = ss.argv.index("-i")
    except ValueError:
        indir = os.getcwd()
    else:
        indir = ss.argv[indir_index + 1]
    try:
        outdir_index = ss.argv.index("-o")
    except ValueError:
        outdir = os.getcwd()
    else:
        outdir = ss.argv[outdir_index + 1]
        outdir_split = outdir.split("/")
        for outdir_split_counter in range(len(outdir_split)):
            try:
                os.mkdir("/".join(outdir_split[:outdir_split_counter + 1]))
            except FileExistsError:
                pass
    try:
        core_percent_threshold_index = ss.argv.index("-c")
    except ValueError:
        core_percent_threshold = 0.9
    else:
        core_percent_threshold = ss.argv[core_percent_threshold_index + 1]
    if "-e" in ss.argv:
        out_excel = False
    else:
        out_excel = True
    table = pd.read_csv(f"{indir}/genes/partition-ungrouped.tsv", sep = "\t")
    table = table.iloc[:, :-1]
    table["gene_end"] = [False]*len(table)
    for gene, loc_table in table.groupby("gene"):
        loc_table_sort = loc_table.sort_values("gene_block_stop", ascending = False)
        table.loc[loc_table_sort.iloc[0].name, "gene_end"] = True
    try:
        os.mkdir(f"{outdir}/mut_split")
    except FileExistsError:
        pass
    with open(f"{indir}/mutations/mut.tsv") as infile:
        line = infile.readline()
        line = infile.readline()
        line_split = line.split("\t")
        while line != "":
            with open(f"{outdir}/mut_split/{line_split[0]}.tsv", "w") as outfile:
                outfile.write(line)
                line = infile.readline()
                line_split = line.split("\t")
                while (line != "") and (line_split[0] == "."):
                    outfile.write(line)
                    line = infile.readline()
                    line_split = line.split("\t")

    def group_table_f(group_table_in, group_id, outdir, mut_path):
        result = {"gene": [], "ori": [], "start": [], "stop": [], "len": [], "start_rf": [], "stop_rf": [], "is_start": [], "is_stop": []}
        for index_counter in range(len(group_table_in)):
            result["gene"].append(group_table_in.iloc[index_counter, 3])
            result["ori"].append(group_table_in.iloc[index_counter, 7])
            result["start"].append(group_table_in.iloc[index_counter, 5] + 1)
            result["stop"].append(group_table_in.iloc[index_counter, 6] + 1)
            result["len"].append(group_table_in.iloc[index_counter, 6] - group_table_in.iloc[index_counter, 5] + 1)
            result["is_start"].append(group_table_in.iloc[index_counter, 8] == 0)
            result["is_stop"].append(group_table_in.iloc[index_counter, 10])
            try:
                with open(f"{mut_path}/{group_table_in.iloc[0, 4]}.tsv") as mutfile:
                    succ = False
                    while not succ:
                        line = mutfile.readline()
                        line_split = line.split("\t")
                        if line == "":
                            succ = True
                        elif line_split[1] != ".":
                            seq_split = line_split[1].split("_")
                            if group_table_in.iloc[index_counter, 1] <= group_table_in.iloc[index_counter, 2]:
                                if (seq_split[0] == group_table_in.iloc[index_counter, 0]) and (group_table_in.iloc[index_counter, 1] >= int(seq_split[1])) and \
                                (group_table_in.iloc[index_counter, 2] <= int(seq_split[2])):
                                    succ = True
                            else:
                                if (seq_split[0] == group_table_in.iloc[index_counter, 0]) and (group_table_in.iloc[index_counter, 1] <= int(seq_split[1])) and \
                                (group_table_in.iloc[index_counter, 2] >= int(seq_split[2])):
                                    succ = True
                    split_sum = 0
                    if line != "":
                        coord = int(line_split[2])
                        if line_split[3].strip() == "-":
                            split_sum += 1
                        elif line_split[3].strip() not in ["A", "T", "G", "C", "N"]:
                            split_sum += int(line_split[3]) - coord + 1
                        line = mutfile.readline()
                        line_split = line.split("\t")
                        if line != "":
                            coord = int(line_split[2])
                    while (line != "") and (coord < group_table_in.iloc[index_counter, 5]) and (line_split[1] == "."):
                        if line_split[3].strip() == "-":
                            split_sum += 1
                        elif line_split[3].strip() not in ["A", "T", "G", "C", "N"]:
                            split_sum += int(line_split[3]) - coord + 1
                        line = mutfile.readline()
                        line_split = line.split("\t")
                        if line != "":
                            coord = int(line_split[2])
                    if group_table_in.iloc[index_counter, 7] == 1:
                        result["start_rf"].append((group_table_in.iloc[index_counter, 5] + 1 - split_sum)%3 + 1)
                    else:
                        result["start_rf"].append(-((group_table_in.iloc[index_counter, 5] + 1 - split_sum)%3 + 1))
                    while (line != "") and (coord < group_table_in.iloc[index_counter, 6]) and (line_split[1] == "."):
                        if line_split[3].strip() == "-":
                            split_sum += 1
                        elif line_split[3].strip() not in ["A", "T", "G", "C", "N"]:
                            split_sum += int(line_split[3]) - coord + 1
                        line = mutfile.readline()
                        line_split = line.split("\t")
                        if line != "":
                            coord = int(line_split[2])
                    if group_table_in.iloc[index_counter, 7] == 1:
                        result["stop_rf"].append((group_table_in.iloc[index_counter, 6] - 1 - split_sum)%3 + 1)
                    else:
                        result["stop_rf"].append(-((group_table_in.iloc[index_counter, 6] - 1 - split_sum)%3 + 1))
            except FileNotFoundError:
                if group_table.iloc[index_counter, 7] == 1:
                    result["start_rf"].append((group_table_in.iloc[index_counter, 5] + 1)%3 + 1)
                    result["stop_rf"].append((group_table_in.iloc[index_counter, 6] - 1)%3 + 1)
                else:
                    result["start_rf"].append(-((group_table_in.iloc[index_counter, 5] + 1)%3 + 1))
                    result["stop_rf"].append(-((group_table_in.iloc[index_counter, 6] - 1)%3 + 1))
        group_table_out = pd.DataFrame(result)
        if out_excel:
            group_table_out.to_excel(f"{outdir}/{group_id}.xlsx")
        else:
            group_table_out.to_csv(f"{outdir}/{group_id}.csv")
        if (group_table_out["ori"] == 1).all():
            result_ori = 1
        elif (group_table_out["ori"] == -1).all():
            result_ori = -1
        else:
            result_ori = 0
        return((group_table_out["is_start"].all() and group_table_out["is_stop"].all()), result_ori)

    result = {"block": [], "block_type": [], "seq_num": [], "block_len": [], "ori": [], "start": [], "stop": [], "len": [], "core_start": [], "core_stop": [], "core_len": [], "genes_num": [], "is_inside": [], "core_%": [], "group_type": []}
    seqs = table["sequence"].unique()
    for seq in seqs:
        result[seq] = []
    for block_ori, loc_table in table.groupby(["npg_block", "npg_block_ori"]):
        starts = sorted(loc_table["npg_block_min"])
        stops = sorted(loc_table["npg_block_max"])
        starts_stops = [[0, True]]*(len(starts) + len(stops))
        starts_counter = 0
        stops_counter = 0
        while starts_counter < len(starts):
            if starts[starts_counter] <= stops[stops_counter]:
                starts_stops[starts_counter + stops_counter] = [starts[starts_counter], True]
                starts_counter += 1
            else:
                starts_stops[starts_counter + stops_counter] = [stops[stops_counter], False]
                stops_counter += 1
        for stops_final_counter in range(stops_counter, len(stops)):
            starts_stops[starts_counter + stops_final_counter] = [stops[stops_final_counter], False]
        last_coord = -1
        for starts_stops_counter in range(len(starts_stops)):
            if starts_stops[starts_stops_counter][1]:
                last_coord = starts_stops[starts_stops_counter][0]
            elif last_coord != -1:
                core_start = last_coord
                core_stop = starts_stops[starts_stops_counter][0]
                group_table = loc_table[(loc_table["npg_block_min"] <= core_start) & (loc_table["npg_block_max"] >= core_stop)]
                group_start = min(group_table["npg_block_min"])
                group_stop = max(group_table["npg_block_max"])
                block_split = block_ori[0].split("x")
                result["block"].append(block_ori[0])
                result["block_type"].append(block_split[0][0])
                result["seq_num"].append(int(block_split[0][1:]))
                try:
                    result["block_len"].append(int(block_split[1]))
                except ValueError:
                    result["block_len"].append(np.nan)
                result["ori"].append(block_ori[1])
                result["start"].append(group_start + 1)
                result["stop"].append(group_stop + 1)
                result["len"].append(group_stop - group_start + 1)
                result["core_start"].append(core_start + 1)
                result["core_stop"].append(core_stop + 1)
                result["core_len"].append(core_stop - core_start + 1)
                result["genes_num"].append(len(group_table))
                seqs_counter = 0
                for seq in seqs:
                    if seq in group_table["sequence"].unique():
                        result[seq].append(True)
                        seqs_counter += 1
                    else:
                        result[seq].append(False)
                last_coord = -1
                try:
                    os.mkdir(f"{outdir}/groups")
                except FileExistsError:
                    pass
                inside = group_table_f(group_table, len(result["block"]) - 1, f"{outdir}/groups", f"{outdir}/mut_split")[0]
                result["is_inside"].append(inside)
                core_percent = round((core_stop - core_start + 1)/(group_stop - group_start + 1), 4)*100
                result["core_%"].append(core_percent)
                if inside:
                    if core_percent/100 >= core_percent_threshold:
                        if seqs_counter == result["genes_num"][-1]:
                            result["group_type"].append(0)
                        else:
                            result["group_type"].append(1)
                    else:
                        result["group_type"].append(2)
                else:
                    result["group_type"].append(3)
    groups = pd.DataFrame(result)
    if out_excel:
        groups.to_excel(f"{outdir}/groups.xlsx")
    else:
        groups.to_csv(f"{outdir}/groups.csv")
    result_no_ori = {"block": [], "block_type": [], "seq_num": [], "block_len": [], "ori": [], "start": [], "stop": [], "len": [], "core_start": [], "core_stop": [], "core_len": [], "genes_num": [], "is_inside": [], "core_%": [], "group_type": []}
    for seq in seqs:
        result_no_ori[seq] = []
    for block, loc_table in table.groupby("npg_block"):
        starts = sorted(loc_table["npg_block_min"])
        stops = sorted(loc_table["npg_block_max"])
        starts_stops = [[0, True]]*(len(starts) + len(stops))
        starts_counter = 0
        stops_counter = 0
        while starts_counter < len(starts):
            if starts[starts_counter] <= stops[stops_counter]:
                starts_stops[starts_counter + stops_counter] = [starts[starts_counter], True]
                starts_counter += 1
            else:
                starts_stops[starts_counter + stops_counter] = [stops[stops_counter], False]
                stops_counter += 1
        for stops_final_counter in range(stops_counter, len(stops)):
            starts_stops[starts_counter + stops_final_counter] = [stops[stops_final_counter], False]
        last_coord = -1
        for starts_stops_counter in range(len(starts_stops)):
            if starts_stops[starts_stops_counter][1]:
                last_coord = starts_stops[starts_stops_counter][0]
            elif last_coord != -1:
                core_start = last_coord
                core_stop = starts_stops[starts_stops_counter][0]
                group_table = loc_table[(loc_table["npg_block_min"] <= core_start) & (loc_table["npg_block_max"] >= core_stop)]
                group_start = min(group_table["npg_block_min"])
                group_stop = max(group_table["npg_block_max"])
                block_split = block.split("x")
                result_no_ori["block"].append(block)
                result_no_ori["block_type"].append(block_split[0][0])
                result_no_ori["seq_num"].append(int(block_split[0][1:]))
                try:
                    result_no_ori["block_len"].append(int(block_split[1]))
                except ValueError:
                    result_no_ori["block_len"].append(np.nan)
                result_no_ori["start"].append(group_start + 1)
                result_no_ori["stop"].append(group_stop + 1)
                result_no_ori["len"].append(group_stop - group_start + 1)
                result_no_ori["core_start"].append(core_start + 1)
                result_no_ori["core_stop"].append(core_stop + 1)
                result_no_ori["core_len"].append(core_stop - core_start + 1)
                result_no_ori["genes_num"].append(len(group_table))
                seqs_counter = 0
                for seq in seqs:
                    if seq in group_table["sequence"].unique():
                        result_no_ori[seq].append(True)
                        seqs_counter += 1
                    else:
                        result_no_ori[seq].append(False)
                last_coord = -1
                try:
                    os.mkdir(f"{outdir}/groups_no_ori")
                except FileExistsError:
                    pass
                inside, ori = group_table_f(group_table, len(result_no_ori["block"]) - 1, f"{outdir}/groups_no_ori", f"{outdir}/mut_split")
                result_no_ori["is_inside"].append(inside)
                result_no_ori["ori"].append(ori)
                core_percent = round((core_stop - core_start + 1)/(group_stop - group_start + 1), 4)*100
                result_no_ori["core_%"].append(core_percent)
                if inside:
                    if core_percent/100 >= core_percent_threshold:
                        if seqs_counter == result_no_ori["genes_num"][-1]:
                            result_no_ori["group_type"].append(0)
                        else:
                            result_no_ori["group_type"].append(1)
                    else:
                        result_no_ori["group_type"].append(2)
                else:
                    result_no_ori["group_type"].append(3)
    groups_no_ori = pd.DataFrame(result_no_ori)
    if out_excel:
        groups_no_ori.to_excel(f"{outdir}/groups_no_ori.xlsx")
    else:
        groups_no_ori.to_csv(f"{outdir}/groups_no_ori.csv")