ids = []
chains = []
pdb_dict = {}

with open("/home/larajuneb/Honours/PROJECT_(721)/Coding/SeqPredNN/SeqPredNN-main/Hons-Project/cullpdb_pc90.0_res0.0-2.2_len40-10000_R0.25_Xray_d2023_07_18_chains31896") as file:
    for line in file:
        row = line.rstrip()
        if not row.startswith("PDBchain"):
            identifier = row[:row.find(" ")]
            id = identifier[:4]
            chain = identifier[4:]
            ids.append(id)
            chains.append(chain)
            pdb_dict[id] = chain

unique_ids = set(ids)

with open("batch_PDB_id_list.txt", mode="w") as file:
    file.write(", ".join(unique_ids))