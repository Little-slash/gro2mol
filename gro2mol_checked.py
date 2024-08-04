from collections import deque


def read_itp_atomtypes(filename):
    print("读取")
    atoms_type = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        read_atomtype = False
        for line in lines:
            if line.startswith('[ atoms ]'):
                read_atomtype = True
            elif line.startswith('['):
                read_atomtype = False
            elif read_atomtype:
                parts = line.split()
                if len(parts) > 6:
                    if parts[0] != ';':
                        index = parts[0]
                        atoms_type.append([index, parts[1], parts[6]])

        file.close()
    print("reads charges and atomtype successfully")
    return atoms_type


def read_itp_bonds(file_path):
    bonds = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        read_bonds = False
        for line in lines:
            if line.startswith('[ bonds ]'):
                read_bonds = True
            elif line.startswith('['):
                read_bonds = False
            elif read_bonds:
                parts = line.split()
                if len(parts) >= 2:
                    if parts[0] != ';':
                        atom1, atom2 = int(parts[0]), int(parts[1])
                        bonds.append([atom1, atom2])
    return bonds


def read_gro_xyz(file_path):
    print("正在读取gro文件中的原子坐标----------------------")
    atoms_coord = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        index = 0
        for line in lines:
            parts = line.split()
            if line.startswith('#'):
                break
            elif len(parts) >= 3:
                # parts = line.split()
                if parts[0] != ';':
                    mol = parts[0]
                    atoms_id = 'MOL' + '_' + mol + '_' + str(index)
                    dx, dy, dz = float(parts[len(parts) - 3])*10, float(parts[len(parts) - 2])*10, float(
                        parts[len(parts) - 1])*10
                    atoms_coord.append([mol, atoms_id, dx, dy, dz])
                    index += 1

        file.close()
    print("reads charges and mass successfully")
    return atoms_coord


# 定义一个函数，将队列中的元素按顺序写入文件
def write_queue_to_file(queue, file_path):
    with open(file_path, 'w') as file:
        while queue:
            element = queue.popleft()
            file.write(element + '\n')


def read_top(filename):
    molecule_sort = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        read_atomtype = False
        for line in lines:
            if line.startswith('[ molecules ]'):
                read_atomtype = True
            elif line.startswith('['):
                read_atomtype = False
            elif read_atomtype:
                parts = line.split()
                if parts[0] != ';':
                    molecule_num = int(parts[1])
                    molecule_sort.append([parts[0], molecule_num])
        file.close()
    print("reads molecules sorts successfully")
    return molecule_sort


# 创建一个空队列
out_force_queue = deque()
filename_top = './sp_12.top'
filename_gro = './sp_12.gro'
force_field_name = "GAFF_gro_cnt"
force_field_filename = "GAFF_gro_cnt.lt"
system_name = "Try"
gro_head = f"import \"{force_field_filename}\""
gro_last = "}"
out_gro_queue = deque()
out_gro_queue.append(gro_head)

atoms_list = read_gro_xyz(filename_gro)

molecules = read_top(filename_top)
# molecules = [['sp10',28]]
ide = 0
out_atom_queue = []
out_bond_queue = []
out_charges_queue = []
out_group_queue = []
data_atom_head = "write(\"Data Atoms\") {"
data_bond_head = "write(\"Data Bond List\") {"
data_charge_head = "write(\"In Charges\") {"
data_group_head = "write(\"In Settings\") {"
index_bond = 0

for mol_i in range(0, len(molecules)):
    molecule_itp_name = molecules[mol_i][0]
    out_atom_queue.append(f"{molecule_itp_name} inherits {force_field_name} {'{'}")
    out_charges_queue.append(data_charge_head)
    out_bond_queue.append(data_bond_head)
    out_atom_queue.append(data_atom_head)
    out_group_queue.append(data_group_head)
    filename_itp = molecule_itp_name + '.itp'

    atoms_single = read_itp_atomtypes(filename_itp)
    bond_single = read_itp_bonds(filename_itp)

    ide_last = ide + len(atoms_single) * molecules[mol_i][1]
    atoms_part = atoms_list[ide:ide_last]
    ide = ide_last

    for i_part in range(0, molecules[mol_i][1]):
        molecule_id = 'MOL'+atoms_part[0+i_part * len(atoms_single)][0]
        out_group_queue.append(f"group\t{molecule_itp_name}\tmolecule\t$mol:{molecule_id}")
        for i_atom in range(0, len(atoms_single)):
            # atom_id = atoms_part[i_atom + i_part * len(atoms_single)][1]

            atom_mol = 'MOL'+atoms_part[i_atom + i_part * len(atoms_single)][0]
            atom_type = atoms_single[i_atom][1]
            atom_id = atom_mol +str(i_part) + '-' + str(i_atom)
            atomx = atoms_part[i_atom + i_part * len(atoms_single)][2]
            atomy = atoms_part[i_atom + i_part * len(atoms_single)][3]
            atomz = atoms_part[i_atom + i_part * len(atoms_single)][4]
            out_atom_queue.append(f"\t$atom:{atom_id}\t$mol:{atom_mol}\t@atom:{atom_type}"
                                  f"\t0.0000\t{atomx}\t{atomy}\t{atomz}")
            out_charges_queue.append(f"\t\t set atom $atom:{atom_id} charge {atoms_single[i_atom][2]}")

        for i_bond in range(0, len(bond_single)):
            bond_id = 'bond' + '-' + str(index_bond)
            index_bond += 1
            atom_i1 = bond_single[i_bond][0] + i_part * len(atoms_single)
            atom_i2 = bond_single[i_bond][1] + i_part * len(atoms_single)
            atom_mol1 = atoms_part[atom_i1 - 1][0]
            atom_mol2 = atoms_part[atom_i2 - 1][0]
            atom_id1 = 'MOL' + atom_mol1 + str(i_part) + '-' + str(bond_single[i_bond][0] - 1)
            atom_id2 = 'MOL' + atom_mol2 + str(i_part) + '-' + str(bond_single[i_bond][1] - 1)
            out_bond_queue.append(f"\t$bond:{bond_id}\t$atom:{atom_id1}\t$atom:{atom_id2}")
    out_group_queue.append("}")
    out_atom_queue.append("}")
    out_bond_queue.append("}")
    out_charges_queue.append("}")
    for i in out_atom_queue:
        out_gro_queue.append(i)
    for i in out_bond_queue:
        out_gro_queue.append(i)
    for i in out_group_queue:
        out_gro_queue.append(i)
    for i in out_charges_queue:
        out_gro_queue.append(i)
    out_atom_queue.clear()
    out_bond_queue.clear()
    out_charges_queue.clear()
    out_group_queue.clear()
    out_gro_queue.append("}")

    file_path = f'{system_name}{mol_i}.lt'

    # 将队列中的元素写入文件
    write_queue_to_file(out_gro_queue, file_path)
    out_gro_queue.clear()
    print("队列中的元素已按顺序写入文件:", file_path)
