from collections import deque


def read_itp_atomtypes(filename):
    print("读取")
    set_type = []
    with open(filename, 'r') as f:
        flag = [1]
        while True:
            line = f.readline()
            if not line:
                break
            if line.find("atomtypes") != -1:
                flag[0] = 2
                continue
            elif line.find("[") != -1:
                break
            if line.find(';') != -1:
                continue
            elif flag[0] == 2:
                line.strip()
                line.replace('\n', '')
                sp = line.split()
                if len(sp) > 4:
                    sigma = float(sp[5]) * 10.000
                    epsilon = float(sp[6]) / 4.184
                    set_type.append([sp[0], sp[3], sp[2], epsilon, sigma])

        f.close()
    print("reads charges and mass successfully")
    return set_type


def puts_atomtypes(set_types):
    string_atom = []
    s_head = "write_once(\"In Charges\") {"
    s_last = "}\n"
    string_atom.append(s_head)
    for i in range(0, len(set_types)):
        string_s = f"\tset type @atom:{set_types[i][0]}     charge {set_types[i][1]}"
        string_atom.append(string_s)

    string_atom.append(s_last)
    s_head2 = "write_once(\"Data Masses\") {"
    string_atom.append(s_head2)

    for i in range(0, len(set_types)):
        string_s = f"\t@atom:{set_types[i][0]}    {set_types[i][2]}"
        string_atom.append(string_s)

    string_atom.append(s_last)

    s_head3 = "write_once(\"In Settings\") {"
    string_atom.append(s_head3)

    for i in range(0, len(set_types)):
        string_s = (f"\tpair_coeff @atom:{set_types[i][0]}  @atom:{set_types[i][0]} "
                    f"lj/charmm/coul/long {set_types[i][3]}  {set_types[i][4]}")
        string_atom.append(string_s)

    string_atom.append(s_last)

    return string_atom


def read_itp_atoms(filename):
    atoms = []
    with open(filename, 'r') as f:
        flag = [1]
        while True:
            line = f.readline()
            if not line:
                break
            if line.find("atoms") != -1:
                flag[0] = 2
                continue
            elif line.find("bond") != -1:
                break
            elif line.find("[") != -1:
                continue
            if line.find(';') != -1:
                continue
            elif flag[0] == 2:
                line.strip()
                line.replace('\n', '')
                sp = line.split()
                if len(sp) > 4:
                    atoms.append([sp[1], sp[6]])
        f.close()
    print("reads atoms successfully")
    return atoms


def read_itp_bonds(filename, atom_list):
    bond_coeff = []
    origin_bond = []
    with open(filename, 'r') as f:
        flag = [1]
        while True:
            line = f.readline()
            if not line:
                break
            if line.find("bonds") != -1:
                flag[0] = 2
                continue
            elif line.find("angle") != -1:
                break
            elif line.find("[") != -1:
                continue
            if flag[0] == 2:
                line.strip()
                line.replace('\n', '')
                sp = line.split()
                if len(sp) > 4:
                    if sp[0] == ';' or sp[1] == ';':
                        continue
                    r0 = float(sp[3]) * 10
                    k = float(sp[4]) / 418.4 / 2
                    b1 = int(sp[0]) - 1
                    b2 = int(sp[1]) - 1
                    spq = 1
                    bond_name_i = atom_list[b1][0] + '-' + atom_list[b2][0]
                    origin_bond.append([sp[0], sp[1]])
                    for i in range(0, len(bond_coeff)):
                        if bond_coeff[i][4] == bond_name_i:
                            spq = 2
                    if spq == 1:
                        bond_coeff.append([atom_list[b1][0], atom_list[b2][0],
                                           k, r0, bond_name_i, b1, b2])
        f.close()
    print("reads bonds successfully")
    return bond_coeff, origin_bond


def puts_band_coeff(bands):
    string_bond = []
    s_head = "write_once(\"In Settings\") {"
    s_last = "}\n"
    string_bond.append(s_head)
    for i in range(0, len(bands)):
        string_s = f"\tbond_coeff @bond:{bands[i][4]}  harmonic  {bands[i][2]}  {bands[i][3]} "
        string_bond.append(string_s)
    string_bond.append(s_last)

    s_head2 = "write_once(\"Data Bonds By Type\") {"
    string_bond.append(s_head2)
    for i in range(0, len(bands)):
        string_s = f"\t@bond:{bands[i][4]} @atom:{bands[i][0]} @atom:{bands[i][1]}"
        string_bond.append(string_s)
    string_bond.append(s_last)

    return string_bond


def read_itp_angles(filename, atom_list):
    angle_coeff = []
    with open(filename, 'r') as f:
        flag = [1]
        while True:
            line = f.readline()
            if not line:
                break
            if line.find("angles") != -1:
                flag[0] = 2
                continue
            elif line.startswith("[ dihedrals ]"):
                break
            elif line.find("[") != -1:
                continue
            if flag[0] == 2:
                line.strip()
                line.replace('\n', '')
                sp = line.split()
                if len(sp) > 6:
                    if sp[0] == ';' or sp[1] == ';':
                        continue
                    k = float(sp[5]) / 4.184 / 2
                    b1 = int(sp[0]) - 1
                    b2 = int(sp[1]) - 1
                    b3 = int(sp[2]) - 1
                    angle_name_i = atom_list[b1][0] + '-' + atom_list[b2][0] + '-' + atom_list[b3][0]
                    spq = 1
                    for i in range(0, len(angle_coeff)):
                        if angle_coeff[i][5] == angle_name_i:
                            spq = 2
                    if spq == 1:
                        angle_coeff.append([atom_list[b1][0], atom_list[b2][0], atom_list[b3][0],
                                            k, sp[4], angle_name_i])
        f.close()
    print("reads angles successfully")
    return angle_coeff


def puts_angle_coeff(angles):
    string_angle = []
    s_head = "write_once(\"In Settings\") {"
    s_last = "}\n"
    string_angle.append(s_head)
    for i in range(0, len(angles)):
        string_s = f"\tangle_coeff @angle:{angles[i][5]}  harmonic  {angles[i][3]}  {angles[i][4]} "
        string_angle.append(string_s)
    string_angle.append(s_last)

    s_head2 = "write_once(\"Data Angles By Type\") {"
    string_angle.append(s_head2)
    for i in range(0, len(angles)):
        string_s = f"\t@angle:{angles[i][5]} @atom:{angles[i][0]} @atom:{angles[i][1]}  @atom:{angles[i][2]}"
        string_angle.append(string_s)
    string_angle.append(s_last)

    return string_angle


def read_itp_diherals(filename, atom_list):
    diherals_coeff = []
    with (open(filename, 'r') as f):
        flag = [1]
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("[ dihedrals ]"):
                flag[0] = 2
                continue
            elif line.startswith("[ pairs ]"):
                break
            elif line.find("[") != -1:
                continue
            if flag[0] == 2:
                line.strip()
                line.replace('\n', '')
                sp = line.split()
                if len(sp) > 6:
                    if sp[0] == ';' or sp[1] == ';':
                        continue
                    k = float(sp[6]) / 4.184
                    b1 = int(sp[0]) - 1
                    b2 = int(sp[1]) - 1
                    b3 = int(sp[2]) - 1
                    b4 = int(sp[3]) - 1
                    diheral_name_i = atom_list[b1][0] + '-' + atom_list[b2][0] + '-' + atom_list[b3][0] + '-' + \
                                     atom_list[b4][0]
                    spq = 1
                    for i in range(0, len(diherals_coeff)):
                        if diherals_coeff[i][7] == diheral_name_i:
                            spq = 2
                    if spq == 1:
                        sp_deg = int(float(sp[5]))
                        diherals_coeff.append([atom_list[b1][0], atom_list[b2][0], atom_list[b3][0], atom_list[b4][0],
                                               k, sp_deg, sp[7], diheral_name_i])
        f.close()
    print("reads dihedral successfully")
    return diherals_coeff


def puts_dihedral_coeff(dihedrals):
    string_dihedral = []
    s_head = "write_once(\"In Settings\") {"
    s_last = "}\n"
    string_dihedral.append(s_head)
    for i in range(0, len(dihedrals)):
        string_s = (f"\tdihedral_coeff @dihedral:{dihedrals[i][7]}  charmm  {dihedrals[i][4]}  {dihedrals[i][6]}"
                    f" {dihedrals[i][5]} 0.5")
        string_dihedral.append(string_s)
    string_dihedral.append(s_last)

    s_head2 = "write_once(\"Data Dihedrals By Type\") {"
    string_dihedral.append(s_head2)
    for i in range(0, len(dihedrals)):
        string_s = (f"\t@dihedral:{dihedrals[i][7]} @atom:{dihedrals[i][0]} @atom:{dihedrals[i][1]}  "
                    f"@atom:{dihedrals[i][2]} @atom:{dihedrals[i][3]}")
        string_dihedral.append(string_s)
    string_dihedral.append(s_last)

    return string_dihedral


def read_itp_impropers(filename, atom_list):
    improper_coeff = []
    with (open(filename, 'r') as f):
        flag = [1]
        while True:
            line = f.readline()
            if not line:
                break
            if line.find("improper") != -1:
                flag[0] = 2
                continue
            elif line.find("[") != -1:
                continue
            if flag[0] == 2:
                line.strip()
                line.replace('\n', '')
                sp = line.split()
                if len(sp) > 6:
                    if sp[0] == ';' or sp[1] == ';':
                        continue
                    k = float(sp[6]) / 4.184
                    b1 = int(sp[0]) - 1
                    b2 = int(sp[1]) - 1
                    b3 = int(sp[3]) - 1
                    b4 = int(sp[4]) - 1
                    improper_name_i = atom_list[b1][0] + '-' + atom_list[b2][0] + '-' + atom_list[b3][0] + '-' + \
                                      atom_list[b4][0]
                    spq = 1
                    for i in range(0, len(improper_coeff)):
                        if improper_coeff[i][6] == improper_name_i:
                            spq = 2
                    if spq == 1:
                        improper_coeff.append([atom_list[b1][0], atom_list[b2][0], atom_list[b3][0], atom_list[b4][0],
                                               k, sp[7], improper_name_i])
        f.close()
    print("reads improper successfully")
    return improper_coeff


def puts_improper_coeff(impropers):
    string_improper = []
    s_head = "write_once(\"In Settings\") {"
    s_last = "}\n"
    string_improper.append(s_head)
    for i in range(0, len(impropers)):
        string_s = (f"\timproper_coeff @improper:{impropers[i][6]}  cvff  {impropers[i][4]}  -1\t"
                    f"{impropers[i][5]}")
        string_improper.append(string_s)
    string_improper.append(s_last)

    s_head2 = "write_once(\"Data Impropers By Type\") {"
    string_improper.append(s_head2)
    for i in range(0, len(impropers)):
        string_s = (f"\t@improper:{impropers[i][6]} @atom:{impropers[i][0]} @atom:{impropers[i][1]}  "
                    f"@atom:{impropers[i][2]} @atom:{impropers[i][3]}")
        string_improper.append(string_s)
    string_improper.append(s_last)

    return string_improper


# 定义一个函数，将队列中的元素按顺序写入文件
def write_queue_to_file(queue, file_path):
    with open(file_path, 'w') as file:
        while queue:
            element = queue.popleft()
            file.write(element + '\n')


out_force_queue = deque()
filename_itp = './sp10.itp'

force_field_name = "GAFF_gro_sp"
all_head = "GAFF_gro_sp {\n"
all_last = "}"
out_force_queue.append(all_head)
atomtypes = read_itp_atomtypes(filename_itp)
strings_atom = puts_atomtypes(atomtypes)
atoms = read_itp_atoms(filename_itp)

# ---------------------------------------------------#
# 读取 bond_coeff,angle_coeff,diheral_coeff,improper_coeff,输出信息
# ---------------------------------------------------#

bonds, origin_bonds = read_itp_bonds(filename_itp, atoms)
strings_bonds = puts_band_coeff(bonds)
angles = read_itp_angles(filename_itp, atoms)
string_angles = puts_angle_coeff(angles)

dihedral = read_itp_diherals(filename_itp, atoms)
string_dihedral = puts_dihedral_coeff(dihedral)
improper = read_itp_impropers(filename_itp, atoms)
string_impropers = puts_improper_coeff(improper)

for i in strings_atom:
    out_force_queue.append(i)
for i in strings_bonds:
    out_force_queue.append(i)
for i in string_angles:
    out_force_queue.append(i)
for i in string_dihedral:
    out_force_queue.append(i)
for i in string_impropers:
    out_force_queue.append(i)

string_init = ("write_once(\"In Init\") {\n# Default styles and settings for AMBER based force-fields:\n"
               "\tunits           real\n"
               "\tatom_style      full\n"
               "\tbond_style      hybrid harmonic\n"
               "\tangle_style     hybrid harmonic\n"
               "\tdihedral_style  hybrid charmm\n"
               "\timproper_style  hybrid cvff\n"
               "\tpair_style      hybrid lj/charmm/coul/long 9.0 10.0 10.0\n"
               "\tkspace_style    pppm 0.0001\n"
               "\tpair_modify     mix arithmetic\n"
               "\tspecial_bonds   charmm\n"
               "}\n")
out_force_queue.append(string_init)
out_force_queue.append(all_last)
print('pause')

# 向队列中添加元素

# 指定文件路径
file_path = f'{force_field_name}.lt'

# 将队列中的元素写入文件
write_queue_to_file(out_force_queue, file_path)

print("队列中的元素已按顺序写入文件:", file_path)