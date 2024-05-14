from rdkit import Chem
import numpy as np
import pandas as pds
import tkinter as tk
from tkinter import filedialog
from enum import Enum
import pubchempy as pcp

import sys
from pathlib import Path

from PySide6.QtCore import QObject, Slot
from PySide6.QtGui import QGuiApplication
from PySide6.QtQml import QQmlApplicationEngine, QmlElement
from PySide6.QtQuickControls2 import QQuickStyle

class PeakType(Enum):
    identical = 1


# 自定义的数据结构PeakUnit,记录了一个峰可能的结果的数据，包括化合物的IUPAC名称，CAS号和是该化合物的概率等重要信息
class PeakUnit:

    def __init__(self, name: str, ID: str, CAS: str, pp: str, peak_seq: int):
        self.name = name
        self.ID = ID
        self.CAS = CAS
        self.pp = pp
        self.peak_seq = peak_seq
        return

    pass


# 自定义的数据结构Peaks,用于储存一个主要含有停留时间、峰序号和占比的数据结构，其中还有一个存储PeakUnit结构的集合
# 此外，还有一个指向当前选择的PeakUnit的int指针select_num
# 有两个用于评判可信度的数据：可信度，一致度
class Peaks:
    reliability = 0
    consistency = 0
    select_num = 0
    peak_unit_collections = []

    def __init__(self, peak_seq: int, stay_time: str, portion: str, pl: list[PeakUnit]):
        self.peak_seq = peak_seq
        self.portion = portion
        self.stay_time = stay_time
        self.peak_unit_collections = pl
        self.consistency, self.reliability = PeakMachine.GetRCVal(self)
        return

    def Select_Peak(self, select_type: int):
        return

    pass


class TxtMachine:

    # 此方法用于输入文件并删除多余的行，会输出一个去掉无用行后的原始的txt文件，名为output_file.txt
    @staticmethod
    def delete_useless(fp: str):
        with open(fp, 'r') as input_file:
            input_file_content = input_file.readlines()
        with open('bin/output_file.txt', 'w') as output_file:
            line_count = sum(1 for line in input_file_content)
            for i in range(17, line_count - 1):
                if input_file_content[i].isspace():
                    continue
                else:
                    new_line = input_file_content[i]
                    output_file.write(new_line)

    # 此方法根据输入一个删除多余后的txt文件，返回一个切分后的2-D array数据结构
    @staticmethod
    def txt_to_numpy_array(fp: str):
        import numpy as np

        with open(fp, 'r') as input_file:
            input_file_content = input_file.readlines()
        line_count = sum(1 for line in input_file_content)
        txt_array = np.zeros((line_count, 4), dtype='<U50')

        cur_raw = 0
        for i in range(0, len(input_file_content)):
            print(input_file_content[i])
            output_line_list = TxtLineMachine.line_to_list(input_file_content[i])
            print(output_line_list)
            cur_column = 0
            for sole_str in output_line_list:
                txt_array[cur_raw, cur_column] = sole_str
                cur_column += 1
            cur_raw += 1

        #
        deleted_row = []
        for i in range(line_count - 1, 1, -1):
            if txt_array[i, 1] == '':
                txt_array[i - 1, 0] = txt_array[i - 1, 0] + txt_array[i, 0]
                txt_array[i, 0] = ''
                deleted_row.append(i)

        txt_array = np.delete(txt_array, deleted_row, 0)


        return txt_array

    # 此方法根据处理后的2D array文件，返回一个自定义的数据结构Peaks的集合
    @staticmethod
    def get_peak_collection(para_array: np.ndarray):
        import re
        peaks_collection = []
        length = para_array.shape[0]
        for i in range(0, length):
            isNumber = re.fullmatch('\d*', para_array[i, 0])
            peak_unit_collection = []
            if isNumber:
                # ！！！注意，这里设置了range的右终点为4，如果更改了出峰的个数的话需要更改这一位置的数字
                for j in range(1, 4):
                    para_array[i + j, 2] = TxtLineMachine.delete_zero(para_array[i + j, 2])
                    peak_unit = PeakUnit(para_array[i + j, 0], para_array[i + j, 1],
                                         para_array[i + j, 2], para_array[i + j, 3],
                                         para_array[i, 0])
                    print('Add Peak: ' + peak_unit.name+' :' + peak_unit.peak_seq)
                    peak_unit_collection.append(peak_unit)
                peaks = Peaks(para_array[i, 0], para_array[i, 1],
                              para_array[i, 2], peak_unit_collection)
                peaks_collection.append(peaks)

        return peaks_collection

    pass


class PeakMachine:

    @staticmethod
    def SelectThePeak(peaks_collection: list[Peaks]):

        return

    # 该方法用于计算单个峰内不同选择的最高可信度以及一致性
    @staticmethod
    def GetRCVal(pk: Peaks):
        max_reliability = pk.peak_unit_collections[0].pp
        name_with_mr = pk.peak_unit_collections[0].name
        consistency = 0
        for i in range(1, len(pk.peak_unit_collections)):
            if pk.peak_unit_collections[i].name != name_with_mr:
                consistency += float(pk.peak_unit_collections[i].pp)

        consistency = round(1 - consistency / 200, 2)
        return consistency, max_reliability

    # 此方法用于将txt文件转化得到的array写入excel文件
    @staticmethod
    def PeaksToExcel(fn: str, pl: list[Peaks]):
        import pandas as pds
        header = ['NO.', 'ST.', 'N%', 'IUPAC NAME', 'CAS', 'Reliability', 'Consistency']
        df = pds.DataFrame(index=range(len(pl)), columns=header)
        i = 0
        for peak in pl:
            df.iloc[i, 0] = int(peak.peak_seq)
            df.iloc[i, 1] = float(peak.stay_time)
            df.iloc[i, 2] = float(peak.portion)
            df.iloc[i, 3] = peak.peak_unit_collections[peak.select_num].name
            df.iloc[i, 4] = peak.peak_unit_collections[peak.select_num].CAS
            df.iloc[i, 5] = int(peak.reliability)
            df.iloc[i, 6] = float(peak.consistency)
            i += 1
        df.to_excel(fn, index=False)
        return

    pass


class TxtLineMachine:

    # 此方法用于去掉逗号后的空格，是为了进一步处理的子方法
    @staticmethod
    def delete_commas_blank(line: str):
        for i in range(1, len(line)):
            try:
                if line[i] == ' ' and line[i - 1] == ',':
                    line = line.replace(', ', ',')
            except:
                break
        return line

    # 此方法用于将单个txt文件的行切割为一个多个字符串对象组成的列表，以空格为切分符，是方法txt_to_numpy_array的依赖子方法
    @staticmethod
    def line_to_list(line: str):
        import re
        start_position = []
        end_position = []
        str_list = []
        for i in range(0, len(line)):
            if i == 0:
                if line[i] == ' ':
                    continue
                else:
                    start_position.append(i)
            elif i == len(line) - 1:
                if line[i] == ' ':
                    continue
                else:
                    end_position.append(i)
            else:
                if line[i] != ' ' and line[i - 1] == ' ':
                    start_position.append(i)
                if line[i] != ' ' and line[i + 1] == ' ':
                    end_position.append(i)

        for i in range(0, len(start_position)):
            if start_position[i] == end_position[i]:
                added_str = line[start_position[i]]
            else:
                added_str = line[start_position[i]:end_position[i] + 1]
            added_str.rstrip()
            str_list.append(added_str)

        while len(str_list) > 1:
            is_not_num_seq = not re.fullmatch('\d*', str_list[1]) or len(str_list[1]) == 1
            is_not_decimal = not re.fullmatch('\d*[.]\d*', str_list[1])
            if is_not_decimal and is_not_num_seq:
                str_list[0] = str_list[0] + ' ' + str_list[1]
                del str_list[1]
            else:
                break

        return str_list

    @staticmethod
    def delete_zero(cas: str):
        for i in cas:
            if i == '0':
                cas = cas[1:]
            else:
                break
        return cas

    pass


class MolecularCoping:
    @staticmethod
    def find_molecular():
        # 加载Excel文件

        df = pds.read_excel('T to E results/test excel.xlsx')
        df.insert(loc=7, column='Molecular Formula', value='Not Found', allow_duplicates=False)
        df.insert(loc=8, column='SMILES', value='Not Found', allow_duplicates=False)
        df.insert(loc=9, column='Carbon Number', value='Not Found', allow_duplicates=False)
        df.insert(loc=10, column='CateGory', value='Not Found', allow_duplicates=False)
        df.insert(loc=11, column='Sub CateGory', value='Not Found', allow_duplicates=False)
        filename = '705'
        last_formula = 'C1H1'
        last_smiles = 'C'
        name = 'Default'
        for i in range(0, len(df)):
            print('start searching MOLECULAR STRUCTURE...:' + df.iat[i, 3])
            try:
                name = df.iat[i, 3]
                compound = pcp.get_compounds(identifier=df.iat[i, 3], namespace='name')[0]
                print('add: ' + compound.molecular_formula)
                print('add: ' + compound.canonical_smiles)

                df.iat[i, 7] = compound.molecular_formula
                df.iat[i, 8] = compound.canonical_smiles
                df.iat[i, 9] = MolecularCoping.get_carbon_number(compound.molecular_formula)[0]
                df.iat[i, 10] = MolecularCoping.get_category(compound.molecular_formula,
                                                             compound.canonical_smiles,
                                                             name)[0]
                df.iat[i, 11] = MolecularCoping.get_category(compound.molecular_formula,
                                                             compound.canonical_smiles,
                                                             name)[1]

                print('add carbon number: ' + str(df.iat[i, 9]))
                print('add category: ' + str(df.iat[i, 10]))
                last_formula = compound.molecular_formula
                last_smiles = compound.canonical_smiles

                df.to_excel('searching results/' + filename + '.xlsx', sheet_name=filename, index=False)

            except:

                print('search err, return not found')
                df.iat[i, 7] = last_formula
                df.iat[i, 8] = last_smiles
                df.iat[i, 9] = MolecularCoping.get_carbon_number(last_formula)[0]
                df.iat[i, 10] = MolecularCoping.get_category(compound.molecular_formula,
                                                             compound.canonical_smiles,
                                                             name)[0]
                df.iat[i, 11] = MolecularCoping.get_category(compound.molecular_formula,
                                                             compound.canonical_smiles,
                                                             name)[1]

                print('add carbon number: ' + str(df.iat[i, 9]))
                print('add category: ' + str(df.iat[i, 10]))
                df.to_excel('searching results/' + filename + '.xlsx', sheet_name=filename, index=False)
        return

    @staticmethod
    def get_carbon_number(formula: str):
        import re
        carbon_num = 0
        hydro_num = 0
        C1H1 = re.fullmatch('C\dH\d', formula) or re.fullmatch('C\dH\d\D.*', formula)
        C1H2 = re.fullmatch('C\dH\d\d', formula) or re.fullmatch('C\dH\d\d\D.*', formula)
        C2H2 = re.fullmatch('C\d\dH\d\d', formula) or re.fullmatch('C\d\dH\d\d\D.*', formula)

        if C1H1:
            carbon_num = formula[1]
            hydro_num = formula[3]
        elif C1H2:
            carbon_num = formula[1]
            hydro_num = formula[3:5]
        elif C2H2:
            carbon_num = formula[1:3]
            hydro_num = formula[4:6]
        return carbon_num, hydro_num

    @staticmethod
    def get_category(formula: str, SMILES: str, IUPAC: str):
        import re
        category = ''
        carbon_num = MolecularCoping.get_carbon_number(formula)[0]
        hydro_num = MolecularCoping.get_carbon_number(formula)[1]
        # 以下是分子类型判断部分

        # 不是烃类的都被分到了杂原子Heterocyclic中，炔烃类的也会被分到杂原子中
        re1 = re.fullmatch('C\d*H\d*\D.*', formula)
        if re1:
            isHeterocyclic = True
        else:
            isHeterocyclic = False
        isAlkyne = re.fullmatch('.*#.*', SMILES)


        # 饱和度计算,如果是饱和烃类则设置属性isSaturated为真
        isSaturated = False
        saturation = (2 * int(carbon_num) + 2 == int(hydro_num))
        if saturation:
            isSaturated = True

        # 饱和的烃类判断为烷烃
        isParaffin = isSaturated and not isHeterocyclic

        # SMILES分子式中没有双键三键，且并非饱和的烃类判断为环烷烃
        isNaphthene = False
        hasDoubleBond = re.fullmatch('.*=.*|.*#.*', SMILES)
        isNaphthene = not isSaturated and not hasDoubleBond and not isHeterocyclic

        # 芳香烃判据:借用RDKit库实现
        mol: Chem.Mol
        mol = Chem.MolFromSmiles(SMILES)
        ring_info = mol.GetRingInfo()
        aromatic_info = mol.GetAromaticAtoms()
        atom_info = mol.GetAtoms()

        isAromatic = False
        if aromatic_info and not isHeterocyclic:
            isAromatic = True

        # 烯烃判据:没有芳香原子而且不饱和的分子

        if re.fullmatch('.*=.*', SMILES) and not isAromatic and not isHeterocyclic:
            isOlefin = True
        else:
            isOlefin = False

        if isHeterocyclic or isAlkyne:
            category = 'Heteroatomic'
        if isAromatic:
            category = 'Aromatic'
        if isParaffin:
            category = 'Paraffin'
        if isOlefin:
            category = 'Olefin'
        if isNaphthene:
            category = 'Naphthene'

        sub_category_list = []

        # 判断是不是单烯烃：
        if category == 'Olefin' and re.fullmatch('[^=]*=[^=]', SMILES):
            sub_category_list.append('mono-olefin')

        # 判断是不是α烯烃：
        if category == 'Olefin' and re.fullmatch('C=[^=]*|[^=]*=C]', SMILES):
            sub_category_list.append('α-Olefin')

        # 判断是不是环状的烯烃：
        if category == 'Olefin' and ring_info.NumRings() != 0:
            sub_category_list.append('r-Olefin')

        # 判断烷烃是不是直链的：
        if category == 'Paraffin' and re.fullmatch('.*[(].*[)].*', formula):
            sub_category_list.append('n-Paraffin')
        elif category == 'Paraffin' and not re.fullmatch('.*[(].*[)].*', formula):
            sub_category_list.append('i-Paraffin')

        # 判断芳烃是否属于BTX：
        if category == 'Aromatic' and IUPAC == 'Benzene' or IUPAC == 'Toluene' or re.fullmatch('.*xylene', IUPAC):
            sub_category_list.append('BTX')

        # 判断杂原子的种类：
        if category == 'Heteroatomic':
            for atom in formula:
                if atom != 'C' and atom != 'H' and not re.match('\d', atom):
                    print(atom)
                    print('Find the one:' + atom)
                    sub_category_list.append(atom)
        print(category)
        print(sub_category_list)

        return category, sub_category_list


class OilMaggotGUI:

    def activate_ui(self):


        return
    pass


class RingHydrocarbonMachine:

    @staticmethod
    def get_ring_list(SMILES: str):
        ring_list = []
        ring_struc = ''
        return ring_list

    @staticmethod
    def isAromatic(ring_struc: str):
        isAromatic = False

        return isAromatic

    pass

if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename()
    TxtMachine.delete_useless(file_path)
    file_path = 'bin/output_file.txt'
    a = TxtMachine.txt_to_numpy_array(file_path)
    pk = TxtMachine.get_peak_collection(a)
    PeakMachine.PeaksToExcel('T to E results/test excel.xlsx', pk)
    MolecularCoping.find_molecular()
