# OilMaggot 1.0 - Chemical Analysis Toolkit 

[English](#english) | [中文](#chinese)

---

<a name="english"></a>
## Introduction

**OilMaggot** is a chemical data analysis tool based on Python, using the public database PubChem and open-source API RDKit. It is designed for researchers to process chromatography data, identify compounds via PubChem, and visualize molecular distributions.

### Key Features
* **ANALYZE**: Automated pipeline for raw data analysis from Angient GC-MS. You can just input the raw data and get the final result! (also a GUI to make this process more clear)
* 
<img width="600" height="500" alt="image" src="https://github.com/user-attachments/assets/e0fbab9b-c160-4b87-81b4-e1be07338619" />

Input data example:

<img width="500" height="400" alt="image" src="https://github.com/user-attachments/assets/96784d7e-1b6e-45bd-a103-4812bdd67d86" />

......

Ouput result be like: 

<img width="900" height="150" alt="image" src="https://github.com/user-attachments/assets/aaec6d69-0bad-47ed-8a0d-65de4155b41e" />

<img width="900" height="200" alt="image" src="https://github.com/user-attachments/assets/15ba5fa8-2c93-412f-94d5-ff0d2c82adbe" />



* **Steps**:
  1. Click the **ANALYZE** button on the left panel.
  2. Select your raw **.txt** data file in the file dialog.
* **Output**: Results are saved in the `saving file` directory.
* **File Descriptions**:
  * **`*(pivot).xlsx` (Pivot Table)**: The primary summary file. It organizes data by **Carbon Number** (rows) and **Molecular Category** (columns).
  * **`*.xlsx` (Detailed Analysis)**:
    * **Category Column**: The family classification of the molecule.
    * **HIC, HC, LIC, LC Columns**: Represent calculation results from different resolutions and heteroatom handling logic. **These can generally be ignored for standard analysis.**





### 📥 Download
You can download the latest pre-compiled executable for Windows here:
> [**Download OilMaggot 1.0 (.zip)**](https://github.com/ranwuge/OilMaggot1.0/releases/latest)

---

<a name="chinese"></a>
## 中文简介

**OilMaggot 1.0** 是一款基于 Python 和 PySide6 开发的化学数据分析工具。用于处理色谱数据、通过 PubChem 数据库鉴定化合物并实现分子结构的统计可视化。

### 🚀 核心功能
* **数据分析 (ANALYZE)**：全自动流程，用于安捷特GC-MS色谱的数据结果的自动分析。


### 📥 软件下载
点击下方链接下载适用于 Windows 的预编译可执行文件压缩包：
> [**下载 OilMaggot 1.0 执行程序 (.zip)**](https://github.com/ranwuge/OilMaggot1.0/releases/latest)

---

## 🛠 Tech Stack / 技术栈
* **GUI**: PySide6 (Qt for Python)
* **Chemistry logic**: RDKit, PubChemPy
* **Data Process**: Pandas, Openpyxl
* **Plotting**: Matplotlib
