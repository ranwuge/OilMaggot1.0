import functions
import tkinter
from tkinter import filedialog
if __name__ == '__main__':
    fp = filedialog.askopenfilename()
    with open(fp, 'r') as input_file:
        input_file_content = input_file.readlines()
    line_count = sum(1 for line in input_file_content)
    for i in input_file_content:
        print(i)