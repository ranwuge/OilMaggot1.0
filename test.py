import sys
import time
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import PySide6.QtGui
from PySide6.QtCore import Slot, Signal, QThread
from PySide6.QtWidgets import QPushButton, QWidget, QVBoxLayout, QLabel, \
    QApplication, QMainWindow, QHBoxLayout, QSizePolicy, QFrame, QFileDialog
from PySide6.QtGui import QPixmap
from matplotlib.figure import Figure
import SlotFunction
import functions

'''Define the UI elements in this part of py file.'''


class MPLCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)


class DoSomething(QThread):
    def __init__(self):
        super().__init__()
        self.new_label = None

    def run(self):
        time.sleep(10)
        self.new_label.setText('Over????')


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.thread = None
        self.button_test = QPushButton('Push')
        self.button_test.setParent(self)
        self.button_test.clicked.connect(self.doSomething)
        self.canvas = MPLCanvas()
        self.pic = self.canvas.fig.add_subplot(111)
        x = range(2, 26, 2)
        y = [15, 13, 14.5, 17, 20, 25, 26, 26, 27, 22, 18, 15]
        self.pic.plot(x, y, color='red')
        self.new_label = QLabel('Fuck!!!!')
        self.new_label.setParent(self)

    @Slot()
    def doSomething(self):
        self.new_label.setText('Start????')
        self.thread = DoSomething()
        self.thread.start()


if __name__ == '__main__':
    app = QApplication()
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec())
