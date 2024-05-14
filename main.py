import sys
from pathlib import Path

import PySide6.QtWidgets
from PySide6.QtCore import QObject, Slot
from PySide6.QtGui import QGuiApplication
from PySide6.QtQml import QQmlApplicationEngine, QmlElement
from PySide6.QtQuickControls2 import QQuickStyle
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel, QApplication

'''Define the UI elements in this part of py file.'''


if __name__ == '__main__':
    app = QApplication()


    gallery = QLabel('Hello My Dick')
    gallery.setAlignment(Qt.AlignmentFlag.AlignCenter)
    gallery.setStyleSheet("""
           background-color: #262626;
           color: #FFFFFF;
           font-family: Titillium;
           font-size: 18px;
           """)
    gallery.show()

    with open("uifile/stylesheet.qss", "r") as f:
        _style = f.read()
        app.setStyleSheet(_style)

    sys.exit(app.exec())

