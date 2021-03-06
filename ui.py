import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from beam_solver import beam_solver


class Ui_MainWindow(object):

    def __init__(self):
        self.point_loads = np.empty((0,3), float)
        self.point_moments = np.empty((0,2), float)
        self.linear_loads = np.empty((0,3), float)

        self.L = 0
        self.Xa = 0
        self.Xb = 0


    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")

        # ______________________INPUT FORCES______________________
        self.label_forces = QtWidgets.QLabel(self.centralwidget)
        self.label_forces.setGeometry(QtCore.QRect(50, 260, 291, 61))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_forces.setFont(font)
        self.label_forces.setObjectName("label_forces")

        font = QtGui.QFont()
        font.setPointSize(10)

        self.label_force_x = QtWidgets.QLabel(self.centralwidget)
        self.label_force_x.setGeometry(QtCore.QRect(50, 320, 50, 20))
        self.label_force_x.setObjectName("label_force_x")

        self.Force_x = QtWidgets.QLineEdit(self.centralwidget)
        self.Force_x.setGeometry(QtCore.QRect(120, 320, 70, 20))
        self.Force_x.setObjectName("Force_x")

        self.label_force_magnitude = QtWidgets.QLabel(self.centralwidget)
        self.label_force_magnitude.setGeometry(QtCore.QRect(50, 420, 50, 20))
        self.label_force_magnitude.setObjectName("label_force_magnitude")

        self.Force_magnitude = QtWidgets.QLineEdit(self.centralwidget)
        self.Force_magnitude.setGeometry(QtCore.QRect(120, 420, 70, 20))
        self.Force_magnitude.setObjectName("Force_magnitude")

        self.buttonForce = QtWidgets.QPushButton(self.centralwidget)
        self.buttonForce.setGeometry(QtCore.QRect(85, 470, 120, 30))
        self.buttonForce.setObjectName("buttonForce")
        self.buttonForce.clicked.connect(self.get_loads)


        #______________________INPUT MOMENTS______________________
        self.label_moments = QtWidgets.QLabel(self.centralwidget)
        self.label_moments.setGeometry(QtCore.QRect(210, 260, 291, 61))
        self.label_moments.setFont(font)
        self.label_moments.setObjectName("label_moments")
        font = QtGui.QFont()
        font.setPointSize(10)

        self.label_moment_x = QtWidgets.QLabel(self.centralwidget)
        self.label_moment_x.setGeometry(QtCore.QRect(210, 320, 50, 20))
        self.label_moment_x.setObjectName("label_moment_x")

        self.Moment_x = QtWidgets.QLineEdit(self.centralwidget)
        self.Moment_x.setGeometry(QtCore.QRect(280, 320, 70, 20))
        self.Moment_x.setObjectName("Moment_x")

        self.label_moment_magnitude = QtWidgets.QLabel(self.centralwidget)
        self.label_moment_magnitude.setGeometry(QtCore.QRect(210, 420, 65, 20))
        self.label_moment_magnitude.setObjectName("label_moment_magnitude")

        self.Moment_magnitude = QtWidgets.QLineEdit(self.centralwidget)
        self.Moment_magnitude.setGeometry(QtCore.QRect(280, 420, 70, 20))
        self.Moment_magnitude.setObjectName("Moment_magnitude")

        self.buttonMoment = QtWidgets.QPushButton(self.centralwidget)
        self.buttonMoment.setGeometry(QtCore.QRect(245, 470, 120, 30))
        self.buttonMoment.setObjectName("buttonMoment")
        self.buttonMoment.clicked.connect(self.get_moments)


        # ______________________INPUT LINEAR LOADS______________________
        self.label_linear_loads = QtWidgets.QLabel(self.centralwidget)
        self.label_linear_loads.setGeometry(QtCore.QRect(370, 260, 291, 61))
        self.label_linear_loads.setFont(font)
        self.label_linear_loads.setObjectName("label_linear_loads")

        self.label_linear_load_start = QtWidgets.QLabel(self.centralwidget)
        self.label_linear_load_start.setGeometry(QtCore.QRect(370, 320, 90, 20))
        self.label_linear_load_start.setObjectName("label_linear_load_start")

        self.LinearLoad_start_x = QtWidgets.QLineEdit(self.centralwidget)
        self.LinearLoad_start_x.setGeometry(QtCore.QRect(470, 320, 70, 20))
        self.LinearLoad_start_x.setObjectName("LinearLoad_start_x")

        self.label_linear_load_stop = QtWidgets.QLabel(self.centralwidget)
        self.label_linear_load_stop.setGeometry(QtCore.QRect(370, 370, 80, 20))
        self.label_linear_load_stop.setObjectName("label_linear_load_stop")

        self.LinearLoad_end_x = QtWidgets.QLineEdit(self.centralwidget)
        self.LinearLoad_end_x.setGeometry(QtCore.QRect(470, 370, 70, 20))
        self.LinearLoad_end_x.setObjectName("LinearLoad_end_x")

        self.label_linear_load_magnitude = QtWidgets.QLabel(self.centralwidget)
        self.label_linear_load_magnitude.setGeometry(QtCore.QRect(370, 420, 70, 20))
        self.label_linear_load_magnitude.setObjectName("label_linear_load_magnitude")

        self.LinearLoad_magnitude = QtWidgets.QLineEdit(self.centralwidget)
        self.LinearLoad_magnitude.setGeometry(QtCore.QRect(470, 420, 70, 20))
        self.LinearLoad_magnitude.setObjectName("LinearLoad_magnitude")

        self.buttonLinearLoad = QtWidgets.QPushButton(self.centralwidget)
        self.buttonLinearLoad.setGeometry(QtCore.QRect(405, 470, 120, 30))
        self.buttonLinearLoad.setObjectName("buttonLinearLoad")
        self.buttonLinearLoad.clicked.connect(self.get_linear_loads)


        # ______________________INPUT DIMENSIONS______________________
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(30, 40, 551, 211))
        self.label.setObjectName("label")

        self.label_L = QtWidgets.QLabel(self.centralwidget)
        self.label_L.setGeometry(QtCore.QRect(600, 70, 70, 20))
        self.label_L.setObjectName("label_L")

        self.lineEdit_L = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_L.setGeometry(QtCore.QRect(680, 70, 70, 20))
        self.lineEdit_L.setObjectName("lineEdit_L")

        self.label_Xa = QtWidgets.QLabel(self.centralwidget)
        self.label_Xa.setGeometry(QtCore.QRect(600, 110, 70, 20))
        self.label_Xa.setObjectName("label_Xa")

        self.lineEdit_Xa = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_Xa.setGeometry(QtCore.QRect(680, 110, 70, 20))
        self.lineEdit_Xa.setObjectName("lineEdit_Xa")

        self.label_Xb = QtWidgets.QLabel(self.centralwidget)
        self.label_Xb.setGeometry(QtCore.QRect(600, 150, 70, 20))
        self.label_Xb.setObjectName("label_Xb")

        self.lineEdit_Xb = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_Xb.setGeometry(QtCore.QRect(680, 150, 70, 20))
        self.lineEdit_Xb.setObjectName("lineEdit_Xb")

        self.buttonInputDimensions = QtWidgets.QPushButton(self.centralwidget)
        self.buttonInputDimensions.setGeometry(QtCore.QRect(620, 200, 120, 30))
        self.buttonInputDimensions.setObjectName("buttonLinearLoad")
        self.buttonInputDimensions.clicked.connect(self.get_dimensions)

        # ______________________INPUT Jz______________________
        self.label_Jz = QtWidgets.QLabel(self.centralwidget)
        self.label_Jz.setGeometry(QtCore.QRect(600, 260, 70, 20))
        self.label_Jz.setObjectName("label_Jz")

        self.lineEdit_Jz = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_Jz.setGeometry(QtCore.QRect(680, 260, 70, 20))
        self.lineEdit_Jz.setObjectName("lineEdit_Jz")

        self.buttonInputJz = QtWidgets.QPushButton(self.centralwidget)
        self.buttonInputJz.setGeometry(QtCore.QRect(620, 300, 120, 30))
        self.buttonInputJz.setObjectName("buttonLinearLoad")
        self.buttonInputJz.clicked.connect(self.get_Jz)

        # ______________________SOLVE BUTTON______________________
        self.buttonSolve = QtWidgets.QPushButton(self.centralwidget)
        self.buttonSolve.setGeometry(QtCore.QRect(300, 520, 160, 50))
        self.buttonSolve.setObjectName("buttonSolve")
        self.buttonSolve.clicked.connect(self.run_solver)

        # ______________________RESET BUTTON______________________
        self.buttonReset = QtWidgets.QPushButton(self.centralwidget)
        self.buttonReset.setGeometry(QtCore.QRect(100, 520, 160, 50))
        self.buttonReset.setObjectName("buttonSolve")
        self.buttonReset.clicked.connect(self.reset_loads)

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Beam Solver"))
        self.label_forces.setText(_translate("MainWindow", "Add force"))
        self.label_moments.setText(_translate("MainWindow", "Add Moment"))
        self.label_linear_loads.setText(_translate("MainWindow", "Add linear load"))
        self.label_force_x.setText(_translate("MainWindow", "x [m]"))
        self.label_force_magnitude.setText(_translate("MainWindow", "F [kN]"))
        self.label_moment_x.setText(_translate("MainWindow", "x [m]"))
        self.label_moment_magnitude.setText(_translate("MainWindow", "M [kNm]"))
        self.label_linear_load_start.setText(_translate("MainWindow", "x_start [m]"))
        self.label_linear_load_stop.setText(_translate("MainWindow", "x_end [m]"))
        self.label_linear_load_magnitude.setText(_translate("MainWindow", "q [kN/m]"))

        self.label_Jz.setText(_translate("MainWindow", "Jz [m4]"))


        self.buttonForce.setText(_translate("MainWindow", "Add force"))
        self.buttonMoment.setText(_translate("MainWindow", "Add moment"))
        self.buttonLinearLoad.setText(_translate("MainWindow", "Add linear load"))
        self.buttonSolve.setText(_translate("MainWindow", "Solve!"))
        self.buttonReset.setText(_translate("MainWindow", "Reset loads"))
        self.buttonInputDimensions.setText(_translate("MainWindow", "Set dimensions"))
        self.buttonInputJz.setText(_translate("MainWindow", "Set Jz"))
        self.label.setText(_translate("MainWindow", "<html><head/><body><p><img src=\"src/beam.png\"/></p></body></html>"))
        self.label_L.setText(_translate("MainWindow", "L [m]"))
        self.label_Xa.setText(_translate("MainWindow", "Xa [m]"))
        self.label_Xb.setText(_translate("MainWindow", "Xb [m]"))

    def get_loads(self):

        if self.Force_x.text() != '' and self.Force_magnitude.text() != '':
            try:
                x = float(self.Force_x.text())
                if x >= 0 and x <= self.L:
                    y = 0
                    F = float(self.Force_magnitude.text())
                    pointLoad = np.array([[x, y, F]])
                    self.point_loads = np.append(self.point_loads, pointLoad, axis=0)
                    print("point_loads = ")
                    print(self.point_loads)
                    print("point_moments = ")
                    print(self.point_moments)
                    print("linear_loads = ")
                    print((self.linear_loads))

                    self.Force_x.clear()
                    self.Force_magnitude.clear()
                else:
                    print("Force outside of beam!")
            except:
                print("Input not a number")
        else:
            print("Values missing!")

    def get_moments(self):

        if self.Moment_x.text() != '' and self.Moment_magnitude.text() != '':
            try:
                x = float(self.Moment_x.text())
                if x >= 0 and x <= self.L:
                    M = float(self.Moment_magnitude.text())
                    pointMoment = np.array([[x, M]])
                    self.point_moments = np.append(self.point_moments, pointMoment, axis=0)
                    print("point_loads = ")
                    print(self.point_loads)
                    print("point_moments = ")
                    print(self.point_moments)
                    print("linear_loads = ")
                    print((self.linear_loads))

                    self.Moment_x.clear()
                    self.Moment_magnitude.clear()
                else:
                    print("Moment outside of beam!")
            except:
                print("Input not a number")
        else:
            print("Values missing!")

    def get_linear_loads(self):

        if self.LinearLoad_start_x.text() != '' and self.LinearLoad_end_x.text() != '' and self.LinearLoad_magnitude.text() != '':
            try:
                x_start = float(self.LinearLoad_start_x.text())
                x_end = float(self.LinearLoad_end_x.text())

                if x_start >= 0 and x_start <= self.L and x_end >= 0 and x_end <= self.L:
                    q = float(self.LinearLoad_magnitude.text())
                    if x_start < x_end:
                        linearLoad = np.array([[x_start, x_end, q]])
                        self.linear_loads = np.append(self.linear_loads, linearLoad, axis=0)

                        print("point_loads = ")
                        print(self.point_loads)
                        print("point_moments = ")
                        print(self.point_moments)
                        print("linear_loads = ")
                        print((self.linear_loads))

                        self.LinearLoad_start_x.clear()
                        self.LinearLoad_end_x.clear()
                        self.LinearLoad_magnitude.clear()
                    else:
                        print("x_start > x_end !")
                else:
                    print("Linear force outside of beam!")
            except:
                print("Input not a number")
        else:
            print("Values missing!")

    def get_dimensions(self):

        if self.lineEdit_L.text() != '' and self.lineEdit_Xa.text() != '' and self.lineEdit_Xb.text() != '':
            try:
                L = float(self.lineEdit_L.text())
                Xa = float(self.lineEdit_Xa.text())
                Xb = float(self.lineEdit_Xb.text())

                print("L = ", L)
                print("Xa = ", Xa)
                print("Xb = ", Xb)

                if Xb < Xa:
                    print("Support dimensions invalid, dimensions not set!")
                elif Xb > L:
                    print("Beam lenght invalid, dimensions not set!")
                elif Xa < 0 or Xb < 0 or L <= 0:
                    print("Value < 0, dimensions not set!")
                else:
                    self.L = L
                    self.Xa = Xa
                    self.Xb = Xb
                    print("Dimensions set!")

                    self.lineEdit_L.clear()
                    self.lineEdit_Xa.clear()
                    self.lineEdit_Xb.clear()
            except:
                print("Input not a number")
        else:
            print("Values missing!")

    def get_Jz(self):

        if self.lineEdit_Jz.text() != '':
            try:
                Jz = float(self.lineEdit_Jz.text())
                self.Jz = Jz

                print("Jz = ", Jz)
                self.lineEdit_Jz.clear()

            except:
                print("Input not a number")
        else:
            print("Values missing!")


    def reset_loads(self):
        self.point_loads = np.empty((0, 3), float)
        self.point_moments = np.empty((0, 2), float)
        self.linear_loads = np.empty((0, 3), float)

        print("point_loads = ")
        print(self.point_loads)
        print("point_moments = ")
        print(self.point_moments)
        print("linear_loads = ")
        print((self.linear_loads))

    def run_solver(self):
        if self.Xa > self.L or self.Xb > self.L:
            print("Support outside of beam!")
        elif len(self.point_loads) == 0 and len(self.point_moments) == 0 and len(self.linear_loads) == 0:
            print("Loads missing!")
        elif (self.Jz is None):
            print("Jz missing")
        else:
            beam_solver(self.point_loads, self.point_moments, self.linear_loads, self.L, self.Xa, self.Xb, self.Jz)