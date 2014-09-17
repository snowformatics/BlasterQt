# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'wizard.ui'
#
# Created: Thu Feb 13 14:26:17 2014
#      by: PyQt4 UI code generator 4.10.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_wizard(object):
    def setupUi(self, wizard):
        wizard.setObjectName(_fromUtf8("wizard"))
        wizard.resize(676, 387)
        self.wizardPage1 = QtGui.QWizardPage()
        self.wizardPage1.setObjectName(_fromUtf8("wizardPage1"))
        self.layoutWidget = QtGui.QWidget(self.wizardPage1)
        self.layoutWidget.setGeometry(QtCore.QRect(20, 20, 401, 181))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.checkBox = QtGui.QCheckBox(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.checkBox.setFont(font)
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.verticalLayout.addWidget(self.checkBox)
        self.checkBox_2 = QtGui.QCheckBox(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.checkBox_2.setFont(font)
        self.checkBox_2.setObjectName(_fromUtf8("checkBox_2"))
        self.verticalLayout.addWidget(self.checkBox_2)
        self.checkBox_3 = QtGui.QCheckBox(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.checkBox_3.setFont(font)
        self.checkBox_3.setObjectName(_fromUtf8("checkBox_3"))
        self.verticalLayout.addWidget(self.checkBox_3)
        self.checkBox_4 = QtGui.QCheckBox(self.layoutWidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.checkBox_4.setFont(font)
        self.checkBox_4.setObjectName(_fromUtf8("checkBox_4"))
        self.verticalLayout.addWidget(self.checkBox_4)
        wizard.addPage(self.wizardPage1)
        self.wizardPage_2 = QtGui.QWizardPage()
        self.wizardPage_2.setObjectName(_fromUtf8("wizardPage_2"))
        self.layoutWidget1 = QtGui.QWidget(self.wizardPage_2)
        self.layoutWidget1.setGeometry(QtCore.QRect(20, 20, 601, 141))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout_3.setMargin(0)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.commandLinkButton_5 = QtGui.QCommandLinkButton(self.layoutWidget1)
        self.commandLinkButton_5.setMaximumSize(QtCore.QSize(200, 16777215))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/Icons/Filesystem-folder-grey-open-icon.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.commandLinkButton_5.setIcon(icon)
        self.commandLinkButton_5.setIconSize(QtCore.QSize(48, 48))
        self.commandLinkButton_5.setObjectName(_fromUtf8("commandLinkButton_5"))
        self.verticalLayout_3.addWidget(self.commandLinkButton_5)
        self.label_5 = QtGui.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_5.setFont(font)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.verticalLayout_3.addWidget(self.label_5)
        wizard.addPage(self.wizardPage_2)
        self.wizardPage_3 = QtGui.QWizardPage()
        self.wizardPage_3.setObjectName(_fromUtf8("wizardPage_3"))
        self.layoutWidget2 = QtGui.QWidget(self.wizardPage_3)
        self.layoutWidget2.setGeometry(QtCore.QRect(20, 20, 401, 91))
        self.layoutWidget2.setObjectName(_fromUtf8("layoutWidget2"))
        self.verticalLayout_4 = QtGui.QVBoxLayout(self.layoutWidget2)
        self.verticalLayout_4.setMargin(0)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.label_6 = QtGui.QLabel(self.layoutWidget2)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_6.setFont(font)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.verticalLayout_4.addWidget(self.label_6)
        self.lineEdit = QtGui.QLineEdit(self.layoutWidget2)
        self.lineEdit.setMaximumSize(QtCore.QSize(200, 16777215))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.lineEdit.setFont(font)
        self.lineEdit.setObjectName(_fromUtf8("lineEdit"))
        self.verticalLayout_4.addWidget(self.lineEdit)
        wizard.addPage(self.wizardPage_3)
        self.wizardPage2 = QtGui.QWizardPage()
        self.wizardPage2.setObjectName(_fromUtf8("wizardPage2"))
        self.layoutWidget3 = QtGui.QWidget(self.wizardPage2)
        self.layoutWidget3.setGeometry(QtCore.QRect(20, 10, 351, 131))
        self.layoutWidget3.setObjectName(_fromUtf8("layoutWidget3"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.layoutWidget3)
        self.verticalLayout_2.setMargin(0)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label_2 = QtGui.QLabel(self.layoutWidget3)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.label_2.setFont(font)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_2.addWidget(self.label_2)
        self.progressBar = QtGui.QProgressBar(self.layoutWidget3)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.progressBar.setFont(font)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.verticalLayout_2.addWidget(self.progressBar)
        self.label_3 = QtGui.QLabel(self.layoutWidget3)
        self.label_3.setText(_fromUtf8(""))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_2.addWidget(self.label_3)
        wizard.addPage(self.wizardPage2)

        self.retranslateUi(wizard)
        QtCore.QMetaObject.connectSlotsByName(wizard)

    def retranslateUi(self, wizard):
        wizard.setWindowTitle(_translate("wizard", "Create new database", None))
        self.label.setText(_translate("wizard", "Please choose which database type you want to create:", None))
        self.checkBox.setText(_translate("wizard", "BLAST DNA Database", None))
        self.checkBox_2.setText(_translate("wizard", "BLAST Protein Database", None))
        self.checkBox_3.setText(_translate("wizard", "BOWTIE Database", None))
        self.checkBox_4.setText(_translate("wizard", "Fasta Database", None))
        self.commandLinkButton_5.setText(_translate("wizard", "Choose source file", None))
        self.label_5.setText(_translate("wizard", "No source selected", None))
        self.label_6.setText(_translate("wizard", "Enter a database name:", None))
        self.label_2.setText(_translate("wizard", "Progress:", None))

import blaster_rc

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    wizard = QtGui.QWizard()
    ui = Ui_wizard()
    ui.setupUi(wizard)
    wizard.show()
    sys.exit(app.exec_())

