import sys
import os
import tempfile
import shutil
import logging
import time
import platform
import webbrowser

import threads

from PyQt4 import QtCore, QtGui, QtDeclarative
from blaster_ui import Ui_MainWindow

import database_tools
import sequence_tools
import db_wizard

# pyrcc4 blaster.qrc > blaster_rc.py
# pyuic4 blaster.ui > blaster_ui.py

class MyPopup(QtGui.QWidget):
    def __init__(self, images_location):
        QtGui.QWidget.__init__(self)

        self.setWindowTitle('BlasterQt v. 1.1.1')
        self.images_location = images_location
        label = QtGui.QLabel(self)
        pixmap = QtGui.QPixmap(self.images_location + "about.png")
        label.setPixmap(pixmap)


class MyMainWindow(QtGui.QMainWindow):
    def __init__(self, *args):
        QtGui.QMainWindow.__init__(self, *args)

        self.ui = Ui_MainWindow()

        # Main paths
        self.data_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.DataLocation))
        self.data_location = self.data_location.split('Local')[0] + '/Local/'
        self.home_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.HomeLocation))
        self.temp_location = str(QtGui.QDesktopServices.storageLocation(QtGui.QDesktopServices.TempLocation))
        self.app_location = self.data_location + '/BlasterQt/'
        self.db_location = self.app_location + '/databases/'
        self.images_location = self.app_location + '/images/'

        #rint self.temp_location
        #print self.data_location.split('Local')[0] + '/Local/'

        try:
            self.create_folders()
            self.copying_files()
            # Logging
            self.log_file = self.app_location + '/logging_blaster.txt'
            logging.basicConfig(filename=self.log_file, level=logging.DEBUG)

        except (IOError, OSError):
            logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
            logging.debug(str(platform.system()+platform.release()))
            logging.exception('Got exception on main handler')
            raise

        self.ui.setupUi(self)
        self.change_method_type()
        self.create_connects()
        self.ui.plainTextEdit_2.setFocus()
        self.get_all_dbs()


        # Hide some widgets (might be used later)
        self.ui.label_21.hide()
        self.ui.comboBox_2.hide()
        self.ui.spinBox.hide()
        self.ui.statusbar.showMessage("System Status | Normal | " + self.app_location)

    def create_connects(self):
        """Connect events."""

        # Tab Similarity Search
        self.ui.pushButton.clicked.connect(self.open_file_and_insert_seq)
        self.ui.comboBox_6.currentIndexChanged['QString'].connect(self.change_method_type)
        self.ui.pushButton_3.clicked.connect(self.start_similarity_search)
        #self.ui.plainTextEdit_2.textChanged.connect(self.clear_labels)

        # Tab Align
        self.ui.pushButton_10.clicked.connect(self.open_file_and_insert_seq)
        self.ui.pushButton_11.clicked.connect(self.start_alignment)
        
        # Tab Bowtie
        #self.ui.pushButton_2.clicked.connect(self.open_file_and_insert_seq)
        #self.ui.pushButton_4.clicked.connect(self.start_bowtie)
        
        # Tab RC
        self.ui.pushButton_8.clicked.connect(self.open_file_and_insert_seq)
        self.ui.pushButton_6.clicked.connect(self.start_sequence_converter)
        
        # Tab Extract
        self.ui.pushButton_5.clicked.connect(self.extract_ids)
        
        # Tab DB
        self.ui.commandLinkButton_2.clicked.connect(self.delete_dbs)
        self.ui.commandLinkButton.clicked.connect(self.start_db_wizard)
        self.ui.tableWidget.setColumnWidth(3, 200)

        # TB Converter
        self.ui.pushButton_9.clicked.connect(self.open_file_and_insert_seq)
        self.ui.pushButton_7.clicked.connect(self.start_format_convert)

        # Menu
        self.ui.actionAbout.triggered.connect(self.show_about_message)
        self.ui.actionDocumentation.triggered.connect(self.show_help)
        self.ui.actionClear_all.triggered.connect(self.clear_all_fields)



    def create_folders(self):
        """Create important folders of Blaster in the data path."""
        location_folders = [self.app_location, self.db_location, self.images_location, self.db_location + '/Fasta']
        for folder in location_folders:
            if not os.path.exists(folder):
                os.mkdir(folder)

    def copying_files(self):
        """Copy important files of Blaster in the data path."""
        to_copy_folders = os.listdir(os.getcwd() + '/to_copy/')
        for folder in to_copy_folders:
            to_copy_files = os.listdir(os.getcwd() + '/to_copy/' + folder)
            for files in to_copy_files:
                if not folder.startswith('.'):
                    if not files.startswith('.'):
                        if not os.path.exists(self.app_location + '/' + folder + '/' + files):
                            shutil.copyfile(os.getcwd() + '/to_copy/' + folder + '/' + files, self.app_location + '/' + folder + '/' + files)

    #================================================================================================================
    ### File open management
    def open_file_and_insert_seq(self):
        """Open a file, validate the input and insert the sequence into plaintext and label."""
        # Get widgets

        sending_button = self.sender()
        file_format = u"Sequence (*.fasta *.fas *.txt)"
        if str(sending_button.objectName()) == "pushButton":
            widgets = [self.ui.label_6, self.ui.plainTextEdit_2]
        #elif str(sending_button.objectName()) == "pushButton_2":
            #widgets = [self.ui.label_11, self.ui.plainTextEdit_3]
        elif str(sending_button.objectName()) == "pushButton_8":
            widgets = [self.ui.label_3, self.ui.plainTextEdit_5]
        elif str(sending_button.objectName()) == "pushButton_9":
            widgets = [self.ui.label_15, self.ui.plainTextEdit_6]
            file_format = u"Sequence or comma separated values(*.fasta *.fas *.txt *.csv)"
        elif str(sending_button.objectName()) == "pushButton_10":
            widgets = [self.ui.label_24, self.ui.plainTextEdit_7]

        # Clear plainTextEdit
        widgets[1].clear()
        # Get path for file path

        path = str(self.open_sequence_file(file_format))
        if path != "None":
            widgets[0].setText(path)
            self.insert_seq(path, widgets[1])

    def insert_seq(self, file_location, plain_textedit):
        """Insert a sequence file content into a text edit widget."""
        f_seq = open(file_location, 'r')
        sequence = f_seq.read()
        plain_textedit.clear()
        plain_textedit.insertPlainText(sequence)
        f_seq.close()
    
    def open_sequence_file(self, file_format):
        """Open a sequence file and return the path."""
        sequence_file_location = QtGui.QFileDialog.getOpenFileName(
                                    self,
                                    u"Open a sequence file",
                                    self.home_location,
                                    file_format)
        if not sequence_file_location.isNull():
            if os.path.exists(sequence_file_location):
                return str(sequence_file_location)

    ### File save management
    def save_validate_query_sequence(self, lable, plaintext, validate_format):
        """Save and validate sequence and fasta format."""
        # Check whether a file or plain text is used.
        if os.path.exists(lable.text()):
            sequence_file_location = lable.text()
            f_query = open(sequence_file_location, 'r')
            sequence = f_query.read()
            f_query.close()
        elif str(plaintext.toPlainText()) != '':
            sequence = str(plaintext.toPlainText())
        else:
            sequence = None
            self.show_info_message('Please enter a DNA sequence or upload a file!')

        # Validate and create temp file
        if sequence is not None:
            # For sequences allowed multiple fasta format
            if validate_format == 'multi':
                validate = sequence_tools.validate_fasta_seq(sequence)
            # For sequences allowed single fasta format
            elif validate_format == 'single':
                validate = sequence_tools.validate_single_fasta_seq(sequence)
            elif validate_format == 'csv':
                validate = sequence_tools.validate_csv_seq(sequence)

            validate_plain = sequence_tools.validate_seq(sequence)

            # Create a temp file with the query sequence
            if validate_plain is True or validate is True:
                temp_seq_file = tempfile.mkstemp()
                f_seq = open(temp_seq_file[1], 'w')
                if validate_plain:
                    f_seq.write('>My_sequence' + '\n')
                f_seq.write(sequence)
                f_seq.close()
                sequence_file_location = temp_seq_file[1]
                return sequence_file_location
            else:
                self.show_info_message('Please enter a valid DNA sequence!')
            
    #================================================================================================================
    ### Database stuff
    def get_all_dbs(self):
        """Get all DBs and update Combobox."""
        #self.ui.comboBox_12.clear()
        db_dict, table_dict = database_tools.get_all_db(str(self.db_location))
        db_type = str(self.ui.comboBox_6.currentText())
    
        for db_data in db_dict.iteritems():
            if db_data[0] == 'BOWTIE':
                self.ui.comboBox_7.clear()
                self.ui.comboBox_7.addItems(db_data[1])
            if db_data[0] == 'FASTA':
                self.ui.comboBox_14.clear()
                self.ui.comboBox_14.addItems(db_data[1])
            if db_data[0] == 'BLASTN':
                if db_type == 'blastn' or db_type == 'tblastn' or db_type == 'tblastx':
                    self.ui.comboBox_7.clear()
                    self.ui.comboBox_7.addItems(db_data[1])
            if db_data[0] == 'BLASTP':
                if db_type == 'blastp' or db_type == 'blastx':
                    self.ui.comboBox_7.clear()
                    self.ui.comboBox_7.addItems(db_data[1])
        
        # Table of databases   
        self.ui.tableWidget.clear()       
        for values in table_dict.values():
            max_rows = len(values)
        self.ui.tableWidget.setRowCount(max_rows)
        self.ui.tableWidget.setColumnCount(4) 
        
        horHeaders = ['Database type', 'Database name', 'Database size (KB)', 'Created'] 
        for n, key in enumerate(table_dict.keys()): 
            horHeaders.append(key) 
            for m, item in enumerate(table_dict[key]):
                newitem = QtGui.QTableWidgetItem(str(item)) 
                if n == 0:
                    newitem.setCheckState(QtCore.Qt.Unchecked)
                self.ui.tableWidget.setItem(m, n, newitem) 
        self.ui.tableWidget.setHorizontalHeaderLabels(horHeaders)
        self.ui.tableWidget.sortItems(0)    
        
    def delete_dbs(self):
        """Get all selected DBs for deleting."""
        row_count = self.ui.tableWidget.rowCount()
        column_count = self.ui.tableWidget.rowCount()
        db_to_delete = []
        for i in range(row_count):
            for x in range(column_count):
                item = self.ui.tableWidget.item(i, x)
                if item is not None:
                    if item.checkState() == QtCore.Qt.Checked:
                        checked_type = self.ui.tableWidget.item(i, x).text()
                        checked_name = self.ui.tableWidget.item(i, x+1).text()
                        if not (checked_type, checked_name) in db_to_delete:
                            db_to_delete.append((str(checked_type), str(checked_name)))
        
        if db_to_delete:
            delete = QtGui.QMessageBox.question(self, u"Delete database?",
                                                u"Are you sure that you want to delete the database(s)?",
                                                QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if delete == QtGui.QMessageBox.Yes:
                info_message = database_tools.delete_databases(db_to_delete, self.db_location)
                time.sleep(1)
                self.show_info_message(info_message[0])
        else:
            self.show_info_message("Please select a database!")
        # Update DB table
        self.get_all_dbs()


    def start_db_wizard(self):
        """Starts DB wizard."""
        self.d = db_wizard.DBWizard(self)
        self.d.show()
                
    #===================================================================================================================
    ### Program start
    def start_format_convert(self):
        """Converts a fasta to excel sequence file or vice versa amd multi/single fasta."""
        if self.ui.radioButton_3.isChecked():
            sequence_tmp = self.save_validate_query_sequence(self.ui.label_15, self.ui.plainTextEdit_6, 'multi')
        elif self.ui.radioButton_4.isChecked():
            sequence_tmp = self.save_validate_query_sequence(self.ui.label_15, self.ui.plainTextEdit_6, 'csv')
        elif self.ui.radioButton_6.isChecked():
            file_directory = QtGui.QFileDialog.getExistingDirectory(self, 'Select a directory')
            sequence_tmp = None
            sequence_tools.convert_multi_to_single(file_directory)
        elif self.ui.radioButton_5.isChecked():
            file_directory = QtGui.QFileDialog.getExistingDirectory(self, 'Please select a directory')
            sequence_tmp = None
            sequence_tools.convert_single_to_multi(file_directory)

        if sequence_tmp is not None:
            self.ui.statusbar.showMessage("Starting Conversation | Please wait...")
            if self.ui.radioButton_3.isChecked():
                converted_file = sequence_tools.convert_seq_fomats(sequence_tmp, 'fasta_to_excel')
                webbrowser.open(converted_file)
            else:
                converted_file = sequence_tools.convert_seq_fomats(sequence_tmp, 'excel_to_fasta')
                webbrowser.open(converted_file)
        elif file_directory is not None:
            webbrowser.open(file_directory)
        self.ui.statusbar.showMessage("Conversation Done")

    def start_sequence_converter(self):
        """Converts a sequence."""
        self.setCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        self.ui.statusbar.showMessage("Starting to convert sequences | Please wait...")
        sequence_tmp = self.save_validate_query_sequence(self.ui.label_3, self.ui.plainTextEdit_5, 'multi')

        if sequence_tmp is not None:
            if self.ui.radioButton_2.isChecked():
                sequence_out = sequence_tools.reverse_complement(sequence_tmp, True)
            elif self.ui.radioButton.isChecked():
                sequence_out = sequence_tools.reverse_complement(sequence_tmp, False)
            elif self.ui.radioButton_7.isChecked():
                sequence_out = sequence_tools.translate_sequence(sequence_tmp, False)
            webbrowser.open(sequence_out)

        self.ui.statusbar.showMessage("Done")
        self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))

    def start_alignment(self):
        """Start a ClustalOmega alignment."""
        self.setCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        self.ui.statusbar.showMessage("Starting Alignment | Please wait...")
        sequence_tmp = self.save_validate_query_sequence(self.ui.label_24, self.ui.plainTextEdit_7, 'multi')

        if sequence_tmp is not None:
            sequence_tmp = sequence_tools.clustal_omega(sequence_tmp, self.db_location)
            webbrowser.open(sequence_tmp)

        self.ui.statusbar.showMessage("Alignment done")
        self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        
    def start_bowtie(self):
        """Starts Bowtie alignment."""
        if self.ui.comboBox_7.currentText() != "":
            self.sequence_tmp = self.save_validate_query_sequence(self.ui.label_6, self.ui.plainTextEdit_2, 'multi')

            # Check weather sequences are too long
            if not sequence_tools.validate_bowtie_seq(self.sequence_tmp):
                if self.sequence_tmp is not None:
                    self.setCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
                    self.ui.statusbar.showMessage("Starting BOWTIE | Please wait...")
                    self.searchThread = threads.BowtieThread(str(self.sequence_tmp),
                                                                     str(self.db_location),
                                                                     str(self.db_location + self.ui.comboBox_7.currentText()),
                                                                     str(self.ui.comboBox_13.currentText()))
                    self.connect(self.searchThread, QtCore.SIGNAL("threadDone(QString, QString)"), self.thread_bowtie_done)
                    self.searchThread.start()
            else:
                self.show_info_message('Sequence too long (max. 1000 bases)!')
        else:
            self.show_info_message('Please choose a database.')
        self.ui.label_6.setText('')

        # if self.ui.comboBox_12.currentText() != "":
        #     self.sequence_tmp = self.save_validate_query_sequence(self.ui.label_11, self.ui.plainTextEdit_3, 'multi')
        #     if self.sequence_tmp is not None:
        #         self.setCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
        #         self.ui.statusbar.showMessage("Starting BOWTIE | Please wait...")
        #         self.searchThread = threads.BowtieThread(str(self.sequence_tmp),
        #                                                          str(self.db_location),
        #                                                          str(self.db_location + self.ui.comboBox_12.currentText()),
        #                                                          str(self.ui.comboBox_11.currentText()))
        #         self.connect(self.searchThread, QtCore.SIGNAL("threadDone(QString, QString)"), self.thread_bowtie_done)
        #         self.searchThread.start()
        # else:
        #     self.show_info_message('Please choose a database.')
        # self.ui.label_11.setText('')

    def thread_bowtie_done(self, info_message, bowtie_file_path):
        """Show message after thread has finished."""
        if os.path.exists(bowtie_file_path):
                statinfo = os.stat(bowtie_file_path)

                if statinfo.st_size > 0:
                    # Add header
                    final_bowtie_file_path = sequence_tools.add_header(bowtie_file_path, 'bowtie')
                    webbrowser.open(final_bowtie_file_path)
                    self.ui.statusbar.showMessage("BOWTIE search done")
                else:
                    self.ui.statusbar.showMessage("BOWTIE search failed")
                    self.show_info_message('Sorry, no matches found.') 
        self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
                    
    def start_similarity_search(self):
        """Starts BLAST alignment."""

        if str(self.ui.comboBox_6.currentText()) == 'Bowtie':
            self.start_bowtie()
        else:

            settings = [str(self.ui.lineEdit_6.text()),
                     str(self.ui.lineEdit.text()),
                     '',
                     str(self.ui.lineEdit_7.text()),
                     str(self.ui.comboBox.currentText()),
                     str(self.ui.lineEdit_8.text()),
                     str(self.ui.comboBox_9.currentText()).split(' ')[0],
                     str(self.ui.comboBox_10.currentText()),
                     str(self.ui.lineEdit_9.text()),
                     str(self.ui.lineEdit_10.text()),
                     str(self.ui.checkBox.checkState())]

            if self.ui.comboBox_7.currentText() != "":
                sequence_tmp = self.save_validate_query_sequence(self.ui.label_6, self.ui.plainTextEdit_2, 'multi')
                if sequence_tmp is not None:
                    self.ui.statusbar.showMessage("Starting BLAST | Please wait...")
                    self.setCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
                    self.searchThreadBlast = threads.BlastThread(str(sequence_tmp),
                                                                     str(self.db_location),
                                                                     str(self.ui.comboBox_7.currentText()),
                                                                     str(self.ui.comboBox_6.currentText()),
                                                                     settings)
                    self.connect(self.searchThreadBlast, QtCore.SIGNAL("threadDone(QString, QString)"), self.thread_similarity_search_done)
                    self.searchThreadBlast.start()
            else:
                self.show_info_message('Please choose a database.')
            self.ui.label_6.setText('')

    def thread_similarity_search_done(self, info_message, blast_file_path):
        """Show message after thread has finished."""
        if os.path.exists(blast_file_path):
                statinfo = os.stat(blast_file_path)
                if statinfo.st_size > 0:
                    webbrowser.open(blast_file_path)
                    self.ui.statusbar.showMessage("BLAST done")
                else:
                    self.ui.statusbar.showMessage("BLAST failed")
                    self.show_info_message('Sorry, no matches found.')  
        self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
                    
    def extract_ids(self):
        """Extract sequences from IDs."""
        if self.ui.comboBox_14.currentText() != "":
            self.setCursor(QtGui.QCursor(QtCore.Qt.WaitCursor))
            self.ui.statusbar.showMessage("Extracting sequences | Please wait...")

            ids = str(self.ui.plainTextEdit_4.toPlainText())
            ignore_chars = int(self.ui.spinBox.value())
            ignore_pos = str(self.ui.comboBox_2.currentText())
            final_ids = []
            for id in ids.split('\n'):
                if ignore_chars == 0:
                    id = id
                else:
                    if ignore_pos == 'Beginning':
                        id = id[ignore_chars:]
                    elif ignore_pos == 'End':
                        id = id[0:(len(id)-ignore_chars)]
                if id != '':
                    final_ids.append(id)

            if len(final_ids) > 0:
                extracted_tmp_file = sequence_tools.id_extract(final_ids, str(self.db_location), str(self.ui.comboBox_14.currentText())+'.fasta')
                if os.path.exists(extracted_tmp_file):
                    statinfo = os.stat(extracted_tmp_file)
                    if statinfo.st_size > 0:
                        webbrowser.open(extracted_tmp_file)
                        self.ui.statusbar.showMessage("Sequences extracted.")
                    else:
                        self.ui.statusbar.showMessage("Sequence extraction failed.")
                        self.show_info_message('Sorry, no matches found.')
            else:
                self.show_info_message('Please enter an ID.')
            self.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        else:
            self.show_info_message('Please choose a database.')
        
    #================================================================================================================
    ### BLAST settings
    def change_method_type(self):
        """Sets the similarity search method type settings (Blast and Bowtie)."""

        method_type = self.ui.comboBox_6.currentText()
        settings = sequence_tools.set_method_settings(method_type)
        
        # evalue
        if settings[0] == 'enabled':
            self.ui.lineEdit_6.setEnabled(False)
        else:
            self.ui.lineEdit_6.setText(settings[0])
        # min_length
        #self.ui.lineEdit_2.setText(settings[2])
        # wordsize
        if settings[3] == 'enabled':
            self.ui.lineEdit_7.setEnabled(False)
        else:
            self.ui.lineEdit_7.setText(settings[3])
        # processor
        if settings[5] == 'enabled':
             self.ui.lineEdit_8.setEnabled(False)
        else:
            self.ui.lineEdit_8.setText(settings[5])
        # min_percent
        if settings[1] == 'enabled':
            self.ui.lineEdit.setEnabled(False)
        else:
            self.ui.lineEdit.setEnabled(True)
            self.ui.lineEdit.setText(settings[1]) 
        # strand
        if settings[4] == 'enabled':
            self.ui.comboBox.setEnabled(False)
        else:
            self.ui.comboBox.setEnabled(True)
            self.ui.comboBox.setCurrentIndex(settings[4])
        # genetic_code
        if settings[6] == 'enabled':
            self.ui.comboBox_9.setEnabled(False)
        else:
            self.ui.comboBox_9.setEnabled(True)
            self.ui.comboBox_9.setCurrentIndex(settings[6])
        # matrix
        if settings[7] == 'enabled':
            self.ui.comboBox_10.setEnabled(False)
        else:
            self.ui.comboBox_10.setEnabled(True)
            self.ui.comboBox_10.setCurrentIndex(settings[7])
        # gap_open
        if settings[8] == 'enabled':
            self.ui.lineEdit_9.setEnabled(False)
        else:
            self.ui.lineEdit_9.setEnabled(True)
            self.ui.lineEdit_9.setText(settings[8])
        # gap_ext
        if settings[9] == 'enabled':
            self.ui.lineEdit_10.setEnabled(False)
        else:
            self.ui.lineEdit_10.setEnabled(True)
            self.ui.lineEdit_10.setText(settings[9])
        # ungapped
        if settings[10] == 'enabled':
            self.ui.checkBox.setEnabled(False)
        else:
            self.ui.checkBox.setEnabled(True)
            self.ui.checkBox.setCheckState(settings[10])
        # Mismmatches for Bowtie
        if settings[11] == 'enabled':
            self.ui.comboBox_13.setEnabled(False)
        else:
            self.ui.comboBox_13.setEnabled(True)

        self.get_all_dbs()

    def clear_all_fields(self):
        """Clears all sequence and ID fields."""
        fields = [self.ui.plainTextEdit_2, self.ui.plainTextEdit_4,
                  self.ui.plainTextEdit_5, self.ui.plainTextEdit_6]

        for field in fields:
            field.clear()

    def show_info_message(self, message):
        QtGui.QMessageBox.information(self,
                    u"Information",
                    message)

    def show_about_message(self):
        """Shows about box."""
        self.w = MyPopup(self.images_location)
        self.w.setGeometry(QtCore.QRect(100, 100, 640, 350))
        self.w.show()


    def show_help(self):
        """Show help."""
        webbrowser.open("labtools.ipk-gatersleben.de/help/BlasterQt/BlasterQt.html")


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    QtGui.QApplication.setStyle(QtGui.QStyleFactory.create("gtk"))
    QtGui.QApplication.setPalette(QtGui.QApplication.style().standardPalette())
    myapp = MyMainWindow()
    myapp.show()
    sys.exit(app.exec_())



   

