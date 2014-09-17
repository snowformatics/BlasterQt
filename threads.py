import tempfile
import subprocess
import os
import logging
import platform
import time
from PyQt4 import QtCore, QtGui, QtDeclarative
from types import *

import database_tools

class BowtieThread(QtCore.QThread):
    def __init__(self, sequence_path, db_location, db_name, missmatches):
        QtCore.QThread.__init__(self)

        self.sequence = sequence_path
        self.db_location = db_location
        self.db_name = db_name
        self.missmatches = missmatches

    def run(self):
        """Bowtie search thread."""
        # Start Bowtie analyse
        info_message, bowtie_file_path = self.start_bowtie(self.sequence,
                                                            self.db_location,
                                                            self.db_name,
                                                            self.missmatches)
        
        self.emit(QtCore.SIGNAL("threadDone(QString, QString)"), info_message, bowtie_file_path)
        
    def start_bowtie(self, sequence, db_location, database_file_location, missmatches):
        """Start Bowtie alignment."""

        # Temp file for results
        temp_bowtie_file = tempfile.mkstemp()

        assert type(missmatches) is StringType
        assert type(database_file_location) is StringType
        assert type(sequence) is StringType
        assert type(temp_bowtie_file[1]) is StringType
    
        os.chdir(db_location)
        
        try:
            process = subprocess.Popen(["bowtie", "-a", "-v", missmatches,  "-y", database_file_location, "-f",
                                        sequence, temp_bowtie_file[1]+'.txt'])
            process.wait()

            return 'Search done!', temp_bowtie_file[1]+'.txt'
    
        except (IOError, OSError):
            logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
            logging.debug(str(platform.system()+platform.release()))
            logging.exception('Bowtie. Got exception on main handler')
            return 'Bowtie was not found! Please install <a href=\"http://bowtie-bio.sourceforge.net/index.shtml\">Bowtie</a>'
    
    
class BlastThread(QtCore.QThread):
    def __init__(self, sequence_path, db_location, db_name, blast_type, settings):
        QtCore.QThread.__init__(self)

        self.sequence = sequence_path
        self.db_location = db_location
        self.db_name = db_name
        self.blast_type = blast_type
        self.settings = settings
        
    def run(self):
        """Blast thread."""
        # Start BLAST
        info_message, blast_file_path = self.start_blast(self.sequence,
                                                            self.db_location,
                                                            self.db_name,
                                                            self.blast_type,
                                                            self.settings)
         
        self.emit(QtCore.SIGNAL("threadDone(QString, QString)"), info_message, blast_file_path)
#         
    def start_blast(self, sequence_file, db_location, db_name, blast_type, settings):
        """Start Blast alignment."""
 
        temp_blast_file = tempfile.mkstemp()
        
        assert type(settings) is ListType
        assert type(db_name) is StringType
        assert type(sequence_file) is StringType
        assert type(blast_type) is StringType
        assert type(temp_blast_file[1]) is StringType
        
        os.chdir(db_location)
        fixed = blast_type + ' -query ' + sequence_file + ' -db ' + db_name + ' -outfmt 7 -out ' + temp_blast_file[1]+'.txt '
        
        evalue = '-evalue ' + str(settings[0]) 
        min_percent = ' -perc_identity ' + str(settings[1])
        #min_length = '-evalue ' + str(settings[2] 
        wordsize = ' -word_size ' + str(settings[3]) 
        strand = ' -strand ' + str(settings[4]) 
        processor = ' -num_threads ' + str(settings[5]) 
        genetic_code = ' -query_gencode ' + str(settings[6]) 
        matrix = ' -matrix ' + str(settings[7]) 
        gap_open = ' -gapopen ' + str(settings[8]) 
        gap_ext = ' -gapextend ' + str(settings[9])
        if str(settings[10]) == '2':
            ungapped = ' -ungapped '
        else:
            ungapped = ''
        
        if blast_type == 'blastn':
            blast_cmd = fixed + evalue + min_percent + wordsize + strand + processor + gap_open + gap_ext + ungapped
        elif blast_type == 'blastx':
            blast_cmd = fixed + evalue + wordsize + strand + processor + genetic_code + matrix + gap_open + gap_ext
        elif blast_type == 'blastp':
            blast_cmd = fixed + evalue + wordsize + processor + matrix + gap_open + gap_ext
        elif blast_type == 'tblastn':
            blast_cmd = fixed + evalue + wordsize + processor + matrix + gap_open + gap_ext
        elif blast_type == 'tblastx':
            blast_cmd = fixed + evalue + wordsize + strand + processor + genetic_code + matrix
        
        try:
            os.system(blast_cmd)
            time.sleep(1)
            return 'BLAST done!', temp_blast_file[1]+'.txt'
        except (IOError, OSError):
            logging.debug(time.strftime("%d.%m.%Y um %H:%M:%S Uhr"))
            logging.debug(str(platform.system()+platform.release()))
            logging.exception('BLAST. Got exception on main handler')
            return 'BLAST failed', None
        
class CreateDBThread(QtCore.QThread):
    def __init__(self, database_file_location, db_location, db_name, db_type):
        QtCore.QThread.__init__(self)

        self.database_file_location = database_file_location
        self.db_location = db_location
        self.db_name = db_name
        self.db_type = db_type

    def run(self):
        """Create Bowtie DB."""
        
        if self.db_type == 'BOWTIE':
            info_message, bowtie_path = database_tools.create_bowtie_database(self.db_name,
                                                                                self.database_file_location,
                                                                                self.db_location)
            
            
            self.emit(QtCore.SIGNAL("threadDone(QString, QString)"), info_message[0], 'Bowtie ok')
        
        elif self.db_type == 'BLASTN' or self.db_type == 'BLASTP':
            
            info_message, bowtie_path = database_tools.create_blast_database(self.db_name,
                                                        self.database_file_location,
                                                        self.db_location,
                                                        self.db_type)
            
            
            self.emit(QtCore.SIGNAL("threadDone(QString, QString)"), info_message[0], 'BLAST ok')
        
        elif self.db_type == 'FASTA':
            info_message, bowtie_path = database_tools.create_fasta_db(self.db_name,
                                                        self.database_file_location,
                                                        self.db_location)
            
            
            self.emit(QtCore.SIGNAL("threadDone(QString, QString)"), info_message[0], 'FASTA ok')