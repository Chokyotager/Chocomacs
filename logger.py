from datetime import datetime

class Logger ():

    def __init__ (self, log_level, log_file):

        self.log_level = log_level
        self.log_file = open(log_file, "a+")

    def log (self, level=0, message=None):

        if level < self.log_level:
            return None

        # Do stuff
        loggable = str(datetime.now()) + ": " + message

        self.log_file.write(loggable + "\n")
        print(loggable)
