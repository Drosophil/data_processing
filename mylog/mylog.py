import logging

logging.basicConfig(
    level=logging.INFO, # root log will show info, warning, error and critical messages
    filename='mainlog.log',
    format='%(asctime)s %(name)s -> %(levelname)s -> %(message)s'
)

my_formatter = logging.Formatter(
    "%(asctime)s %(filename)s, line %(lineno)d => %(levelname)s: %(message)s "
)

process_formatter = logging.Formatter(
    "%(message)s"
)

data_formatter = logging.Formatter(
    "%(message)s"
)


# implement filters
class FilterAllButProcesses(logging.Filter):
    def filter(self, record):
        return (record.msg.startswith("PROCESS:"))

class FilterAllButData(logging.Filter):
    def filter(self, record):
        return (record.msg.startswith("DATA:"))


# implement custom handler
class MyCustomHandler(logging.Handler):
    def __init__(self, filename: str):
        super().__init__()
        self.file = open(filename, "a")

    def emit(self, record: logging.LogRecord):
        self.file.write(self.format(record))  #  .swapcase())
        self.file.write('\n')

    def __del__(self):
        self.file.close()

logger = logging.getLogger()

# default handlers: console, file

console = logging.StreamHandler()
console.setLevel(logging.INFO)  # console will show info, warning, error and critical msgs
console.setFormatter(my_formatter)

process_file = logging.FileHandler(filename="process_log.log")
process_file.setLevel(logging.INFO)  #  additional log file
process_file.setFormatter(process_formatter)

# data_file = MyCustomHandler("data_log.log")
data_file = logging.FileHandler(filename="data_log.log")
data_file.setLevel(logging.INFO) #  custom log file
data_file.setFormatter(data_formatter)

logger.addHandler(console)
logger.addHandler(process_file)
logger.addHandler(data_file)

# connect filters
process_file.addFilter(FilterAllButProcesses())  # filter PROCESS prefixes
data_file.addFilter(FilterAllButData())  # filter DATA prefixes

