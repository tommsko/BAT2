import configparser
import os

core_config: configparser.ConfigParser = configparser.ConfigParser()
core_config.read(os.path.join(os.path.dirname(os.path.dirname(__file__)), "config.ini"))
