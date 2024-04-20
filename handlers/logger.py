import os
import pathlib
import logging
from typing import Dict

class Formatter(logging.Formatter):
    """Uses custom colors for traceback to format excepetion"""

    def __init__(self, *args, color=True, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.color = color

    def format_exception(self, ei):
        _, exc_value, traceback = ei
        
class MultiFormatter():
    pass

