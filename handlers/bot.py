from abc import ABC, abstractclassmethod

class Bot(ABC):

    def __init__(self, official):
        self.official = official

    @abstractclassmethod
    def post_message(self):
        pass

    @abstractclassmethod
    def post_image(self):
        pass
