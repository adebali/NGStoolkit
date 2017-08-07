
class voidArgs:
	def __init__(self):
		self.__dict__ = {}

class args:
	def __init__(self, args = voidArgs()):
		self.args = args
		self.dict = vars(args)

	def get(self, key):
		return self.dict.get(key, False)