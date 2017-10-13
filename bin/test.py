import inspect

def hello(name='World'):
    f = inspect.currentframe().f_back
    mod = f.f_code.co_filename
    lineno = f.f_lineno
    print('Hi, %s. You called this from %s at line # %d.' %
          (name, mod, lineno))

hello()
