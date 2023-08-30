import os
#install packages using pip
dependencies_list = [
    'bootstrapped', 
    'scipy', 
    'matplotlib', 
    'pandas', 
    'numpy', 
    'statsmodels',
    'seaborn',
    'scikit-posthocs',
    'openpyxl']

for package in dependencies_list:
    #check if package is installed
    try:
        __import__(package)
        print(package + ' is installed')
    except ImportError:
        print(package + ' is not installed')
        #if not installed, install
        os.system('pip install ' + package)

