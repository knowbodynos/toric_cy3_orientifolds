# toric_cy3_orientifolds

## Requirements
1. Sage with OpenSSL enabled for installing external python packages
    a. For MacOS, you can install it via: https://github.com/3-manifolds/Sage_macOS/releases
2. PyMongo: `sage -pip install pymongo`
3. SymPy: `sage -pip install sympy`
4. python-dotenv: `sage -pip python-dotenv`

## Usage
Set the environment variable `MONGO_URI` to the URI of a local or remote MongoDB database.

### Script
```
scripts/compare_methods.py --query "H11=int(3),VOLFORMPARITY=int(-1)"
```

### Notebook
```
sage -n jupyterlab
```
