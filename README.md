# toric_cy3_orientifolds
Scripts for computing toric Calabi-Yau 3-fold orientifolds

## Requirements
1. Sage with OpenSSL enabled for installing external python packages
    a. For MacOS, you can install it via: https://github.com/3-manifolds/Sage_macOS/releases
2. PyMongo: `sage -pip install pymongo`
3. SymPy: `sage -pip install sympy`
4. python-dotenv: `sage -pip python-dotenv`

## Usage
1. Set the environment variable `MONGO_URI` to the URI of a local or remote MongoDB database.
```
echo "export MONGO_URI=mongodb://username:password@hostname:port/database_name" > ~/.bashrc
```
2. Clone the git repo
```
git clone git@github.com:knowbodynos/toric_cy3_orientifolds.git
cd toric_cy3_orientifolds
```

### Run Script
1. Choose a random sample from the database (based on a query) for comparison:
```
scripts/compare_methods.py --query "H11=int(3),VOLFORMPARITY=int(-1)" sample
```
2. Run comparisons for all database records (based on a query):
```
mkdir -p data
for i in {1..6}; do
    scripts/compare_methods.py -q "H11=int(${i})" full > "data/h11_${i}.json"
done
```

### Run Notebook
```
sage -n jupyterlab
```
